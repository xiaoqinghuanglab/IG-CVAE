#!/usr/bin/env python3
"""Leakage-controlled conditional VAE: MRI MFV128 -> generated gene GFV128."""
from __future__ import annotations
import argparse, hashlib, importlib.util, json, logging, random, sys
from dataclasses import asdict, dataclass
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (accuracy_score, average_precision_score,
    balanced_accuracy_score, f1_score, matthews_corrcoef, precision_score,
    r2_score, recall_score, roc_auc_score)
from sklearn.preprocessing import StandardScaler
from torch.utils.data import DataLoader, TensorDataset

ROOT = Path(__file__).resolve().parents[1]
SPLIT = ROOT / "Output/Intermediate/Splits/model_split_manifest.csv"
RUN_TAG = "dual_cohort_independent_v1"
GENE_TAG = "dual_cohort_top500_v1"
N_FOLDS = 5

def load_local_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod

GENE = load_local_module("xig_gene_branch", ROOT / "scripts/04_run_gene_branch.py")
MRI = load_local_module("xig_mri_branch", ROOT / "scripts/05_run_mri_branch.py")

@dataclass(frozen=True)
class Config:
    gfv_dim: int = 128
    mfv_dim: int = 128
    latent_dim: int = 16
    batch_size: int = 32
    epochs: int = 300
    patience: int = 40
    lr: float = 3e-4
    weight_decay: float = 1e-4
    beta_kl: float = 0.01
    kl_warmup_epochs: int = 40
    seed: int = 20260716

CFG = Config()

def outdir(kind: str, cohort: str) -> Path:
    p = ROOT / "Output" / kind / "CVAE" / RUN_TAG / cohort
    p.mkdir(parents=True, exist_ok=True)
    return p

def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda: f.read(1024 * 1024), b""):
            h.update(b)
    return h.hexdigest()

def save_json(path: Path, obj) -> None:
    def convert(x):
        if isinstance(x, (np.integer,)): return int(x)
        if isinstance(x, (np.floating,)): return float(x)
        if isinstance(x, Path): return str(x)
        raise TypeError(type(x).__name__)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as f:
        json.dump(obj, f, indent=2, sort_keys=True, default=convert)

def set_seed(seed: int) -> None:
    random.seed(seed); np.random.seed(seed); torch.manual_seed(seed)
    if torch.cuda.is_available(): torch.cuda.manual_seed_all(seed)
    torch.backends.cudnn.benchmark = False
    torch.backends.cudnn.deterministic = True

def setup_logger(cohort: str) -> logging.Logger:
    lg = logging.getLogger(f"cvae_{cohort}")
    lg.handlers.clear(); lg.setLevel(logging.INFO)
    fmt = logging.Formatter("%(asctime)s | %(levelname)s | %(message)s")
    for h in (logging.StreamHandler(),
              logging.FileHandler(outdir("Log", cohort) / "development.log")):
        h.setFormatter(fmt); lg.addHandler(h)
    return lg

def cvae_rows(cohort: str) -> pd.DataFrame:
    d = pd.read_csv(SPLIT)
    x = d[(d.cohort == cohort) & (d.model == "cvae") &
          (d.role == "development")].copy()
    if x.empty or x.subject_id.duplicated().any():
        raise RuntimeError(f"{cohort}: invalid CVAE development manifest")
    x["subject_id"] = x.subject_id.astype(str)
    x["diagnosis_binary"] = pd.to_numeric(x.diagnosis_binary).astype(int)
    x["cv_fold"] = pd.to_numeric(x.cv_fold).astype(int)
    if set(x.cv_fold) != set(range(1, N_FOLDS + 1)):
        raise RuntimeError(f"{cohort}: expected fixed folds 1..5")
    return x.sort_values("subject_id").reset_index(drop=True)

class ConditionalVAE(nn.Module):
    """q(z|GFV,MFV), p(z|MFV), and p(GFV|z,MFV)."""
    def __init__(self):
        super().__init__()
        self.prior = nn.Sequential(nn.Linear(CFG.mfv_dim, 128),
            nn.LayerNorm(128), nn.GELU(), nn.Linear(128, 128), nn.GELU())
        self.posterior = nn.Sequential(nn.Linear(CFG.mfv_dim + CFG.gfv_dim, 256),
            nn.LayerNorm(256), nn.GELU(), nn.Linear(256, 128), nn.GELU())
        self.prior_mu = nn.Linear(128, CFG.latent_dim)
        self.prior_logvar = nn.Linear(128, CFG.latent_dim)
        self.post_mu = nn.Linear(128, CFG.latent_dim)
        self.post_logvar = nn.Linear(128, CFG.latent_dim)
        self.decoder = nn.Sequential(nn.Linear(CFG.mfv_dim + CFG.latent_dim, 256),
            nn.LayerNorm(256), nn.GELU(), nn.Linear(256, 128),
            nn.LayerNorm(128), nn.GELU(), nn.Linear(128, CFG.gfv_dim))

    def prior_stats(self, mfv):
        h = self.prior(mfv)
        return self.prior_mu(h), self.prior_logvar(h).clamp(-8, 8)

    def posterior_stats(self, mfv, gfv):
        h = self.posterior(torch.cat([mfv, gfv], dim=1))
        return self.post_mu(h), self.post_logvar(h).clamp(-8, 8)

    def decode(self, mfv, z):
        return self.decoder(torch.cat([mfv, z], dim=1))

    def forward(self, mfv, gfv):
        pm, plv = self.prior_stats(mfv)
        qm, qlv = self.posterior_stats(mfv, gfv)
        z = qm + torch.randn_like(qm) * torch.exp(0.5 * qlv)
        return self.decode(mfv, z), pm, plv, qm, qlv

    def generate(self, mfv):
        pm, _ = self.prior_stats(mfv)
        return self.decode(mfv, pm)

def load_gene_matrix(cohort: str, meta: pd.DataFrame) -> np.ndarray:
    records = [
        {"subject_id": str(r.subject_id),
         "diagnosis_binary": str(int(r.diagnosis_binary)),
         "cv_fold": str(int(r.cv_fold))}
        for r in meta.itertuples(index=False)
    ]
    matrix_path, order_path = GENE.paths_for(cohort)
    X, _, ids, _ = GENE.read_matrix_for_ids(
        matrix_path, GENE.read_gene_order(order_path), records)
    if ids != meta.subject_id.tolist():
        raise RuntimeError(f"{cohort}: gene matrix order mismatch")
    return X.astype(np.float32)

def gene_checkpoint(cohort: str, fold: int) -> tuple[Path, Path]:
    base = ROOT / "Output/Model/Gene_branch" / GENE_TAG / cohort / f"fold_{fold}"
    model_path = base / "seed_1_gene_encoder_classifier.pt"
    scaler_path = base / "standard_scaler.joblib"
    if not model_path.is_file() or not scaler_path.is_file():
        raise FileNotFoundError(f"Missing gene fold-{fold} assets for {cohort}")
    return model_path, scaler_path

def extract_gene_features(cohort: str, fold: int, X: np.ndarray,
                          device: torch.device) -> tuple[np.ndarray, dict]:
    model_path, scaler_path = gene_checkpoint(cohort, fold)
    checkpoint = torch.load(model_path, map_location="cpu", weights_only=False)
    scaler = joblib.load(scaler_path)
    model = GENE.GeneBranch(GENE.TrainConfig(**checkpoint["config"])).to(device)
    model.load_state_dict(checkpoint["state_dict"])
    _, gfv = GENE.predict(model, scaler.transform(X).astype(np.float32), device)
    if gfv.shape != (len(X), CFG.gfv_dim) or not np.isfinite(gfv).all():
        raise RuntimeError(f"{cohort}: invalid fold-{fold} GFV extraction")
    return gfv.astype(np.float32), {
        "gene_model": str(model_path.relative_to(ROOT)),
        "gene_model_sha256": sha256(model_path),
        "gene_scaler": str(scaler_path.relative_to(ROOT)),
        "gene_scaler_sha256": sha256(scaler_path),
        "gene_seed": int(checkpoint["seed"]),
    }

def mri_checkpoint(cohort: str, fold: int) -> Path:
    p = ROOT / "Output/Model/MRI_branch" / RUN_TAG / cohort / f"fold_{fold}_best.pt"
    if not p.is_file():
        raise FileNotFoundError(f"Missing MRI fold-{fold} checkpoint for {cohort}")
    return p

def extract_mri_features(cohort: str, fold: int, meta: pd.DataFrame,
                         device: torch.device) -> tuple[np.ndarray, dict]:
    model_path = mri_checkpoint(cohort, fold)
    checkpoint = torch.load(model_path, map_location="cpu", weights_only=False)
    dropout = float(checkpoint["config"]["dropout"])
    model = MRI.MRIBranch(dropout).to(device)
    model.load_state_dict(checkpoint["state_dict"])
    cache = outdir("Intermediate", cohort) / "mri_preprocessing_cache"
    dl = MRI.loader(meta, cache, False, device)
    _, _, mfv, ids = MRI.predict(model, dl, device)
    by_id = {str(s): v for s, v in zip(ids, mfv)}
    ordered = np.vstack([by_id[s] for s in meta.subject_id]).astype(np.float32)
    if ordered.shape != (len(meta), CFG.mfv_dim) or not np.isfinite(ordered).all():
        raise RuntimeError(f"{cohort}: invalid fold-{fold} MFV extraction")
    return ordered, {
        "mri_model": str(model_path.relative_to(ROOT)),
        "mri_model_sha256": sha256(model_path),
        "mri_best_epoch": int(checkpoint["best_epoch"]),
        "mri_target_shape": list(checkpoint["target_shape"]),
    }

def fold_features(cohort: str, fold: int, meta: pd.DataFrame,
                  gene_matrix: np.ndarray, device: torch.device):
    """Features share one frozen encoder coordinate system within this CVAE fold."""
    gfv, gene_info = extract_gene_features(cohort, fold, gene_matrix, device)
    mfv, mri_info = extract_mri_features(cohort, fold, meta, device)
    return mfv, gfv, {
        "cohort": cohort, "cvae_fold": fold,
        "representation_rule": (
            "All paired CVAE development subjects were embedded with the "
            "matched fold-specific frozen gene and MRI encoders. Fixed-fold "
            "agreement verifies the CVAE validation subjects were excluded "
            "from those branch encoders' training sets."),
        **gene_info, **mri_info,
    }

def cvae_loss(generated, target, pm, plv, qm, qlv, beta):
    reconstruction = F.mse_loss(generated, target)
    kl = 0.5 * torch.mean(torch.sum(
        plv - qlv + (qlv.exp() + (qm - pm).pow(2)) / plv.exp() - 1.0, dim=1))
    return reconstruction + beta * kl, reconstruction, kl

def reconstruction_metrics(target, generated):
    cosine = np.sum(target * generated, axis=1) / np.clip(
        np.linalg.norm(target, axis=1) * np.linalg.norm(generated, axis=1), 1e-12, None)
    return {
        "mse": float(np.mean((target - generated) ** 2)),
        "mae": float(np.mean(np.abs(target - generated))),
        "cosine_similarity": float(np.mean(cosine)),
        "r2": float(r2_score(target.reshape(-1), generated.reshape(-1))),
    }

def predict_generated(model, mfv, device):
    model.eval(); out = []
    with torch.no_grad():
        for start in range(0, len(mfv), CFG.batch_size):
            x = torch.from_numpy(mfv[start:start + CFG.batch_size]).to(device)
            out.append(model.generate(x).cpu().numpy())
    return np.concatenate(out).astype(np.float32)

def class_metrics(y, probability):
    pred = (probability >= 0.5).astype(int)
    return {
        "roc_auc": float(roc_auc_score(y, probability)),
        "pr_auc": float(average_precision_score(y, probability)),
        "accuracy": float(accuracy_score(y, pred)),
        "balanced_accuracy": float(balanced_accuracy_score(y, pred)),
        "sensitivity": float(recall_score(y, pred, zero_division=0)),
        "precision": float(precision_score(y, pred, zero_division=0)),
        "f1": float(f1_score(y, pred, zero_division=0)),
        "mcc": float(matthews_corrcoef(y, pred)),
    }

def verify_fold_agreement(cohort: str, meta: pd.DataFrame) -> None:
    d = pd.read_csv(SPLIT)
    z = meta[["subject_id", "cv_fold"]].copy()
    for model in ["gene_branch", "mri_branch"]:
        q = d[(d.cohort == cohort) & (d.model == model) &
              (d.role == "development")][["subject_id", "cv_fold"]].copy()
        q["subject_id"] = q.subject_id.astype(str)
        col = f"{model}_fold"
        q[col] = pd.to_numeric(q["cv_fold"]).astype(int)
        q = q[["subject_id", col]]
        z = z.merge(q, on="subject_id", how="left", validate="one_to_one")
        if z[col].isna().any() or not (z["cv_fold"] == z[col]).all():
            raise RuntimeError(f"{cohort}: CVAE/{model} fold agreement failed")

def fit_fold(mfv_train, gfv_train, mfv_val, gfv_val, device, seed):
    set_seed(seed)
    model = ConditionalVAE().to(device)
    opt = torch.optim.AdamW(model.parameters(), lr=CFG.lr,
                            weight_decay=CFG.weight_decay)
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
        opt, mode="min", factor=0.5, patience=12, min_lr=1e-5)
    loader = DataLoader(
        TensorDataset(torch.from_numpy(mfv_train), torch.from_numpy(gfv_train)),
        batch_size=CFG.batch_size, shuffle=True)
    best_mse, best_epoch, stale, best_state = np.inf, 0, 0, None

    for epoch in range(1, CFG.epochs + 1):
        model.train()
        beta = CFG.beta_kl * min(1.0, epoch / CFG.kl_warmup_epochs)
        for mfv, gfv in loader:
            mfv, gfv = mfv.to(device), gfv.to(device)
            opt.zero_grad(set_to_none=True)
            generated, pm, plv, qm, qlv = model(mfv, gfv)
            loss, _, _ = cvae_loss(generated, gfv, pm, plv, qm, qlv, beta)
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 5.0)
            opt.step()

        generated = predict_generated(model, mfv_val, device)
        val_mse = float(np.mean((gfv_val - generated) ** 2))
        scheduler.step(val_mse)
        if val_mse < best_mse - 1e-7:
            best_mse, best_epoch, stale = val_mse, epoch, 0
            best_state = {k: v.detach().cpu().clone()
                          for k, v in model.state_dict().items()}
        else:
            stale += 1
            if stale >= CFG.patience:
                break

    model.load_state_dict(best_state)
    return model, best_epoch

def development(cohort: str, device: torch.device) -> None:
    log = setup_logger(cohort)
    meta = cvae_rows(cohort)
    verify_fold_agreement(cohort, meta)
    gene_matrix = load_gene_matrix(cohort, meta)
    folds = meta.cv_fold.to_numpy()

    log.info("%s CVAE development-only run started on %s.", cohort, device)
    log.info("Locked-test IDs, features, labels, predictions, and metrics are not loaded.")
    model_dir, int_dir, table_dir, qc_dir = (
        outdir("Model", cohort), outdir("Intermediate", cohort),
        outdir("Table", cohort), outdir("QC", cohort))
    oof_rows, fold_rows, epochs = [], [], []

    for fold in range(1, N_FOLDS + 1):
        tr, va = folds != fold, folds == fold
        mfv_raw, gfv_raw, branch_info = fold_features(
            cohort, fold, meta, gene_matrix, device)

        mfv_scaler, gfv_scaler = StandardScaler(), StandardScaler()
        mfv_tr = mfv_scaler.fit_transform(mfv_raw[tr]).astype(np.float32)
        gfv_tr = gfv_scaler.fit_transform(gfv_raw[tr]).astype(np.float32)
        mfv_va = mfv_scaler.transform(mfv_raw[va]).astype(np.float32)
        gfv_va = gfv_scaler.transform(gfv_raw[va]).astype(np.float32)

        model, best_epoch = fit_fold(
            mfv_tr, gfv_tr, mfv_va, gfv_va, device, CFG.seed + fold)
        generated_scaled = predict_generated(model, mfv_va, device)
        generated = gfv_scaler.inverse_transform(generated_scaled)
        recon = reconstruction_metrics(gfv_raw[va], generated)

        classifier = LogisticRegression(
            C=1.0, class_weight="balanced", max_iter=2000,
            random_state=CFG.seed + fold)
        classifier.fit(gfv_tr, meta.loc[tr, "diagnosis_binary"])
        probability = classifier.predict_proba(generated_scaled)[:, 1]
        clf = class_metrics(meta.loc[va, "diagnosis_binary"].to_numpy(), probability)

        rows = meta.loc[va, ["subject_id","diagnosis_binary","cv_fold"]].reset_index(drop=True)
        rows["probability_AD_from_generated_GFV"] = probability
        agfv = pd.DataFrame(generated, columns=[
            f"AGFV{j+1:03d}" for j in range(CFG.gfv_dim)])
        oof_rows.append(pd.concat([rows, agfv], axis=1))

        record = {"fold": fold, "n_train": int(tr.sum()),
                  "n_validation": int(va.sum()), "best_epoch": best_epoch,
                  **recon, **clf}
        fold_rows.append(record); epochs.append(best_epoch)

        torch.save({
            "state_dict": model.state_dict(), "best_epoch": best_epoch,
            "mfv_mean": mfv_scaler.mean_, "mfv_scale": mfv_scaler.scale_,
            "gfv_mean": gfv_scaler.mean_, "gfv_scale": gfv_scaler.scale_,
            "logistic_coef": classifier.coef_,
            "logistic_intercept": classifier.intercept_,
            "config": asdict(CFG), "branch_provenance": branch_info,
        }, model_dir / f"fold_{fold}_best.pt")
        save_json(qc_dir / f"fold_{fold}_branch_provenance.json", branch_info)
        log.info("fold=%d epoch=%d MSE=%.4f R2=%.3f AUC=%.3f",
                 fold, best_epoch, recon["mse"], recon["r2"], clf["roc_auc"])

    oof = pd.concat(oof_rows, ignore_index=True).sort_values("subject_id")
    if len(oof) != len(meta) or oof.subject_id.duplicated().any():
        raise RuntimeError(f"{cohort}: CVAE OOF output integrity failure")
    oof.to_csv(int_dir / "development_OOF_generated_GFV128.csv", index=False)

    metrics = pd.DataFrame(fold_rows)
    metrics.to_csv(table_dir / "fold_metrics.csv", index=False)
    numeric = [c for c in metrics.columns if c not in
               {"fold","n_train","n_validation","best_epoch"}]
    summary = {f"{c}_mean": float(metrics[c].mean()) for c in numeric}
    summary.update({f"{c}_sd": float(metrics[c].std(ddof=1)) for c in numeric})
    save_json(table_dir / "fivefold_mean_sd.json", summary)
    save_json(qc_dir / "development_provenance.json", {
        "cohort": cohort, "n_development": len(meta),
        "fold_counts": {str(k): int(v) for k, v in
                        meta.cv_fold.value_counts().sort_index().items()},
        "best_epochs": epochs, "config": asdict(CFG),
        "locked_test_accessed": False,
        "representation_strategy": "matched_frozen_branch_models_per_CVAE_fold",
    })
    log.info("Development complete; locked test remains unopened.")

def locked_rows(cohort: str) -> pd.DataFrame:
    d = pd.read_csv(SPLIT)
    x = d[(d.cohort == cohort) & (d.model == "cvae") &
          (d.role == "locked_test")].copy()
    if x.empty or x.subject_id.duplicated().any():
        raise RuntimeError(f"{cohort}: invalid CVAE locked-test manifest")
    x["subject_id"] = x.subject_id.astype(str)
    x["diagnosis_binary"] = pd.to_numeric(x.diagnosis_binary).astype(int)
    x["cv_fold"] = -1
    return x.sort_values("subject_id").reset_index(drop=True)

def extract_final_gene_features(cohort, X, device):
    base = ROOT / "Output/Model/Gene_branch" / GENE_TAG / cohort / "final_development_models"
    model_path = base / "seed_1_gene_encoder_classifier.pt"
    scaler_path = base / "standard_scaler.joblib"
    checkpoint = torch.load(model_path, map_location="cpu", weights_only=False)
    scaler = joblib.load(scaler_path)
    model = GENE.GeneBranch(GENE.TrainConfig(**checkpoint["config"])).to(device)
    model.load_state_dict(checkpoint["state_dict"])
    _, gfv = GENE.predict(model, scaler.transform(X).astype(np.float32), device)
    return gfv.astype(np.float32), {
        "gene_model": str(model_path.relative_to(ROOT)),
        "gene_model_sha256": sha256(model_path),
        "gene_scaler": str(scaler_path.relative_to(ROOT)),
        "gene_scaler_sha256": sha256(scaler_path),
        "gene_seed": int(checkpoint["seed"]),
    }

def extract_final_mri_features(cohort, meta, device):
    model_path = ROOT / "Output/Model/MRI_branch" / RUN_TAG / cohort / "final_development_fit.pt"
    checkpoint = torch.load(model_path, map_location="cpu", weights_only=False)
    model = MRI.MRIBranch(float(checkpoint["config"]["dropout"])).to(device)
    model.load_state_dict(checkpoint["state_dict"])
    cache = outdir("Intermediate", cohort) / "mri_preprocessing_cache"
    _, _, mfv, ids = MRI.predict(model, MRI.loader(meta, cache, False, device), device)
    by_id = {str(s): v for s, v in zip(ids, mfv)}
    ordered = np.vstack([by_id[s] for s in meta.subject_id]).astype(np.float32)
    return ordered, {
        "mri_model": str(model_path.relative_to(ROOT)),
        "mri_model_sha256": sha256(model_path),
        "mri_epochs": int(checkpoint.get("epochs", checkpoint.get("best_epoch", -1))),
        "mri_target_shape": list(checkpoint["target_shape"]),
    }

def final_features(cohort, meta, device):
    X = load_gene_matrix(cohort, meta)
    gfv, gene_info = extract_final_gene_features(cohort, X, device)
    mfv, mri_info = extract_final_mri_features(cohort, meta, device)
    if not np.isfinite(gfv).all() or not np.isfinite(mfv).all():
        raise RuntimeError(f"{cohort}: non-finite final branch features")
    return mfv, gfv, {**gene_info, **mri_info}

def fit_final(mfv, gfv, device, epochs, seed):
    set_seed(seed)
    model = ConditionalVAE().to(device)
    opt = torch.optim.AdamW(model.parameters(), lr=CFG.lr, weight_decay=CFG.weight_decay)
    loader = DataLoader(TensorDataset(torch.from_numpy(mfv), torch.from_numpy(gfv)),
                        batch_size=CFG.batch_size, shuffle=True)
    for epoch in range(1, epochs + 1):
        model.train()
        beta = CFG.beta_kl * min(1.0, epoch / CFG.kl_warmup_epochs)
        for mb, gb in loader:
            mb, gb = mb.to(device), gb.to(device)
            opt.zero_grad(set_to_none=True)
            generated, pm, plv, qm, qlv = model(mb, gb)
            loss, _, _ = cvae_loss(generated, gb, pm, plv, qm, qlv, beta)
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 5.0)
            opt.step()
    return model

def test_metrics(y, probability):
    out = class_metrics(y, probability)
    out["brier"] = float(np.mean((y - probability) ** 2))
    return out

def stratified_bootstrap_ci(y, probability, observed, generated, n_bootstrap):
    rng = np.random.default_rng(CFG.seed + 99)
    groups = [np.flatnonzero(y == c) for c in [0, 1]]
    names = list(test_metrics(y, probability)) + list(reconstruction_metrics(observed, generated))
    values = {name: [] for name in names}
    for _ in range(n_bootstrap):
        idx = np.concatenate([rng.choice(g, len(g), replace=True) for g in groups])
        for name, value in {
            **test_metrics(y[idx], probability[idx]),
            **reconstruction_metrics(observed[idx], generated[idx]),
        }.items():
            values[name].append(value)
    return {name: [float(np.quantile(v, .025)), float(np.quantile(v, .975))]
            for name, v in values.items()}

def run_locked_test(cohort, device, allow_locked_test, bootstrap_replicates):
    if not allow_locked_test:
        raise PermissionError("TEST requires --allow-locked-test")
    marker = outdir("QC", cohort) / "locked_test_opened_once.json"
    if marker.exists():
        raise RuntimeError(f"{cohort}: locked test already opened: {marker}")

    log = setup_logger(cohort)
    dev_meta = cvae_rows(cohort)
    dev_mfv_raw, dev_gfv_raw, branch_info = final_features(cohort, dev_meta, device)
    mfv_scaler, gfv_scaler = StandardScaler(), StandardScaler()
    dev_mfv = mfv_scaler.fit_transform(dev_mfv_raw).astype(np.float32)
    dev_gfv = gfv_scaler.fit_transform(dev_gfv_raw).astype(np.float32)

    fold_metrics = pd.read_csv(outdir("Table", cohort) / "fold_metrics.csv")
    final_epochs = max(1, int(round(float(fold_metrics["best_epoch"].mean()))))
    model = fit_final(dev_mfv, dev_gfv, device, final_epochs, CFG.seed + 1000)
    classifier = LogisticRegression(
        C=1.0, class_weight="balanced", max_iter=2000, random_state=CFG.seed + 1000)
    classifier.fit(dev_gfv, dev_meta["diagnosis_binary"])

    save_json(marker, {
        "cohort": cohort, "status": "opened_once_before_locked_data_load",
        "final_epochs_from_mean_CV_best_epoch": final_epochs,
    })
    log.warning("Opening %s CVAE locked test exactly once.", cohort)

    test_meta = locked_rows(cohort)
    test_mfv_raw, test_gfv_raw, test_branch_info = final_features(cohort, test_meta, device)
    generated_scaled = predict_generated(
        model, mfv_scaler.transform(test_mfv_raw).astype(np.float32), device)
    generated = gfv_scaler.inverse_transform(generated_scaled)
    probability = classifier.predict_proba(generated_scaled)[:, 1]
    y = test_meta["diagnosis_binary"].to_numpy()

    reconstruction = reconstruction_metrics(test_gfv_raw, generated)
    classification = test_metrics(y, probability)
    ci = stratified_bootstrap_ci(
        y, probability, test_gfv_raw, generated, bootstrap_replicates)

    rows = test_meta[["subject_id","diagnosis_binary"]].copy()
    rows["probability_AD_from_generated_GFV"] = probability
    rows = pd.concat([rows.reset_index(drop=True),
        pd.DataFrame(test_gfv_raw, columns=[f"GFV{i+1:03d}" for i in range(CFG.gfv_dim)]),
        pd.DataFrame(generated, columns=[f"AGFV{i+1:03d}" for i in range(CFG.gfv_dim)])], axis=1)

    int_dir, table_dir, model_dir = (
        outdir("Intermediate", cohort), outdir("Table", cohort), outdir("Model", cohort))
    rows.to_csv(int_dir / "locked_test_generated_GFV128.csv", index=False)
    result = {
        "cohort": cohort, "status": "locked_test_complete",
        "n": int(len(y)), "AD": int(y.sum()), "CN": int((y == 0).sum()),
        "bootstrap_replicates": bootstrap_replicates,
        "final_epochs_from_mean_CV_best_epoch": final_epochs,
        "reconstruction": reconstruction, "classification": classification,
        "stratified_bootstrap_95ci": ci,
        "branch_provenance_development": branch_info,
        "branch_provenance_locked_test": test_branch_info,
    }
    save_json(table_dir / "locked_test_metrics.json", result)
    torch.save({
        "state_dict": model.state_dict(), "final_epochs": final_epochs,
        "mfv_mean": mfv_scaler.mean_, "mfv_scale": mfv_scaler.scale_,
        "gfv_mean": gfv_scaler.mean_, "gfv_scale": gfv_scaler.scale_,
        "logistic_coef": classifier.coef_, "logistic_intercept": classifier.intercept_,
        "config": asdict(CFG), "branch_provenance": branch_info,
    }, model_dir / "final_development_fit.pt")
    log.info("Locked test complete: %s", json.dumps(result, sort_keys=True))

def validate(cohort: str) -> None:
    meta = cvae_rows(cohort)
    verify_fold_agreement(cohort, meta)
    gene_matrix = load_gene_matrix(cohort, meta)
    report = {
        "cohort": cohort,
        "development_n": int(len(meta)),
        "AD": int(meta.diagnosis_binary.sum()),
        "CN": int((meta.diagnosis_binary == 0).sum()),
        "fold_counts": {str(k): int(v) for k, v in
                        meta.cv_fold.value_counts().sort_index().items()},
        "gene_matrix_shape": list(gene_matrix.shape),
        "fold_agreement_verified": True,
        "locked_test_accessed": False,
    }
    print(json.dumps(report, indent=2, sort_keys=True))

def main() -> None:
    ap = argparse.ArgumentParser(
        description="Leakage-controlled CVAE using fold-matched frozen branch encoders.")
    ap.add_argument("command", choices=["VALIDATE", "DEVELOP", "TEST"])
    ap.add_argument("--cohort", choices=["ADNI", "AddNeuroMed", "ALL"], default="ALL")
    ap.add_argument("--device", choices=["auto", "cpu", "cuda"], default="auto")
    ap.add_argument("--allow-locked-test", action="store_true")
    ap.add_argument("--bootstrap-replicates", type=int, default=2000)
    args = ap.parse_args()

    if args.device == "auto":
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    else:
        device = torch.device(args.device)
    if device.type == "cuda" and not torch.cuda.is_available():
        raise RuntimeError("CUDA requested but unavailable.")

    cohorts = ["ADNI", "AddNeuroMed"] if args.cohort == "ALL" else [args.cohort]
    for cohort in cohorts:
        if args.command == "VALIDATE":
            validate(cohort)
        elif args.command == "DEVELOP":
            development(cohort, device)
        else:
            run_locked_test(
                cohort, device, args.allow_locked_test, args.bootstrap_replicates)

if __name__ == "__main__":
    main()
