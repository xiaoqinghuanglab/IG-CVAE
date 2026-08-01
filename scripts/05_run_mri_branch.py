#!/usr/bin/env python3
"""Leakage-controlled dual-cohort structural-MRI branch: MRI -> MFV128 -> AD/CN."""
from __future__ import annotations
import argparse, gzip, hashlib, json, logging, random, re, struct
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F
from sklearn.metrics import (accuracy_score, average_precision_score, brier_score_loss,
                             confusion_matrix, f1_score, matthews_corrcoef,
                             precision_recall_curve, precision_score, roc_auc_score,
                             roc_curve)
from torch.utils.data import DataLoader, Dataset
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
SPLIT = ROOT / "Output/Intermediate/Splits/model_split_manifest.csv"
RUN_TAG = "dual_cohort_independent_v1"
TARGET_SHAPE = (48, 56, 48)
N_FOLDS = 5

@dataclass(frozen=True)
class Config:
    batch_size: int = 8
    epochs: int = 60
    patience: int = 15
    lr: float = 1e-4
    weight_decay: float = 1e-4
    dropout: float = 0.00
    num_workers: int = 4
    seed: int = 20260716

CFG = Config()

def outdir(kind: str, cohort: str) -> Path:
    p = ROOT / "Output" / kind / "MRI_branch" / RUN_TAG / cohort
    p.mkdir(parents=True, exist_ok=True)
    return p

def save_json(path: Path, obj: Any) -> None:
    def convert(x):
        if isinstance(x, np.integer): return int(x)
        if isinstance(x, np.floating): return float(x)
        if isinstance(x, Path): return str(x)
        raise TypeError(type(x).__name__)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as f:
        json.dump(obj, f, indent=2, sort_keys=True, default=convert)

def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda: f.read(1024 * 1024), b""): h.update(b)
    return h.hexdigest()

def set_seed(seed: int) -> None:
    random.seed(seed); np.random.seed(seed); torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    torch.backends.cudnn.benchmark = False
    torch.backends.cudnn.deterministic = True

def setup_logger(cohort: str) -> logging.Logger:
    lg = logging.getLogger(f"mri_{cohort}")
    lg.handlers.clear(); lg.setLevel(logging.INFO)
    fmt = logging.Formatter("%(asctime)s | %(levelname)s | %(message)s")
    for h in (logging.StreamHandler(),
              logging.FileHandler(outdir("Log", cohort) / "development.log")):
        h.setFormatter(fmt); lg.addHandler(h)
    return lg

def load_rows(cohort: str, role: str) -> pd.DataFrame:
    df = pd.read_csv(SPLIT)
    x = df[(df["cohort"] == cohort) & (df["model"] == "mri_branch") &
           (df["role"] == role)].copy()
    if x.empty: raise RuntimeError(f"No {cohort} MRI rows for role={role}")
    x["diagnosis_binary"] = pd.to_numeric(x["diagnosis_binary"], errors="raise").astype(int)
    if not set(x["diagnosis_binary"]).issubset({0, 1}):
        raise RuntimeError("diagnosis_binary must be CN=0 or AD=1")
    if x["mri_path"].isna().any() or not x["mri_path"].map(lambda p: Path(p).is_file()).all():
        raise RuntimeError(f"{cohort} has a missing MRI path")
    if x["subject_id"].duplicated().any():
        raise RuntimeError(f"{cohort} {role} contains duplicate subject IDs")
    return x.sort_values("subject_id").reset_index(drop=True)

def read_nifti(path: str) -> np.ndarray:
    with gzip.open(path, "rb") as f:
        hdr = f.read(348)
        if len(hdr) != 348: raise RuntimeError(f"Short NIfTI header: {path}")
        endian = "<" if struct.unpack("<i", hdr[:4])[0] == 348 else ">"
        if struct.unpack(endian + "i", hdr[:4])[0] != 348:
            raise RuntimeError(f"Invalid NIfTI header: {path}")
        dim = struct.unpack(endian + "8h", hdr[40:56])
        datatype = struct.unpack(endian + "h", hdr[70:72])[0]
        vox_offset = struct.unpack(endian + "f", hdr[108:112])[0]
        slope = struct.unpack(endian + "f", hdr[112:116])[0]
        intercept = struct.unpack(endian + "f", hdr[116:120])[0]
        shape = tuple(int(v) for v in dim[1:4])
        dtypes = {2:"u1", 4:"i2", 8:"i4", 16:"f4", 64:"f8", 512:"u2", 768:"u4"}
        if datatype not in dtypes: raise RuntimeError(f"Unsupported datatype {datatype}: {path}")
        dt = np.dtype(endian + dtypes[datatype])
        f.seek(max(352, int(round(vox_offset))))
        raw = f.read(int(np.prod(shape)) * dt.itemsize)
    arr = np.frombuffer(raw, dtype=dt, count=int(np.prod(shape)))
    if arr.size != int(np.prod(shape)): raise RuntimeError(f"Short NIfTI image: {path}")
    arr = arr.reshape(shape, order="F").astype(np.float32, copy=False)
    if np.isfinite(slope) and slope not in (0.0, 1.0): arr *= slope
    if np.isfinite(intercept) and intercept != 0.0: arr += intercept
    return arr

def preprocess_volume(path: str) -> np.ndarray:
    x = read_nifti(path)
    mask = np.isfinite(x) & (x != 0)
    if mask.sum() < 1000: raise RuntimeError(f"Insufficient brain voxels: {path}")
    vals = x[mask]; sd = float(vals.std())
    if not np.isfinite(sd) or sd <= 1e-8: raise RuntimeError(f"Invalid intensity SD: {path}")
    x = np.where(mask, np.clip((x - float(vals.mean())) / sd, -5.0, 5.0), 0.0).astype(np.float32)
    return F.interpolate(torch.from_numpy(x)[None,None], size=TARGET_SHAPE,
                         mode="trilinear", align_corners=False)[0,0].numpy().astype(np.float32)

def cache_path(cache_dir: Path, row: pd.Series) -> Path:
    safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(row["subject_id"]))
    return cache_dir / f"{safe}.npy"

class MRIDataset(Dataset):
    def __init__(self, rows: pd.DataFrame, cache_dir: Path):
        self.rows = rows.reset_index(drop=True)
        self.cache_dir = cache_dir
        self.cache_dir.mkdir(parents=True, exist_ok=True)
    def __len__(self): return len(self.rows)
    def __getitem__(self, idx: int):
        r = self.rows.iloc[idx]; p = cache_path(self.cache_dir, r)
        if p.exists():
            x = np.load(p, allow_pickle=False)
        else:
            x = preprocess_volume(str(r["mri_path"]))
            tmp = p.with_suffix(".tmp.npy")
            np.save(tmp, x); tmp.replace(p)
        return (torch.from_numpy(np.asarray(x, dtype=np.float32))[None],
                torch.tensor(int(r["diagnosis_binary"]), dtype=torch.float32),
                str(r["subject_id"]))

class MRIBranch(nn.Module):
    """Lightweight 8-convolution MRI CNN; output is MFV128."""
    def __init__(self, dropout: float):
        super().__init__()
        def conv(a, b):
            return nn.Sequential(
                nn.Conv3d(a, b, 3, padding=1, bias=False),
                nn.GroupNorm(8, b),
                nn.GELU(),
            )
        def stage(a, b):
            return nn.Sequential(conv(a, b), conv(b, b), nn.MaxPool3d(2))
        self.features = nn.Sequential(
            stage(1, 16), stage(16, 32), stage(32, 64), stage(64, 64),
            nn.AdaptiveAvgPool3d(1),
        )
        self.mfv = nn.Sequential(
            nn.Flatten(), nn.Linear(64, 128), nn.LayerNorm(128), nn.GELU(),
            nn.Dropout(dropout),
        )
        self.head = nn.Sequential(
            nn.Linear(128, 64), nn.LayerNorm(64), nn.GELU(),
            nn.Dropout(dropout), nn.Linear(64, 1),
        )
    def forward(self, x):
        z = self.mfv(self.features(x))
        return self.head(z).squeeze(1), z

def loader(rows, cache_dir, shuffle, device):
    return DataLoader(MRIDataset(rows, cache_dir), batch_size=CFG.batch_size,
        shuffle=shuffle, num_workers=CFG.num_workers, pin_memory=device.type=="cuda",
        persistent_workers=CFG.num_workers > 0)

def predict(model, dl, device):
    model.eval(); probs=[]; mfvs=[]; ids=[]; ys=[]
    with torch.no_grad():
        for x, y, sid in dl:
            x = x.to(device, non_blocking=True)
            logits, z = model(x)
            probs.extend(torch.sigmoid(logits).cpu().numpy().tolist())
            mfvs.append(z.cpu().numpy()); ids.extend(sid); ys.extend(y.numpy().astype(int).tolist())
    return np.asarray(ys), np.asarray(probs), np.concatenate(mfvs), ids

def ece(y, p, bins=10):
    edges = np.linspace(0, 1, bins + 1); total = len(y); ans = 0.0
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (p >= lo) & ((p < hi) if hi < 1 else (p <= hi))
        if m.any(): ans += m.mean() * abs(float(y[m].mean()) - float(p[m].mean()))
    return float(ans)

def metrics(y, p, threshold=0.5):
    q = (p >= threshold).astype(int)
    tn, fp, fn, tp = confusion_matrix(y, q, labels=[0,1]).ravel()
    return {"roc_auc":float(roc_auc_score(y,p)), "pr_auc":float(average_precision_score(y,p)),
      "accuracy":float(accuracy_score(y,q)), "balanced_accuracy":float((tp/(tp+fn)+tn/(tn+fp))/2),
      "sensitivity":float(tp/(tp+fn)), "specificity":float(tn/(tn+fp)),
      "precision":float(precision_score(y,q,zero_division=0)), "f1":float(f1_score(y,q,zero_division=0)),
      "mcc":float(matthews_corrcoef(y,q)), "brier":float(brier_score_loss(y,p)),
      "ece":ece(y,p), "tn":int(tn), "fp":int(fp), "fn":int(fn), "tp":int(tp),
      "threshold":float(threshold)}

def choose_oof_threshold(y, p):
    candidates = np.unique(np.r_[0.01, np.linspace(0.02,0.98,97), p])
    bas = [metrics(y,p,float(t))["balanced_accuracy"] for t in candidates]
    best = np.flatnonzero(np.isclose(bas, np.max(bas)))
    return float(candidates[best[np.argmin(abs(candidates[best]-0.5))]])

def augment_batch(x: torch.Tensor) -> torch.Tensor:
    """Training-only, separate affine transforms; never used for validation/test."""
    n, _, d, h, w = x.shape
    theta = torch.eye(3, 4, device=x.device, dtype=x.dtype).unsqueeze(0).repeat(n, 1, 1)
    operation = torch.randint(0, 4, (n,), device=x.device)  # identity, zoom, shift, rotation

    zoom = operation == 1
    if zoom.any():
        scale = 0.90 + 0.20 * torch.rand(int(zoom.sum()), device=x.device, dtype=x.dtype)
        theta[zoom, 0, 0] = 1.0 / scale
        theta[zoom, 1, 1] = 1.0 / scale
        theta[zoom, 2, 2] = 1.0 / scale

    shift = operation == 2
    if shift.any():
        theta[shift, 0, 3] = torch.empty(int(shift.sum()), device=x.device).uniform_(-0.08, 0.08)
        theta[shift, 1, 3] = torch.empty(int(shift.sum()), device=x.device).uniform_(-0.08, 0.08)
        theta[shift, 2, 3] = torch.empty(int(shift.sum()), device=x.device).uniform_(-0.08, 0.08)

    rotate = operation == 3
    if rotate.any():
        angle = torch.empty(int(rotate.sum()), device=x.device).uniform_(-10.0, 10.0) * np.pi / 180.0
        c, sn = torch.cos(angle), torch.sin(angle)
        theta[rotate, 0, 0] = c;  theta[rotate, 0, 1] = -sn
        theta[rotate, 1, 0] = sn; theta[rotate, 1, 1] = c

    grid = F.affine_grid(theta, x.size(), align_corners=False)
    return F.grid_sample(x, grid, mode="bilinear", padding_mode="zeros", align_corners=False)

def train_fold(train_rows, val_rows, cache_dir, device, seed, log):
    set_seed(seed)
    tr = loader(train_rows, cache_dir, True, device)
    va = loader(val_rows, cache_dir, False, device)
    model = MRIBranch(CFG.dropout).to(device)
    n0 = int((train_rows.diagnosis_binary == 0).sum())
    n1 = int((train_rows.diagnosis_binary == 1).sum())
    loss_fn = nn.BCEWithLogitsLoss(pos_weight=torch.tensor([n0/n1], device=device))
    opt = torch.optim.AdamW(model.parameters(), lr=CFG.lr, weight_decay=CFG.weight_decay)
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
        opt, mode="max", factor=0.5, patience=4, min_lr=1e-6)
    amp = device.type == "cuda"
    scaler = torch.amp.GradScaler("cuda", enabled=amp)
    best_auc, best_epoch, stale, best_state = -np.inf, 0, 0, None
    for epoch in range(1, CFG.epochs + 1):
        model.train()
        for x, y, _ in tr:
            x=x.to(device,non_blocking=True); y=y.to(device,non_blocking=True)
            x=augment_batch(x)
            opt.zero_grad(set_to_none=True)
            with torch.autocast(device_type=device.type, dtype=torch.float16, enabled=amp):
                logits,_ = model(x); loss = loss_fn(logits,y)
            scaler.scale(loss).backward(); scaler.step(opt); scaler.update()
        vy, vp, _, _ = predict(model, va, device)
        auc = roc_auc_score(vy, vp)
        scheduler.step(auc)
        if auc > best_auc + 1e-7:
            best_auc, best_epoch, stale = auc, epoch, 0
            best_state = {k:v.detach().cpu().clone() for k,v in model.state_dict().items()}
        else:
            stale += 1
            if stale >= CFG.patience: break
    model.load_state_dict(best_state)
    vy, vp, vz, vid = predict(model, va, device)
    log.info("fold validation: best_epoch=%d AUC=%.3f AP=%.3f BA=%.3f",
             best_epoch, roc_auc_score(vy,vp), average_precision_score(vy,vp),
             metrics(vy,vp)["balanced_accuracy"])
    return model, vy, vp, vz, vid, best_epoch

def plot_curves(y, p, cohort, phase, png, pdf):
    fpr,tpr,_ = roc_curve(y,p); prec,rec,_ = precision_recall_curve(y,p)
    fig, ax = plt.subplots(1,2,figsize=(8.2,3.5))
    ax[0].plot(fpr,tpr,lw=2.2,color="#1F5A99",
               label=f"AUC = {roc_auc_score(y,p):.3f}")
    ax[0].plot([0,1],[0,1],"--",color="0.55",lw=1)
    ax[0].set(xlabel="False-positive rate",ylabel="True-positive rate",
              title=f"{cohort}: MRI ROC ({phase})",xlim=(0,1),ylim=(0,1))
    ax[1].plot(rec,prec,lw=2.2,color="#B23A48",
               label=f"AP = {average_precision_score(y,p):.3f}")
    ax[1].axhline(y.mean(),ls="--",color="0.55",lw=1,label=f"Prevalence = {y.mean():.3f}")
    ax[1].set(xlabel="Recall",ylabel="Precision",title=f"{cohort}: MRI PR ({phase})",
              xlim=(0,1),ylim=(0,1))
    for a in ax:
        a.legend(frameon=False,loc="lower right"); a.spines[["top","right"]].set_visible(False)
    fig.tight_layout(); png.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(png,dpi=600,bbox_inches="tight"); fig.savefig(pdf,bbox_inches="tight")
    plt.close(fig)

def development(cohort, device):
    log = setup_logger(cohort)
    rows = load_rows(cohort, "development")
    if rows.cv_fold.isna().any() or set(rows.cv_fold.astype(int)) != set(range(1,N_FOLDS+1)):
        raise RuntimeError(f"{cohort}: invalid fixed development folds")
    log.info("%s MRI development-only run started on %s; locked test is not loaded.",cohort,device)
    cache = outdir("Intermediate",cohort) / "development_cache"
    model_dir, int_dir, table_dir, qc_dir = (outdir("Model",cohort), outdir("Intermediate",cohort),
                                               outdir("Table",cohort), outdir("QC",cohort))
    oof=[]; fold_rows=[]; best_epochs=[]
    for fold in range(1,N_FOLDS+1):
        tr=rows[rows.cv_fold.astype(int)!=fold]; va=rows[rows.cv_fold.astype(int)==fold]
        model,y,p,z,sids,epoch=train_fold(tr,va,cache,device,CFG.seed+fold,log)
        m=metrics(y,p); m.update({"fold":fold,"n_train":len(tr),"n_validation":len(va),
                                  "best_epoch":epoch}); fold_rows.append(m); best_epochs.append(epoch)
        pd.DataFrame(z,columns=[f"MFV{i:03d}" for i in range(1,129)]).assign(
            subject_id=sids, diagnosis_binary=y, fold=fold).to_csv(
            int_dir/f"fold_{fold}_MFV128_validation.csv",index=False)
        oof.append(pd.DataFrame({"subject_id":sids,"diagnosis_binary":y,
                                 "probability_AD":p,"fold":fold}))
        torch.save({"state_dict":model.state_dict(),"fold":fold,"best_epoch":epoch,
                    "config":asdict(CFG),"target_shape":TARGET_SHAPE},
                   model_dir/f"fold_{fold}_best.pt")
        del model; torch.cuda.empty_cache()
    oof=pd.concat(oof,ignore_index=True).sort_values("subject_id")
    if oof.subject_id.duplicated().any() or len(oof)!=len(rows): raise RuntimeError("OOF mismatch")
    threshold=choose_oof_threshold(oof.diagnosis_binary.values,oof.probability_AD.values)
    fold_df=pd.DataFrame(fold_rows); fold_df.to_csv(table_dir/"fold_metrics_threshold_0p50.csv",index=False)
    numeric=[c for c in fold_df if c not in {"fold","n_train","n_validation","best_epoch","tn","fp","fn","tp","threshold"}]
    summary={f"{c}_mean":float(fold_df[c].mean()) for c in numeric}
    summary.update({f"{c}_sd":float(fold_df[c].std(ddof=1)) for c in numeric})
    save_json(table_dir/"fivefold_mean_sd_threshold_0p50.json",summary)
    oof.to_csv(int_dir/"pooled_OOF_predictions.csv",index=False)
    save_json(table_dir/"pooled_OOF_threshold.json",
              metrics(oof.diagnosis_binary.values,oof.probability_AD.values,threshold))
    pd.concat([pd.read_csv(int_dir/f"fold_{f}_MFV128_validation.csv") for f in range(1,6)]).to_csv(
        int_dir/"pooled_OOF_MFV128.csv",index=False)
    figdir=ROOT/"Output/Figure/Main/Figure2_components"; figdir.mkdir(parents=True,exist_ok=True)
    panel="C" if cohort=="ADNI" else "D"
    stem=f"Figure2_component_{panel}_{cohort}_MRI_ROC_PR_development_OOF"
    plot_curves(oof.diagnosis_binary.values,oof.probability_AD.values,cohort,"five-fold OOF",
                figdir/f"{stem}.png",figdir/f"{stem}.pdf")
    save_json(qc_dir/"development_provenance.json",{"cohort":cohort,"run_tag":RUN_TAG,
        "locked_test_accessed":False,"split_manifest_sha256":sha256(SPLIT),
        "n_development":len(rows),"class_counts":rows.diagnosis.value_counts().to_dict(),
        "fold_counts":rows.cv_fold.value_counts().sort_index().to_dict(),
        "best_epochs":best_epochs,"config":asdict(CFG),"target_shape":TARGET_SHAPE})
    log.info("Development complete; locked test remains unopened.")

def train_final(rows, cache_dir, device, epochs, seed):
    """Final development-only fit; no validation or locked-test data."""
    set_seed(seed)
    dl = loader(rows, cache_dir, True, device)
    model = MRIBranch(CFG.dropout).to(device)
    n0 = int((rows.diagnosis_binary == 0).sum())
    n1 = int((rows.diagnosis_binary == 1).sum())
    loss_fn = nn.BCEWithLogitsLoss(
        pos_weight=torch.tensor([n0 / n1], device=device))
    opt = torch.optim.AdamW(model.parameters(), lr=CFG.lr,
                            weight_decay=CFG.weight_decay)
    amp = device.type == "cuda"
    scaler = torch.amp.GradScaler("cuda", enabled=amp)
    model.train()
    for _ in range(epochs):
        for x, y, _ in dl:
            x = x.to(device, non_blocking=True)
            y = y.to(device, non_blocking=True)
            x = augment_batch(x)
            opt.zero_grad(set_to_none=True)
            with torch.autocast(device_type=device.type, dtype=torch.float16,
                                enabled=amp):
                logits, _ = model(x)
                loss = loss_fn(logits, y)
            scaler.scale(loss).backward()
            scaler.step(opt)
            scaler.update()
    return model

def bootstrap_ci(y, p, threshold, n_bootstrap):
    rng = np.random.default_rng(CFG.seed + 5000)
    idx0, idx1 = np.where(y == 0)[0], np.where(y == 1)[0]
    keys = ["roc_auc", "pr_auc", "accuracy", "balanced_accuracy",
            "sensitivity", "specificity", "precision", "f1", "mcc", "brier"]
    values = {k: [] for k in keys}
    for _ in range(n_bootstrap):
        idx = np.r_[rng.choice(idx0, len(idx0), replace=True),
                    rng.choice(idx1, len(idx1), replace=True)]
        m = metrics(y[idx], p[idx], threshold)
        for k in keys:
            values[k].append(m[k])
    return {k: [float(np.quantile(v, 0.025)), float(np.quantile(v, 0.975))]
            for k, v in values.items()}

def locked_test(cohort, device, bootstrap_replicates, allow_locked_test):
    if not allow_locked_test:
        raise RuntimeError("Locked test requires --allow-locked-test.")

    qc_dir = outdir("QC", cohort)
    marker = qc_dir / "locked_test_opened_once.json"
    if marker.exists():
        raise RuntimeError(f"Refusing to reopen locked test: {marker}")

    table_dir = outdir("Table", cohort)
    fold_metrics = pd.read_csv(table_dir / "fold_metrics_threshold_0p50.csv")
    final_epochs = max(1, int(round(float(fold_metrics["best_epoch"].mean()))))
    threshold = float(json.loads(
        (table_dir / "pooled_OOF_threshold.json").read_text())["threshold"])

    save_json(marker, {
        "cohort": cohort,
        "status": "locked_test_opening_authorized",
        "final_epochs_from_mean_CV_best_epoch": final_epochs,
        "threshold_fixed_from_development_OOF": threshold,
    })

    log = setup_logger(cohort)
    log.warning("Opening %s MRI locked test exactly once.", cohort)

    development_rows = load_rows(cohort, "development")
    test_rows = load_rows(cohort, "locked_test")

    int_dir, model_dir = outdir("Intermediate", cohort), outdir("Model", cohort)
    dev_cache = int_dir / "development_cache"
    test_cache = int_dir / "locked_test_cache"

    model = train_final(development_rows, dev_cache, device, final_epochs,
                        CFG.seed + 900)
    y, prob, mfv, subject_ids = predict(
        model, loader(test_rows, test_cache, False, device), device)

    result = metrics(y, prob, threshold)
    ci = bootstrap_ci(y, prob, threshold, bootstrap_replicates)

    pd.DataFrame({
        "subject_id": subject_ids,
        "diagnosis_binary": y,
        "probability_AD": prob,
        "prediction_AD": (prob >= threshold).astype(int),
    }).to_csv(int_dir / "locked_test_predictions.csv", index=False)

    pd.DataFrame(mfv, columns=[f"MFV{i:03d}" for i in range(1, 129)]).assign(
        subject_id=subject_ids, diagnosis_binary=y
    ).to_csv(int_dir / "locked_test_MFV128.csv", index=False)

    torch.save({
        "state_dict": model.state_dict(),
        "final_epochs": final_epochs,
        "threshold_fixed_from_development_OOF": threshold,
        "config": asdict(CFG),
        "target_shape": TARGET_SHAPE,
    }, model_dir / "final_development_fit.pt")

    figdir = ROOT / "Output/Figure/Main/Figure2_components"
    panel = "C" if cohort == "ADNI" else "D"
    stem = f"Figure2_component_{panel}_{cohort}_MRI_ROC_PR_locked_test"
    plot_curves(y, prob, cohort, "locked test",
                figdir / f"{stem}.png", figdir / f"{stem}.pdf")

    report = {
        "status": "locked_test_complete",
        "cohort": cohort,
        "n": int(len(y)),
        "AD": int((y == 1).sum()),
        "CN": int((y == 0).sum()),
        "final_epochs_from_mean_CV_best_epoch": final_epochs,
        "threshold_fixed_from_development_OOF": threshold,
        "bootstrap_replicates": int(bootstrap_replicates),
        "metrics": result,
        "stratified_bootstrap_95ci": ci,
        "figure_components": [
            f"Output/Figure/Main/Figure2_components/{stem}.png",
            f"Output/Figure/Main/Figure2_components/{stem}.pdf",
        ],
    }
    save_json(table_dir / "locked_test_results.json", report)
    save_json(marker, report)
    log.info("Locked test complete: %s", json.dumps(report, sort_keys=True))

def validate(cohort):
    rows=load_rows(cohort,"development")
    report={"cohort":cohort,"development_n":len(rows),"AD":int((rows.diagnosis_binary==1).sum()),
      "CN":int((rows.diagnosis_binary==0).sum()),
      "fold_counts":{str(k):int(v) for k,v in rows.cv_fold.value_counts().sort_index().items()},
      "matrix_shape_after_preprocessing":[len(rows),1,*TARGET_SHAPE],
      "split_manifest_sha256":sha256(SPLIT),"locked_test_accessed":False}
    print(json.dumps(report,indent=2,sort_keys=True))

def main():
    ap=argparse.ArgumentParser(description="Leakage-controlled dual-cohort MRI branch.")
    ap.add_argument("command",choices=["VALIDATE","DEVELOP","TEST"])
    ap.add_argument("--cohort",choices=["ADNI","AddNeuroMed","ALL"],default="ALL")
    ap.add_argument("--device",choices=["auto","cpu","cuda"],default="auto")
    ap.add_argument("--allow-locked-test", action="store_true")
    ap.add_argument("--bootstrap-replicates", type=int, default=2000)
    args=ap.parse_args()
    if args.device=="auto": device=torch.device("cuda" if torch.cuda.is_available() else "cpu")
    else: device=torch.device(args.device)
    if device.type=="cuda" and not torch.cuda.is_available(): raise RuntimeError("CUDA requested but unavailable")
    cohorts=["ADNI","AddNeuroMed"] if args.cohort=="ALL" else [args.cohort]
    for cohort in cohorts:
        if args.command == "VALIDATE":
            validate(cohort)
        elif args.command == "DEVELOP":
            development(cohort, device)
        else:
            locked_test(cohort, device, args.bootstrap_replicates,
                        args.allow_locked_test)

if __name__ == "__main__":
    main()
