#!/usr/bin/env python3
"""Leakage-controlled dual-cohort top-500 gene-branch training."""
from __future__ import annotations

import argparse
import ast
import csv
import hashlib
import json
import logging
import math
import random
import sys
from collections import Counter
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import joblib
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import torch
from sklearn.calibration import calibration_curve
from sklearn.metrics import (
    accuracy_score, average_precision_score, balanced_accuracy_score,
    brier_score_loss, confusion_matrix, f1_score, matthews_corrcoef,
    precision_recall_curve, precision_score, roc_auc_score, roc_curve,
)
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from torch import nn
from torch.utils.data import DataLoader, TensorDataset

ROOT = Path(__file__).resolve().parents[1]
SPLIT_DIR = ROOT / "Output" / "Intermediate" / "Splits"
SPLIT_MANIFEST = SPLIT_DIR / "model_split_manifest.csv"
EXPECTED_SPLIT_VERSION = "dual_cohort_independent_v1"
COHORTS = ("ADNI", "AddNeuroMed")
SEED_OFFSETS = (1, 2, 3)
BASE_SEED = 20260716
N_FEATURES = 500


@dataclass(frozen=True)
class TrainConfig:
    hidden_dim: int = 256
    gfv_dim: int = 128
    classifier_dim: int = 64
    dropout_1: float = 0.20
    dropout_2: float = 0.15
    dropout_3: float = 0.10
    learning_rate: float = 5e-4
    weight_decay: float = 1e-4
    batch_size: int = 32
    max_epochs: int = 300
    patience: int = 35
    min_delta: float = 1e-4


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as fh:
        json.dump(value, fh, indent=2, sort_keys=True)


def setup_logging(path: Path) -> logging.Logger:
    path.parent.mkdir(parents=True, exist_ok=True)
    logger = logging.getLogger("gene_branch")
    logger.handlers.clear()
    logger.setLevel(logging.INFO)
    formatter = logging.Formatter("%(asctime)s | %(levelname)s | %(message)s")
    for handler in (logging.StreamHandler(sys.stdout), logging.FileHandler(path, mode="a")):
        handler.setFormatter(formatter)
        logger.addHandler(handler)
    return logger


def paths_for(cohort: str) -> tuple[Path, Path]:
    base = ROOT / "Output" / "Model_input" / "Gene"
    return (
        base / f"{cohort}_top500_gene_model_input.csv",
        base / f"{cohort}_top500_gene_order.txt",
    )


def run_dir(kind: str, tag: str, cohort: str) -> Path:
    return ROOT / "Output" / kind / "Gene_branch" / tag / cohort


def figure_dir() -> Path:
    return ROOT / "Output" / "Figure" / "Main" / "Figure2_components"


def seed_everything(seed: int) -> None:
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)
    torch.use_deterministic_algorithms(True, warn_only=True)


def device_from_arg(name: str) -> torch.device:
    if name == "cpu":
        return torch.device("cpu")
    if name == "cuda":
        if not torch.cuda.is_available():
            raise RuntimeError("--device cuda requested, but CUDA is unavailable")
        return torch.device("cuda")
    return torch.device("cuda" if torch.cuda.is_available() else "cpu")


def read_manifest() -> list[dict[str, str]]:
    if not SPLIT_MANIFEST.exists():
        raise FileNotFoundError(SPLIT_MANIFEST)
    with SPLIT_MANIFEST.open(newline="", encoding="utf-8-sig") as fh:
        rows = list(csv.DictReader(fh))
    required = {
        "split_version", "cohort", "model", "role", "cv_fold",
        "subject_id", "diagnosis_binary",
    }
    if not rows or not required.issubset(rows[0]):
        raise ValueError("Split manifest is missing required gene-branch columns")
    versions = {r["split_version"] for r in rows}
    if versions != {EXPECTED_SPLIT_VERSION}:
        raise ValueError(f"Unexpected split versions: {sorted(versions)}")
    return rows


def cohort_rows(cohort: str, role: str) -> list[dict[str, str]]:
    rows = [
        r for r in read_manifest()
        if r["cohort"] == cohort and r["model"] == "gene_branch" and r["role"] == role
    ]
    if not rows:
        raise ValueError(f"No gene-branch {role} rows for {cohort}")
    ids = [r["subject_id"] for r in rows]
    if len(ids) != len(set(ids)):
        raise ValueError(f"Duplicate gene-branch subject IDs for {cohort} {role}")
    return rows


def development_rows(cohort: str) -> list[dict[str, str]]:
    rows = cohort_rows(cohort, "development")
    folds = Counter(r["cv_fold"] for r in rows)
    if set(folds) != {"1", "2", "3", "4", "5"}:
        raise ValueError(f"{cohort} development folds are not exactly 1..5: {dict(folds)}")
    if any(r["diagnosis_binary"] not in {"0", "1"} for r in rows):
        raise ValueError(f"{cohort} development labels are not binary")
    return rows


def read_gene_order(path: Path) -> list[str]:
    genes = [line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
    if len(genes) != N_FEATURES or len(set(genes)) != N_FEATURES:
        raise ValueError(f"{path.name} must contain {N_FEATURES} unique genes")
    return genes


def read_matrix_for_ids(matrix_path: Path, expected_genes: list[str], records: list[dict[str, str]]) -> tuple[np.ndarray, np.ndarray, list[str], np.ndarray]:
    wanted = {r["subject_id"] for r in records}
    labels = {r["subject_id"]: int(r["diagnosis_binary"]) for r in records}
    folds = {r["subject_id"]: int(r["cv_fold"]) if r["cv_fold"] else -1 for r in records}
    found: dict[str, list[float]] = {}
    with matrix_path.open(newline="", encoding="utf-8-sig") as fh:
        reader = csv.reader(fh)
        header = next(reader)
        if len(header) != N_FEATURES + 1 or header[1:] != expected_genes:
            raise ValueError(f"{matrix_path.name} does not match its fixed 500-gene order")
        for row in reader:
            if not row:
                continue
            subject_id = row[0]
            if subject_id in wanted:
                if len(row) != N_FEATURES + 1:
                    raise ValueError(f"Malformed row for {subject_id} in {matrix_path.name}")
                found[subject_id] = [float(x) for x in row[1:]]
    missing = wanted.difference(found)
    if missing:
        raise ValueError(f"{matrix_path.name} is missing {len(missing)} required subject IDs")
    ids = [r["subject_id"] for r in records]
    X = np.asarray([found[sid] for sid in ids], dtype=np.float32)
    y = np.asarray([labels[sid] for sid in ids], dtype=np.int64)
    fold_values = np.asarray([folds[sid] for sid in ids], dtype=np.int64)
    if X.shape != (len(records), N_FEATURES) or not np.isfinite(X).all():
        raise ValueError("Selected gene matrix is not finite with the expected shape")
    return X, y, ids, fold_values


class GeneBranch(nn.Module):
    def __init__(self, cfg: TrainConfig) -> None:
        super().__init__()
        self.encoder = nn.Sequential(
            nn.Linear(N_FEATURES, cfg.hidden_dim), nn.LayerNorm(cfg.hidden_dim), nn.GELU(), nn.Dropout(cfg.dropout_1),
            nn.Linear(cfg.hidden_dim, cfg.gfv_dim), nn.LayerNorm(cfg.gfv_dim), nn.GELU(), nn.Dropout(cfg.dropout_2),
        )
        self.head = nn.Sequential(
            nn.Linear(cfg.gfv_dim, cfg.classifier_dim), nn.LayerNorm(cfg.classifier_dim), nn.GELU(), nn.Dropout(cfg.dropout_3),
            nn.Linear(cfg.classifier_dim, 1),
        )

    def forward(self, x: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        gfv = self.encoder(x)
        return self.head(gfv).squeeze(1), gfv


def predict(model: GeneBranch, X: np.ndarray, device: torch.device) -> tuple[np.ndarray, np.ndarray]:
    model.eval()
    probs, gfvs = [], []
    with torch.no_grad():
        for start in range(0, len(X), 256):
            xb = torch.from_numpy(X[start:start + 256]).to(device)
            logits, gfv = model(xb)
            probs.append(torch.sigmoid(logits).cpu().numpy())
            gfvs.append(gfv.cpu().numpy())
    return np.concatenate(probs), np.concatenate(gfvs)


def train_seed(
    X_train: np.ndarray,
    y_train: np.ndarray,
    X_val: np.ndarray,
    y_val: np.ndarray,
    cfg: TrainConfig,
    seed: int,
    device: torch.device,
) -> tuple[GeneBranch, int, list[dict[str, float]]]:
    seed_everything(seed)
    model = GeneBranch(cfg).to(device)
    n_ad = int(y_train.sum())
    n_cn = int((y_train == 0).sum())
    if n_ad == 0 or n_cn == 0:
        raise ValueError("A training fold must contain both AD and CN subjects")

    loss_fn = nn.BCEWithLogitsLoss(
        pos_weight=torch.tensor(n_cn / n_ad, dtype=torch.float32, device=device)
    )
    optimizer = torch.optim.AdamW(
        model.parameters(),
        lr=cfg.learning_rate,
        weight_decay=cfg.weight_decay,
    )
    data = TensorDataset(
        torch.from_numpy(X_train),
        torch.from_numpy(y_train.astype(np.float32)),
    )
    generator = torch.Generator().manual_seed(seed)
    loader = DataLoader(
        data,
        batch_size=min(cfg.batch_size, len(data)),
        shuffle=True,
        generator=generator,
    )

    best_state = None
    best_auc = -np.inf
    best_epoch = 0
    stale = 0
    history = []

    for epoch in range(1, cfg.max_epochs + 1):
        model.train()
        losses = []

        for xb, yb in loader:
            xb = xb.to(device)
            yb = yb.to(device)
            optimizer.zero_grad(set_to_none=True)
            logits, _ = model(xb)
            loss = loss_fn(logits, yb)
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=5.0)
            optimizer.step()
            losses.append(float(loss.item()))

        val_prob, _ = predict(model, X_val, device)
        val_auc = float(roc_auc_score(y_val, val_prob))
        history.append(
            {
                "epoch": epoch,
                "train_loss": float(np.mean(losses)),
                "validation_roc_auc": val_auc,
            }
        )

        if val_auc > best_auc + cfg.min_delta:
            best_state = {
                key: value.detach().cpu().clone()
                for key, value in model.state_dict().items()
            }
            best_auc = val_auc
            best_epoch = epoch
            stale = 0
        else:
            stale += 1

        if stale >= cfg.patience:
            break

    if best_state is None:
        raise RuntimeError("No valid development checkpoint was selected")

    model.load_state_dict(best_state)
    return model, best_epoch, history


def fit_fixed_epochs(
    X: np.ndarray,
    y: np.ndarray,
    cfg: TrainConfig,
    seed: int,
    epochs: int,
    device: torch.device,
) -> GeneBranch:
    seed_everything(seed)
    model = GeneBranch(cfg).to(device)
    n_ad = int(y.sum())
    n_cn = int((y == 0).sum())
    loss_fn = nn.BCEWithLogitsLoss(
        pos_weight=torch.tensor(n_cn / n_ad, dtype=torch.float32, device=device)
    )
    optimizer = torch.optim.AdamW(
        model.parameters(),
        lr=cfg.learning_rate,
        weight_decay=cfg.weight_decay,
    )
    data = TensorDataset(
        torch.from_numpy(X),
        torch.from_numpy(y.astype(np.float32)),
    )
    loader = DataLoader(
        data,
        batch_size=min(cfg.batch_size, len(data)),
        shuffle=True,
        generator=torch.Generator().manual_seed(seed),
    )

    for _ in range(epochs):
        model.train()
        for xb, yb in loader:
            xb = xb.to(device)
            yb = yb.to(device)
            optimizer.zero_grad(set_to_none=True)
            logits, _ = model(xb)
            loss = loss_fn(logits, yb)
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=5.0)
            optimizer.step()

    return model


def metrics_at_threshold(
    y: np.ndarray,
    probability: np.ndarray,
    threshold: float = 0.50,
) -> dict[str, float]:
    pred = (probability >= threshold).astype(int)
    tn, fp, fn, tp = confusion_matrix(y, pred, labels=[0, 1]).ravel()

    bins = np.linspace(0, 1, 11)
    ece = 0.0
    for low, high in zip(bins[:-1], bins[1:]):
        selected = (probability >= low) & (
            (probability < high) if high < 1 else (probability <= high)
        )
        if selected.any():
            ece += selected.mean() * abs(
                float(y[selected].mean()) - float(probability[selected].mean())
            )

    return {
        "roc_auc": float(roc_auc_score(y, probability)),
        "pr_auc": float(average_precision_score(y, probability)),
        "accuracy": float(accuracy_score(y, pred)),
        "balanced_accuracy": float(balanced_accuracy_score(y, pred)),
        "sensitivity": float(tp / (tp + fn)),
        "specificity": float(tn / (tn + fp)),
        "precision": float(precision_score(y, pred, zero_division=0)),
        "f1": float(f1_score(y, pred, zero_division=0)),
        "mcc": float(matthews_corrcoef(y, pred)),
        "brier": float(brier_score_loss(y, probability)),
        "ece": float(ece),
        "tn": int(tn),
        "fp": int(fp),
        "fn": int(fn),
        "tp": int(tp),
    }


def select_development_threshold(
    y: np.ndarray,
    probability: np.ndarray,
) -> dict[str, float]:
    candidates = np.unique(np.r_[0.0, 0.50, 1.0, probability])
    best = None

    for threshold in candidates:
        values = metrics_at_threshold(y, probability, float(threshold))
        candidate = (
            values["balanced_accuracy"],
            -abs(float(threshold) - 0.50),
            float(threshold),
        )
        if best is None or candidate > best[0]:
            best = (candidate, values)

    assert best is not None
    return {"threshold": best[0][2], **best[1]}


def mean_sd_metrics(frame: pd.DataFrame) -> dict[str, float]:
    output: dict[str, float] = {}
    metric_names = [
        "roc_auc",
        "pr_auc",
        "accuracy",
        "balanced_accuracy",
        "sensitivity",
        "specificity",
        "precision",
        "f1",
        "mcc",
        "brier",
        "ece",
    ]

    for name in metric_names:
        output[f"{name}_mean"] = float(frame[name].mean())
        output[f"{name}_sd"] = float(frame[name].std(ddof=1))

    return output


def save_gfv(
    path: Path,
    ids: list[str],
    y: np.ndarray,
    folds: np.ndarray,
    gfv: np.ndarray,
) -> None:
    frame = pd.DataFrame(
        gfv,
        columns=[f"GFV{i:03d}" for i in range(1, 129)],
    )
    frame.insert(0, "diagnosis_binary", y)
    frame.insert(0, "cv_fold", folds)
    frame.insert(0, "subject_id", ids)
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, index=False)


def save_oof_figure(
    cohort: str,
    y: np.ndarray,
    probability: np.ndarray,
) -> list[str]:
    color = "#2166AC" if cohort == "ADNI" else "#B2182B"
    panel = "A" if cohort == "ADNI" else "B"

    plt.rcParams.update(
        {"font.family": "DejaVu Sans", "font.size": 9, "axes.linewidth": 0.8}
    )
    fig, axes = plt.subplots(
        1, 2, figsize=(7.2, 3.1), constrained_layout=True
    )

    fpr, tpr, _ = roc_curve(y, probability)
    auc = roc_auc_score(y, probability)
    axes[0].plot(fpr, tpr, color=color, lw=2, label=f"OOF AUC = {auc:.3f}")
    axes[0].plot([0, 1], [0, 1], "--", color="#777777", lw=0.8)
    axes[0].set(
        xlabel="False-positive rate",
        ylabel="True-positive rate",
        xlim=(0, 1),
        ylim=(0, 1),
        title=f"{cohort} gene branch: ROC",
    )
    axes[0].legend(frameon=False, loc="lower right")

    precision, recall, _ = precision_recall_curve(y, probability)
    ap = average_precision_score(y, probability)
    axes[1].plot(recall, precision, color=color, lw=2, label=f"OOF AP = {ap:.3f}")
    axes[1].axhline(
        float(y.mean()),
        ls="--",
        color="#777777",
        lw=0.8,
        label="Prevalence",
    )
    axes[1].set(
        xlabel="Recall",
        ylabel="Precision",
        xlim=(0, 1),
        ylim=(0, 1),
        title=f"{cohort} gene branch: PR",
    )
    axes[1].legend(frameon=False, loc="lower left")

    output = figure_dir()
    output.mkdir(parents=True, exist_ok=True)
    stem = output / (
        f"Figure2_component_{panel}_{cohort}_gene_ROC_PR_development_OOF"
    )
    png = stem.with_suffix(".png")
    pdf = stem.with_suffix(".pdf")
    fig.savefig(png, dpi=600, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)
    return [str(p.relative_to(ROOT)) for p in (png, pdf)]


def validate_cohort(cohort: str) -> dict[str, Any]:
    records = development_rows(cohort)
    matrix, order = paths_for(cohort)
    X, y, ids, folds = read_matrix_for_ids(
        matrix,
        read_gene_order(order),
        records,
    )
    expected = {"ADNI": 242, "AddNeuroMed": 280}[cohort]
    if len(ids) != expected:
        raise ValueError(
            f"{cohort} development count is {len(ids)}, expected {expected}"
        )

    return {
        "cohort": cohort,
        "development_n": len(ids),
        "AD": int(y.sum()),
        "CN": int((y == 0).sum()),
        "matrix_shape": list(X.shape),
        "fold_counts": {str(f): int((folds == f).sum()) for f in range(1, 6)},
        "matrix_sha256": sha256(matrix),
        "gene_order_sha256": sha256(order),
        "split_manifest_sha256": sha256(SPLIT_MANIFEST),
        "locked_test_accessed": False,
    }


def develop_cohort(
    cohort: str,
    tag: str,
    cfg: TrainConfig,
    device: torch.device,
) -> None:
    validation = validate_cohort(cohort)

    dirs = {
        kind: run_dir(kind, tag, cohort)
        for kind in ("Model", "Intermediate", "Table", "QC", "Log")
    }
    for directory in dirs.values():
        directory.mkdir(parents=True, exist_ok=True)

    marker = dirs["QC"] / "development_complete.json"
    if marker.exists():
        raise FileExistsError(
            f"Completed development run exists: {marker}; use a new --run-tag"
        )

    logger = setup_logging(dirs["Log"] / "development.log")
    logger.info("%s development-only run started on %s", cohort, device)
    logger.info(
        "No locked-test IDs, values, labels, predictions, or metrics are loaded."
    )

    write_json(
        dirs["QC"] / "development_provenance.json",
        {
            "run_tag": tag,
            "cohort": cohort,
            "split_version": EXPECTED_SPLIT_VERSION,
            "architecture": "500 -> 256 -> GFV128 -> 64 -> AD/CN",
            "scaler": "fit only on each four-fold training partition",
            "loss": "class-weighted BCEWithLogitsLoss",
            "seeds_per_fold": 3,
            "config": asdict(cfg),
            "validation": validation,
        },
    )

    records = development_rows(cohort)
    matrix, order = paths_for(cohort)
    X, y, ids, folds = read_matrix_for_ids(
        matrix,
        read_gene_order(order),
        records,
    )

    oof_prob = np.full(len(ids), np.nan)
    oof_gfv = np.full((len(ids), 128), np.nan, dtype=np.float32)
    fold_metrics = []
    seed_summary = []
    best_epochs = []

    for fold in range(1, 6):
        train_idx = folds != fold
        val_idx = folds == fold

        scaler = StandardScaler().fit(X[train_idx])
        X_train = scaler.transform(X[train_idx]).astype(np.float32)
        X_val = scaler.transform(X[val_idx]).astype(np.float32)

        fold_dir = dirs["Model"] / f"fold_{fold}"
        fold_dir.mkdir(parents=True, exist_ok=True)
        joblib.dump(scaler, fold_dir / "standard_scaler.joblib")

        probabilities = []
        representations = []

        for seed_number, offset in enumerate(SEED_OFFSETS, start=1):
            seed = BASE_SEED + fold * 100 + offset
            model, epoch, history = train_seed(
                X_train,
                y[train_idx],
                X_val,
                y[val_idx],
                cfg,
                seed,
                device,
            )
            probability, gfv = predict(model, X_val, device)

            torch.save(
                {
                    "state_dict": model.cpu().state_dict(),
                    "seed": seed,
                    "best_epoch": epoch,
                    "config": asdict(cfg),
                },
                fold_dir / f"seed_{seed_number}_gene_encoder_classifier.pt",
            )
            pd.DataFrame(history).to_csv(
                fold_dir / f"seed_{seed_number}_history.csv",
                index=False,
            )

            probabilities.append(probability)
            representations.append(gfv)
            best_epochs.append(epoch)
            seed_summary.append(
                {
                    "cohort": cohort,
                    "fold": fold,
                    "seed": seed,
                    "best_epoch": epoch,
                }
            )

        ensemble_prob = np.mean(np.stack(probabilities), axis=0)
        ensemble_gfv = np.mean(np.stack(representations), axis=0)
        oof_prob[val_idx] = ensemble_prob
        oof_gfv[val_idx] = ensemble_gfv

        metric = metrics_at_threshold(y[val_idx], ensemble_prob)
        fold_metrics.append(
            {
                "cohort": cohort,
                "fold": fold,
                "n": int(val_idx.sum()),
                "AD": int(y[val_idx].sum()),
                "CN": int((y[val_idx] == 0).sum()),
                **metric,
            }
        )
        logger.info(
            "fold=%d n=%d AUC=%.3f AP=%.3f BA=%.3f",
            fold,
            int(val_idx.sum()),
            metric["roc_auc"],
            metric["pr_auc"],
            metric["balanced_accuracy"],
        )

    if not np.isfinite(oof_prob).all() or not np.isfinite(oof_gfv).all():
        raise RuntimeError("Incomplete pooled OOF output")

    fold_frame = pd.DataFrame(fold_metrics)
    seed_frame = pd.DataFrame(seed_summary)
    threshold = select_development_threshold(y, oof_prob)

    pd.DataFrame(
        {
            "subject_id": ids,
            "cv_fold": folds,
            "diagnosis_binary": y,
            "oof_probability": oof_prob,
        }
    ).to_csv(dirs["Intermediate"] / "development_OOF_predictions.csv", index=False)

    save_gfv(
        dirs["Intermediate"] / "development_OOF_GFV128.csv",
        ids,
        y,
        folds,
        oof_gfv,
    )

    fold_frame.to_csv(
        dirs["Table"] / "fivefold_metrics_threshold_0p50.csv",
        index=False,
    )
    seed_frame.to_csv(
        dirs["Table"] / "seed_training_summary.csv",
        index=False,
    )
    write_json(
        dirs["Table"] / "fivefold_mean_sd_threshold_0p50.json",
        mean_sd_metrics(fold_frame),
    )
    write_json(dirs["Table"] / "pooled_OOF_threshold.json", threshold)

    with pd.ExcelWriter(
        dirs["Table"] / "development_performance.xlsx",
        engine="openpyxl",
    ) as writer:
        fold_frame.to_excel(writer, sheet_name="Fold_metrics", index=False)
        seed_frame.to_excel(writer, sheet_name="Seed_summary", index=False)
        pd.DataFrame([mean_sd_metrics(fold_frame)]).to_excel(
            writer,
            sheet_name="Mean_SD",
            index=False,
        )

    final_epochs = max(1, int(round(float(np.median(best_epochs)))))
    final_scaler = StandardScaler().fit(X)
    X_all = final_scaler.transform(X).astype(np.float32)

    final_dir = dirs["Model"] / "final_development_models"
    final_dir.mkdir(parents=True, exist_ok=True)
    joblib.dump(final_scaler, final_dir / "standard_scaler.joblib")

    final_gfvs = []
    for seed_number, offset in enumerate(SEED_OFFSETS, start=1):
        model = fit_fixed_epochs(
            X_all,
            y,
            cfg,
            BASE_SEED + 900 + offset,
            final_epochs,
            device,
        )
        _, gfv = predict(model, X_all, device)
        torch.save(
            {
                "state_dict": model.cpu().state_dict(),
                "seed": BASE_SEED + 900 + offset,
                "epochs": final_epochs,
                "config": asdict(cfg),
            },
            final_dir / f"seed_{seed_number}_gene_encoder_classifier.pt",
        )
        final_gfvs.append(gfv)

    save_gfv(
        dirs["Intermediate"] / "development_final_GFV128.csv",
        ids,
        y,
        folds,
        np.mean(np.stack(final_gfvs), axis=0),
    )

    state = {
        "cohort": cohort,
        "run_tag": tag,
        "final_epochs": final_epochs,
        "pooled_OOF_threshold": threshold["threshold"],
        "figure_components": save_oof_figure(cohort, y, oof_prob),
        "locked_test_accessed": False,
        "status": "development_complete",
    }
    write_json(marker, state)
    logger.info("Development complete; locked test remains unopened.")


def write_summary(tag: str) -> None:
    summaries = []
    sheets = {}

    for cohort in COHORTS:
        table = run_dir(
            "Table",
            tag,
            cohort,
        ) / "fivefold_metrics_threshold_0p50.csv"

        if table.exists():
            sheets[cohort] = pd.read_csv(table)
            summaries.append(
                {"cohort": cohort, **mean_sd_metrics(sheets[cohort])}
            )

    if not summaries:
        raise FileNotFoundError("No development metrics found for this --run-tag")

    summary = pd.DataFrame(summaries)
    out = (
        ROOT
        / "Output"
        / "Table"
        / "Supplementary"
        / "Supplementary_Table_6_gene_branch_performance.xlsx"
    )
    out.parent.mkdir(parents=True, exist_ok=True)

    with pd.ExcelWriter(out, engine="openpyxl") as writer:
        summary.to_excel(writer, sheet_name="Fivefold_summary", index=False)
        for cohort, frame in sheets.items():
            frame.to_excel(writer, sheet_name=f"{cohort}_folds", index=False)

    print(summary.to_string(index=False))
    print(f"Saved: {out}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "mode",
        choices=("VALIDATE", "DEVELOP", "SUMMARY", "TEST"),
    )
    parser.add_argument(
        "--cohort",
        choices=("ADNI", "AddNeuroMed", "ALL"),
        default="ALL",
    )
    parser.add_argument(
        "--run-tag",
        default="dual_cohort_top500_v1",
    )
    parser.add_argument(
        "--device",
        choices=("auto", "cpu", "cuda"),
        default="auto",
    )
    args = parser.parse_args()

    cohorts = COHORTS if args.cohort == "ALL" else (args.cohort,)

    if args.mode == "VALIDATE":
        report = {cohort: validate_cohort(cohort) for cohort in cohorts}
        print(json.dumps(report, indent=2, sort_keys=True))
        return

    if args.mode == "SUMMARY":
        write_summary(args.run_tag)
        return

    if args.mode == "TEST":
        raise SystemExit(
            "Locked-test evaluation is disabled until development results are reviewed."
        )

    device = device_from_arg(args.device)
    for cohort in cohorts:
        develop_cohort(cohort, args.run_tag, TrainConfig(), device)


def bootstrap_intervals(
    y: np.ndarray,
    probability: np.ndarray,
    threshold: float,
    replicates: int = 2000,
) -> dict[str, list[float]]:
    rng = np.random.default_rng(BASE_SEED + 777)
    zero = np.flatnonzero(y == 0)
    one = np.flatnonzero(y == 1)
    names = (
        "roc_auc",
        "pr_auc",
        "accuracy",
        "balanced_accuracy",
        "sensitivity",
        "specificity",
        "precision",
        "f1",
        "mcc",
        "brier",
    )
    values = {name: [] for name in names}

    for _ in range(replicates):
        idx = np.r_[
            rng.choice(zero, len(zero), replace=True),
            rng.choice(one, len(one), replace=True),
        ]
        metric = metrics_at_threshold(y[idx], probability[idx], threshold)
        for name in names:
            values[name].append(metric[name])

    return {
        name: [
            float(np.percentile(value, 2.5)),
            float(np.percentile(value, 97.5)),
        ]
        for name, value in values.items()
    }


def save_locked_test_figure(
    cohort: str,
    y: np.ndarray,
    probability: np.ndarray,
) -> list[str]:
    color = "#2166AC" if cohort == "ADNI" else "#B2182B"
    panel = "A" if cohort == "ADNI" else "B"

    fig, axes = plt.subplots(
        1,
        2,
        figsize=(7.2, 3.1),
        constrained_layout=True,
    )

    fpr, tpr, _ = roc_curve(y, probability)
    axes[0].plot(
        fpr,
        tpr,
        color=color,
        lw=2,
        label=f"Test AUC = {roc_auc_score(y, probability):.3f}",
    )
    axes[0].plot([0, 1], [0, 1], "--", color="#777777", lw=0.8)
    axes[0].set(
        xlabel="False-positive rate",
        ylabel="True-positive rate",
        xlim=(0, 1),
        ylim=(0, 1),
        title=f"{cohort} gene branch: locked test ROC",
    )
    axes[0].legend(frameon=False, loc="lower right")

    precision, recall, _ = precision_recall_curve(y, probability)
    axes[1].plot(
        recall,
        precision,
        color=color,
        lw=2,
        label=f"Test AP = {average_precision_score(y, probability):.3f}",
    )
    axes[1].axhline(
        float(y.mean()),
        ls="--",
        color="#777777",
        lw=0.8,
        label="Prevalence",
    )
    axes[1].set(
        xlabel="Recall",
        ylabel="Precision",
        xlim=(0, 1),
        ylim=(0, 1),
        title=f"{cohort} gene branch: locked test PR",
    )
    axes[1].legend(frameon=False, loc="lower left")

    output = figure_dir()
    output.mkdir(parents=True, exist_ok=True)
    stem = output / f"Figure2_component_{panel}_{cohort}_gene_ROC_PR_locked_test"
    png = stem.with_suffix(".png")
    pdf = stem.with_suffix(".pdf")
    fig.savefig(png, dpi=600, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)

    return [str(p.relative_to(ROOT)) for p in (png, pdf)]


def run_locked_test(
    cohort: str,
    tag: str,
    device: torch.device,
    replicates: int,
) -> None:
    dirs = {
        kind: run_dir(kind, tag, cohort)
        for kind in ("Model", "Intermediate", "Table", "QC", "Log")
    }
    development = dirs["QC"] / "development_complete.json"
    marker = dirs["QC"] / "locked_test_complete.json"

    if not development.exists():
        raise FileNotFoundError(f"Development is not frozen: {development}")
    if marker.exists():
        raise FileExistsError(f"Locked test already evaluated: {marker}")

    state = json.loads(development.read_text())
    threshold = float(state["pooled_OOF_threshold"])
    logger = setup_logging(dirs["Log"] / "locked_test.log")
    logger.warning("Opening %s locked test exactly once.", cohort)

    records = cohort_rows(cohort, "locked_test")
    matrix, order = paths_for(cohort)
    X, y, ids, folds = read_matrix_for_ids(
        matrix,
        read_gene_order(order),
        records,
    )

    if not np.all(folds == -1):
        raise ValueError("Locked-test rows unexpectedly have CV folds")

    final_dir = dirs["Model"] / "final_development_models"
    scaler = joblib.load(final_dir / "standard_scaler.joblib")
    X = scaler.transform(X).astype(np.float32)

    probabilities = []
    gfvs = []

    for seed in range(1, 4):
        checkpoint = torch.load(
            final_dir / f"seed_{seed}_gene_encoder_classifier.pt",
            map_location=device,
            weights_only=False,
        )
        model = GeneBranch(TrainConfig()).to(device)
        model.load_state_dict(checkpoint["state_dict"])
        probability, gfv = predict(model, X, device)
        probabilities.append(probability)
        gfvs.append(gfv)

    probability = np.mean(np.stack(probabilities), axis=0)
    gfv = np.mean(np.stack(gfvs), axis=0)
    metric = metrics_at_threshold(y, probability, threshold)

    report = {
        "cohort": cohort,
        "n": int(len(ids)),
        "AD": int(y.sum()),
        "CN": int((y == 0).sum()),
        "threshold_fixed_from_development_OOF": threshold,
        "metrics": metric,
        "stratified_bootstrap_95ci": bootstrap_intervals(
            y,
            probability,
            threshold,
            replicates,
        ),
        "bootstrap_replicates": replicates,
        "figure_components": save_locked_test_figure(
            cohort,
            y,
            probability,
        ),
        "status": "locked_test_complete",
    }

    pd.DataFrame(
        {
            "subject_id": ids,
            "diagnosis_binary": y,
            "probability": probability,
            "threshold": threshold,
        }
    ).to_csv(dirs["Intermediate"] / "locked_test_predictions.csv", index=False)

    save_gfv(
        dirs["Intermediate"] / "locked_test_GFV128.csv",
        ids,
        y,
        folds,
        gfv,
    )
    write_json(dirs["Table"] / "locked_test_metrics_95ci.json", report)
    write_json(marker, report)
    logger.info("Locked test complete: %s", json.dumps(report, sort_keys=True))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "mode",
        choices=("VALIDATE", "DEVELOP", "SUMMARY", "TEST"),
    )
    parser.add_argument(
        "--cohort",
        choices=("ADNI", "AddNeuroMed", "ALL"),
        default="ALL",
    )
    parser.add_argument(
        "--run-tag",
        default="dual_cohort_top500_v1",
    )
    parser.add_argument(
        "--device",
        choices=("auto", "cpu", "cuda"),
        default="auto",
    )
    parser.add_argument("--allow-locked-test", action="store_true")
    parser.add_argument("--bootstrap-replicates", type=int, default=2000)
    args = parser.parse_args()

    cohorts = COHORTS if args.cohort == "ALL" else (args.cohort,)

    if args.mode == "VALIDATE":
        report = {cohort: validate_cohort(cohort) for cohort in cohorts}
        print(json.dumps(report, indent=2, sort_keys=True))
        return

    if args.mode == "SUMMARY":
        write_summary(args.run_tag)
        return

    device = device_from_arg(args.device)

    if args.mode == "DEVELOP":
        for cohort in cohorts:
            develop_cohort(cohort, args.run_tag, TrainConfig(), device)
        return

    if not args.allow_locked_test:
        raise SystemExit(
            "TEST requires --allow-locked-test after development review."
        )

    for cohort in cohorts:
        run_locked_test(
            cohort,
            args.run_tag,
            device,
            args.bootstrap_replicates,
        )


if __name__ == "__main__":
    main()
