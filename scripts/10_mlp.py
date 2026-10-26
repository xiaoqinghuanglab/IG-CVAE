#!/usr/bin/env python3
"""
Deterministic MLP baseline for IG-CVAE comparison.

Frozen structural-MRI MFV128 -> deterministic MLP -> generated gene GFV128

The baseline deliberately uses the SAME:
    - participant splits
    - fold assignments
    - frozen gene-branch encoders
    - frozen MRI-branch encoders
    - fold-specific representation strategy
    - train-fold-only MFV/GFV scaling
    - downstream logistic-regression evaluation
    - reconstruction metrics
    - bootstrap procedure

as scripts/06_run_cvae.py.

This allows a fair deterministic MLP comparison against IG-CVAE.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import logging
import os
import random
import sys
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

os.environ.setdefault("CUBLAS_WORKSPACE_CONFIG", ":4096:8")

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F

from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from torch.utils.data import DataLoader, TensorDataset


# ============================================================
# Project paths
# ============================================================

ROOT = Path(__file__).resolve().parents[1]

CVAE_SCRIPT = ROOT / "scripts" / "06_run_cvae.py"

MODEL_NAME = "MLP"
RUN_TAG = "deterministic_mlp_faircv_v1"


# ============================================================
# Import current CVAE script
# ============================================================

spec = importlib.util.spec_from_file_location(
    "igcvae_primary_for_mlp",
    CVAE_SCRIPT,
)

if spec is None or spec.loader is None:
    raise RuntimeError(
        f"Could not import current CVAE script: {CVAE_SCRIPT}"
    )

CVAE = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = CVAE
spec.loader.exec_module(CVAE)


# ============================================================
# Configuration
# ============================================================

@dataclass(frozen=True)
class Config:
    hidden_dim: int = 256
    dropout: float = 0.10

    learning_rate: float = 1e-3
    weight_decay: float = 1e-4

    batch_size: int = 32
    max_epochs: int = 200
    patience: int = 20

    # Preserved from the uploaded baseline.
    #
    # Total loss:
    #     MSE + cosine_weight * (1 - cosine_similarity)
    #
    # Set to 0.0 later if you want a pure-MSE deterministic
    # regression baseline.
    cosine_weight: float = 0.20

    seed: int = 20261010


CFG = Config()

N_FOLDS = 5
MFV_DIM = 128
GFV_DIM = 128


# ============================================================
# Output organization
# ============================================================

def outdir(kind: str, cohort: str) -> Path:
    """
    Match the current project organization:

    Output/
        Model/MLP/<RUN_TAG>/<cohort>/
        Intermediate/MLP/<RUN_TAG>/<cohort>/
        Table/MLP/<RUN_TAG>/<cohort>/
        QC/MLP/<RUN_TAG>/<cohort>/
        Log/MLP/<RUN_TAG>/<cohort>/
    """

    path = (
        ROOT
        / "Output"
        / kind
        / MODEL_NAME
        / RUN_TAG
        / cohort
    )

    path.mkdir(
        parents=True,
        exist_ok=True,
    )

    return path


def dirs_for(cohort: str) -> dict[str, Path]:
    return {
        kind: outdir(kind, cohort)
        for kind in (
            "Model",
            "Intermediate",
            "Table",
            "QC",
            "Log",
        )
    }


# ============================================================
# Utilities
# ============================================================

def set_seed(seed: int) -> None:
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)

    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)

    torch.backends.cudnn.benchmark = False
    torch.backends.cudnn.deterministic = True

    torch.use_deterministic_algorithms(
        True,
        warn_only=True,
    )


def save_json(
    path: Path,
    obj: Any,
) -> None:

    def convert(x):
        if isinstance(x, np.integer):
            return int(x)

        if isinstance(x, np.floating):
            return float(x)

        if isinstance(x, np.ndarray):
            return x.tolist()

        if isinstance(x, Path):
            return str(x)

        raise TypeError(
            f"Cannot serialize {type(x).__name__}"
        )

    path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    with path.open(
        "w",
        encoding="utf-8",
    ) as f:

        json.dump(
            obj,
            f,
            indent=2,
            sort_keys=True,
            default=convert,
        )


def setup_logger(
    cohort: str,
    phase: str,
) -> logging.Logger:

    log_dir = outdir(
        "Log",
        cohort,
    )

    logger = logging.getLogger(
        f"mlp_{cohort}_{phase}"
    )

    logger.handlers.clear()
    logger.setLevel(
        logging.INFO
    )

    formatter = logging.Formatter(
        "%(asctime)s | %(levelname)s | %(message)s"
    )

    stream = logging.StreamHandler(
        sys.stdout
    )

    file_handler = logging.FileHandler(
        log_dir / f"{phase}.log",
        mode="a",
    )

    stream.setFormatter(
        formatter
    )

    file_handler.setFormatter(
        formatter
    )

    logger.addHandler(
        stream
    )

    logger.addHandler(
        file_handler
    )

    return logger


def cpu_state_dict(
    model: nn.Module,
) -> dict[str, torch.Tensor]:

    return {
        key: value.detach().cpu().clone()
        for key, value
        in model.state_dict().items()
    }


# ============================================================
# Deterministic MLP
# ============================================================

class DeterministicMLP(nn.Module):

    def __init__(self) -> None:
        super().__init__()

        self.net = nn.Sequential(

            nn.Linear(
                MFV_DIM,
                CFG.hidden_dim,
            ),

            nn.LayerNorm(
                CFG.hidden_dim
            ),

            nn.GELU(),

            nn.Dropout(
                CFG.dropout
            ),

            nn.Linear(
                CFG.hidden_dim,
                CFG.hidden_dim,
            ),

            nn.LayerNorm(
                CFG.hidden_dim
            ),

            nn.GELU(),

            nn.Dropout(
                CFG.dropout
            ),

            nn.Linear(
                CFG.hidden_dim,
                GFV_DIM,
            ),
        )

    def forward(
        self,
        x: torch.Tensor,
    ) -> torch.Tensor:

        return self.net(x)


# ============================================================
# Loss
# ============================================================

def mlp_loss(
    pred: torch.Tensor,
    target: torch.Tensor,
) -> torch.Tensor:
    """Pure MSE loss for the deterministic MLP baseline."""

    return F.mse_loss(
        pred,
        target,
    )


# ============================================================
# Training
# ============================================================

def train_model(
    X_train: np.ndarray,
    Y_train: np.ndarray,
    X_val: np.ndarray,
    Y_val: np.ndarray,
    device: torch.device,
    seed: int,
) -> tuple[
    DeterministicMLP,
    int,
    list[dict[str, float]],
]:

    set_seed(seed)

    model = DeterministicMLP().to(
        device
    )

    optimizer = torch.optim.AdamW(
        model.parameters(),
        lr=CFG.learning_rate,
        weight_decay=CFG.weight_decay,
    )

    train_dataset = TensorDataset(
        torch.from_numpy(
            X_train
        ),
        torch.from_numpy(
            Y_train
        ),
    )

    train_loader = DataLoader(
        train_dataset,
        batch_size=min(
            CFG.batch_size,
            len(train_dataset),
        ),
        shuffle=True,
        generator=torch.Generator().manual_seed(
            seed
        ),
    )

    xv = torch.from_numpy(
        X_val
    ).to(
        device
    )

    yv = torch.from_numpy(
        Y_val
    ).to(
        device
    )

    best_state = None
    best_epoch = 0
    best_val_mse = np.inf

    stale = 0
    history: list[
        dict[str, float]
    ] = []

    for epoch in range(
        1,
        CFG.max_epochs + 1,
    ):

        model.train()

        batch_mse = []

        for xb, yb in train_loader:

            xb = xb.to(
                device
            )

            yb = yb.to(
                device
            )

            optimizer.zero_grad(
                set_to_none=True
            )

            pred = model(
                xb
            )

            mse = mlp_loss(
                pred,
                yb,
            )

            mse.backward()

            torch.nn.utils.clip_grad_norm_(
                model.parameters(),
                max_norm=5.0,
            )

            optimizer.step()

            batch_mse.append(
                float(
                    mse.item()
                )
            )

        model.eval()

        with torch.no_grad():

            val_pred = model(
                xv
            )

            val_mse = float(
                F.mse_loss(
                    val_pred,
                    yv,
                ).item()
            )

        history.append({
            "epoch": int(epoch),

            "train_mse":
                float(
                    np.mean(
                        batch_mse
                    )
                ),

            "validation_mse":
                val_mse,
        })

        if (
            val_mse
            < best_val_mse - 1e-6
        ):

            best_val_mse = (
                val_mse
            )

            best_epoch = (
                epoch
            )

            best_state = (
                cpu_state_dict(
                    model
                )
            )

            stale = 0

        else:

            stale += 1

        if (
            stale
            >= CFG.patience
        ):
            break

    if best_state is None:

        raise RuntimeError(
            "No MLP checkpoint selected"
        )

    model.load_state_dict(
        best_state
    )

    return (
        model,
        best_epoch,
        history,
    )


def fit_final(
    X: np.ndarray,
    Y: np.ndarray,
    device: torch.device,
    epochs: int,
    seed: int,
) -> DeterministicMLP:

    set_seed(
        seed
    )

    model = DeterministicMLP().to(
        device
    )

    optimizer = torch.optim.AdamW(
        model.parameters(),
        lr=CFG.learning_rate,
        weight_decay=CFG.weight_decay,
    )

    dataset = TensorDataset(
        torch.from_numpy(
            X
        ),
        torch.from_numpy(
            Y
        ),
    )

    loader = DataLoader(
        dataset,
        batch_size=min(
            CFG.batch_size,
            len(dataset),
        ),
        shuffle=True,
        generator=torch.Generator().manual_seed(
            seed
        ),
    )

    for _ in range(
        epochs
    ):

        model.train()

        for xb, yb in loader:

            xb = xb.to(
                device
            )

            yb = yb.to(
                device
            )

            optimizer.zero_grad(
                set_to_none=True
            )

            pred = model(
                xb
            )

            loss = mlp_loss(
                pred,
                yb,
            )

            loss.backward()

            torch.nn.utils.clip_grad_norm_(
                model.parameters(),
                max_norm=5.0,
            )

            optimizer.step()

    return model


# ============================================================
# Prediction
# ============================================================

def predict(
    model: DeterministicMLP,
    X: np.ndarray,
    device: torch.device,
) -> np.ndarray:

    model.eval()

    output = []

    with torch.no_grad():

        for start in range(
            0,
            len(X),
            CFG.batch_size,
        ):

            xb = torch.from_numpy(
                X[
                    start:
                    start + CFG.batch_size
                ]
            ).to(
                device
            )

            pred = model(
                xb
            )

            output.append(
                pred.cpu().numpy()
            )

    return np.concatenate(
        output
    ).astype(
        np.float32
    )


# ============================================================
# Validation / dependency check
# ============================================================

def validate(
    cohort: str,
) -> None:

    meta = CVAE.cvae_rows(
        cohort
    )

    CVAE.verify_fold_agreement(
        cohort,
        meta,
    )

    gene_matrix = (
        CVAE.load_gene_matrix(
            cohort,
            meta,
        )
    )

    fold_assets = {}

    for fold in range(
        1,
        N_FOLDS + 1,
    ):

        gene_model, gene_scaler = (
            CVAE.gene_checkpoint(
                cohort,
                fold,
            )
        )

        mri_model = (
            CVAE.mri_checkpoint(
                cohort,
                fold,
            )
        )

        fold_assets[str(fold)] = {
            "gene_model":
                str(
                    gene_model.relative_to(
                        ROOT
                    )
                ),

            "gene_scaler":
                str(
                    gene_scaler.relative_to(
                        ROOT
                    )
                ),

            "mri_model":
                str(
                    mri_model.relative_to(
                        ROOT
                    )
                ),
        }

    report = {
        "cohort": cohort,

        "development_n":
            int(
                len(meta)
            ),

        "AD":
            int(
                meta[
                    "diagnosis_binary"
                ].sum()
            ),

        "CN":
            int(
                (
                    meta[
                        "diagnosis_binary"
                    ]
                    == 0
                ).sum()
            ),

        "gene_matrix_shape":
            list(
                gene_matrix.shape
            ),

        "fold_counts": {
            str(k): int(v)
            for k, v
            in meta[
                "cv_fold"
            ]
            .value_counts()
            .sort_index()
            .items()
        },

        "fold_agreement_verified":
            True,

        "branch_assets":
            fold_assets,

        "locked_test_accessed":
            False,
    }

    print(
        json.dumps(
            report,
            indent=2,
            sort_keys=True,
        )
    )


# ============================================================
# Five-fold development
# ============================================================

def develop(
    cohort: str,
    device: torch.device,
) -> None:

    d = dirs_for(
        cohort
    )

    logger = setup_logger(
        cohort,
        "development",
    )

    complete_marker = (
        d["QC"]
        / "development_complete.json"
    )

    if complete_marker.exists():

        raise FileExistsError(
            "Completed MLP development "
            f"already exists: "
            f"{complete_marker}"
        )

    logger.info(
        "%s deterministic MLP "
        "development started on %s",
        cohort,
        device,
    )

    logger.info(
        "Using same paired CVAE folds "
        "and matched frozen branch "
        "encoders as IG-CVAE."
    )

    logger.info(
        "Locked-test participants are "
        "not accessed during development."
    )

    meta = CVAE.cvae_rows(
        cohort
    )

    CVAE.verify_fold_agreement(
        cohort,
        meta,
    )

    gene_matrix = (
        CVAE.load_gene_matrix(
            cohort,
            meta,
        )
    )

    folds = (
        meta[
            "cv_fold"
        ]
        .astype(int)
        .to_numpy()
    )

    y_all = (
        meta[
            "diagnosis_binary"
        ]
        .astype(int)
        .to_numpy()
    )

    fold_rows = []
    oof_rows = []
    provenance_by_fold = {}
    best_epochs = []

    for fold in range(
        1,
        N_FOLDS + 1,
    ):

        logger.info(
            "Starting fold %d",
            fold,
        )

        train_idx = (
            folds != fold
        )

        val_idx = (
            folds == fold
        )

        # ----------------------------------------------------
        # Exactly the same fold-specific frozen branch
        # representations used by the CVAE.
        # ----------------------------------------------------

        (
            mfv_raw,
            gfv_raw,
            provenance,
        ) = CVAE.fold_features(
            cohort,
            fold,
            meta,
            gene_matrix,
            device,
        )

        provenance_by_fold[
            str(fold)
        ] = provenance

        # ----------------------------------------------------
        # Leakage-controlled fold-specific scaling
        # ----------------------------------------------------

        mfv_scaler = (
            StandardScaler()
        )

        gfv_scaler = (
            StandardScaler()
        )

        X_train = (
            mfv_scaler
            .fit_transform(
                mfv_raw[
                    train_idx
                ]
            )
            .astype(
                np.float32
            )
        )

        Y_train = (
            gfv_scaler
            .fit_transform(
                gfv_raw[
                    train_idx
                ]
            )
            .astype(
                np.float32
            )
        )

        X_val = (
            mfv_scaler
            .transform(
                mfv_raw[
                    val_idx
                ]
            )
            .astype(
                np.float32
            )
        )

        Y_val = (
            gfv_scaler
            .transform(
                gfv_raw[
                    val_idx
                ]
            )
            .astype(
                np.float32
            )
        )

        # ----------------------------------------------------
        # Fit MLP
        # ----------------------------------------------------

        (
            model,
            best_epoch,
            history,
        ) = train_model(
            X_train,
            Y_train,
            X_val,
            Y_val,
            device,
            CFG.seed + fold,
        )

        best_epochs.append(
            int(
                best_epoch
            )
        )

        # ----------------------------------------------------
        # MRI-only generation
        # ----------------------------------------------------

        generated_scaled = (
            predict(
                model,
                X_val,
                device,
            )
        )

        generated_raw = (
            gfv_scaler
            .inverse_transform(
                generated_scaled
            )
            .astype(
                np.float32
            )
        )

        # ----------------------------------------------------
        # Reconstruction metrics
        # ----------------------------------------------------

        recon = (
            CVAE.reconstruction_metrics(
                gfv_raw[
                    val_idx
                ],
                generated_raw,
            )
        )

        # ----------------------------------------------------
        # Same downstream AD/CN evaluation as CVAE.
        #
        # Train logistic regression on TRUE training GFVs.
        # Evaluate generated validation GFVs.
        # ----------------------------------------------------

        classifier = (
            LogisticRegression(
                C=1.0,
                class_weight="balanced",
                max_iter=2000,
                random_state=(
                    CFG.seed
                    + fold
                ),
            )
        )

        classifier.fit(
            Y_train,
            y_all[
                train_idx
            ],
        )

        probability = (
            classifier
            .predict_proba(
                generated_scaled
            )[:, 1]
        )

        classification = (
            CVAE.test_metrics(
                y_all[
                    val_idx
                ],
                probability,
            )
        )

        # ----------------------------------------------------
        # Fold metrics
        # ----------------------------------------------------

        record = {
            "fold":
                int(fold),

            "n_train":
                int(
                    train_idx.sum()
                ),

            "n_validation":
                int(
                    val_idx.sum()
                ),

            "best_epoch":
                int(
                    best_epoch
                ),

            **recon,
            **classification,
        }

        fold_rows.append(
            record
        )

        # ----------------------------------------------------
        # Save OOF generated features
        # ----------------------------------------------------

        rows = (
            meta.loc[
                val_idx,
                [
                    "subject_id",
                    "diagnosis_binary",
                    "cv_fold",
                ],
            ]
            .reset_index(
                drop=True
            )
        )

        rows[
            "probability_AD_from_generated_GFV"
        ] = probability

        generated_frame = (
            pd.DataFrame(
                generated_raw,
                columns=[
                    f"MLP_GFV{i+1:03d}"
                    for i
                    in range(
                        GFV_DIM
                    )
                ],
            )
        )

        rows = pd.concat(
            [
                rows,
                generated_frame,
            ],
            axis=1,
        )

        oof_rows.append(
            rows
        )

        # ----------------------------------------------------
        # Save fold model checkpoint
        # ----------------------------------------------------

        torch.save(
            {
                "state_dict":
                    cpu_state_dict(
                        model
                    ),

                "fold":
                    int(fold),

                "seed":
                    int(
                        CFG.seed
                        + fold
                    ),

                "best_epoch":
                    int(
                        best_epoch
                    ),

                "config":
                    asdict(
                        CFG
                    ),

                "mfv_mean":
                    mfv_scaler.mean_,

                "mfv_scale":
                    mfv_scaler.scale_,

                "gfv_mean":
                    gfv_scaler.mean_,

                "gfv_scale":
                    gfv_scaler.scale_,

                "logistic_coef":
                    classifier.coef_,

                "logistic_intercept":
                    classifier.intercept_,

                "branch_provenance":
                    provenance,
            },

            d["Model"]
            / f"fold_{fold}_best.pt",
        )

        pd.DataFrame(
            history
        ).to_csv(
            d["Model"]
            / f"fold_{fold}_history.csv",
            index=False,
        )

        logger.info(
            "fold=%d epoch=%d "
            "MSE=%.4f MAE=%.4f "
            "Cos=%.4f R2=%.4f "
            "AUC=%.4f AP=%.4f",
            fold,
            best_epoch,
            recon["mse"],
            recon["mae"],
            recon[
                "cosine_similarity"
            ],
            recon["r2"],
            classification[
                "roc_auc"
            ],
            classification[
                "pr_auc"
            ],
        )

        del model

        if torch.cuda.is_available():
            torch.cuda.empty_cache()

    # ========================================================
    # Save pooled development outputs
    # ========================================================

    fold_df = pd.DataFrame(
        fold_rows
    )

    fold_df.to_csv(
        d["Table"]
        / "fold_metrics.csv",
        index=False,
    )

    oof = (
        pd.concat(
            oof_rows,
            ignore_index=True,
        )
        .sort_values(
            "subject_id"
        )
        .reset_index(
            drop=True
        )
    )

    if (
        len(oof)
        != len(meta)
    ):

        raise RuntimeError(
            f"{cohort}: OOF count mismatch"
        )

    if (
        oof[
            "subject_id"
        ].duplicated().any()
    ):

        raise RuntimeError(
            f"{cohort}: duplicate OOF subjects"
        )

    oof.to_csv(
        d["Intermediate"]
        / "development_OOF_generated_GFV128.csv",
        index=False,
    )

    metric_names = [
        "mse",
        "mae",
        "cosine_similarity",
        "r2",
        "roc_auc",
        "pr_auc",
        "accuracy",
        "balanced_accuracy",
        "sensitivity",
        "precision",
        "f1",
        "mcc",
        "brier",
    ]

    summary = {
        "cohort":
            cohort,

        "baseline":
            "deterministic_MLP_MFV128_to_GFV128",

        "n":
            int(
                len(meta)
            ),

        "config":
            asdict(
                CFG
            ),

        "best_epochs":
            best_epochs,

        "branch_provenance":
            provenance_by_fold,

        "locked_test_accessed":
            False,
    }

    for name in metric_names:

        if name not in fold_df.columns:
            continue

        summary[
            f"{name}_mean"
        ] = float(
            fold_df[
                name
            ].mean()
        )

        summary[
            f"{name}_sd"
        ] = float(
            fold_df[
                name
            ].std(
                ddof=1
            )
        )

    save_json(
        d["Table"]
        / "fivefold_summary.json",
        summary,
    )

    marker = {
        "status":
            "development_complete",

        "cohort":
            cohort,

        "baseline":
            "deterministic_MLP_MFV128_to_GFV128",

        "run_tag":
            RUN_TAG,

        "n":
            int(
                len(meta)
            ),

        "best_epochs":
            best_epochs,

        "locked_test_accessed":
            False,
    }

    save_json(
        complete_marker,
        marker,
    )

    logger.info(
        "%s development complete.",
        cohort,
    )

    logger.info(
        "CV mean: "
        "MSE=%.4f ± %.4f | "
        "MAE=%.4f ± %.4f | "
        "Cos=%.4f ± %.4f | "
        "R2=%.4f ± %.4f | "
        "AUC=%.4f ± %.4f | "
        "AP=%.4f ± %.4f",
        summary["mse_mean"],
        summary["mse_sd"],
        summary["mae_mean"],
        summary["mae_sd"],
        summary[
            "cosine_similarity_mean"
        ],
        summary[
            "cosine_similarity_sd"
        ],
        summary["r2_mean"],
        summary["r2_sd"],
        summary[
            "roc_auc_mean"
        ],
        summary[
            "roc_auc_sd"
        ],
        summary[
            "pr_auc_mean"
        ],
        summary[
            "pr_auc_sd"
        ],
    )


# ============================================================
# Locked independent test
# ============================================================

def locked_test(
    cohort: str,
    device: torch.device,
    bootstrap_replicates: int,
    allow_locked_test: bool,
) -> None:

    if not allow_locked_test:

        raise PermissionError(
            "TEST requires "
            "--allow-locked-test "
            "after development results "
            "have been reviewed and frozen."
        )

    d = dirs_for(
        cohort
    )

    logger = setup_logger(
        cohort,
        "locked_test",
    )

    development_marker = (
        d["QC"]
        / "development_complete.json"
    )

    if not development_marker.exists():

        raise FileNotFoundError(
            "MLP development must be "
            "completed before TEST: "
            f"{development_marker}"
        )

    marker = (
        d["QC"]
        / "locked_test_opened_once.json"
    )

    if marker.exists():

        raise RuntimeError(
            f"{cohort}: MLP locked test "
            "has already been opened: "
            f"{marker}"
        )

    # --------------------------------------------------------
    # Determine final training duration from DEVELOPMENT only.
    # --------------------------------------------------------

    fold_df = pd.read_csv(
        d["Table"]
        / "fold_metrics.csv"
    )

    final_epochs = max(
        1,
        int(
            round(
                float(
                    fold_df[
                        "best_epoch"
                    ].mean()
                )
            )
        ),
    )

    logger.info(
        "%s final MLP epochs = %d "
        "(mean CV best epoch).",
        cohort,
        final_epochs,
    )

    # --------------------------------------------------------
    # Load development paired subjects using CURRENT
    # final frozen branch models.
    # --------------------------------------------------------

    dev_meta = (
        CVAE.cvae_rows(
            cohort
        )
    )

    (
        dev_mfv_raw,
        dev_gfv_raw,
        dev_provenance,
    ) = CVAE.final_features(
        cohort,
        dev_meta,
        device,
    )

    y_dev = (
        dev_meta[
            "diagnosis_binary"
        ]
        .astype(int)
        .to_numpy()
    )

    # --------------------------------------------------------
    # Fit development-only scalers
    # --------------------------------------------------------

    mfv_scaler = (
        StandardScaler()
    )

    gfv_scaler = (
        StandardScaler()
    )

    dev_mfv = (
        mfv_scaler
        .fit_transform(
            dev_mfv_raw
        )
        .astype(
            np.float32
        )
    )

    dev_gfv = (
        gfv_scaler
        .fit_transform(
            dev_gfv_raw
        )
        .astype(
            np.float32
        )
    )

    # --------------------------------------------------------
    # Final deterministic MLP
    # --------------------------------------------------------

    final_seed = (
        CFG.seed
        + 1000
    )

    model = fit_final(
        dev_mfv,
        dev_gfv,
        device,
        final_epochs,
        final_seed,
    )

    # --------------------------------------------------------
    # Same gene-feature AD/CN classifier used in CVAE logic.
    # --------------------------------------------------------

    classifier = (
        LogisticRegression(
            C=1.0,
            class_weight="balanced",
            max_iter=2000,
            random_state=(
                final_seed
            ),
        )
    )

    classifier.fit(
        dev_gfv,
        y_dev,
    )

    # --------------------------------------------------------
    # Mark locked test BEFORE loading locked-test participants.
    # --------------------------------------------------------

    save_json(
        marker,
        {
            "cohort":
                cohort,

            "status":
                "locked_test_opening_authorized",

            "final_epochs_from_mean_CV_best_epoch":
                int(
                    final_epochs
                ),

            "bootstrap_replicates":
                int(
                    bootstrap_replicates
                ),
        },
    )

    logger.warning(
        "Opening %s MLP locked test "
        "exactly once.",
        cohort,
    )

    # --------------------------------------------------------
    # Locked-test data access starts here.
    # --------------------------------------------------------

    test_meta = (
        CVAE.locked_rows(
            cohort
        )
    )

    (
        test_mfv_raw,
        test_gfv_raw,
        test_provenance,
    ) = CVAE.final_features(
        cohort,
        test_meta,
        device,
    )

    test_mfv = (
        mfv_scaler
        .transform(
            test_mfv_raw
        )
        .astype(
            np.float32
        )
    )

    # --------------------------------------------------------
    # MRI -> MLP-generated GFV
    # --------------------------------------------------------

    generated_scaled = (
        predict(
            model,
            test_mfv,
            device,
        )
    )

    generated_raw = (
        gfv_scaler
        .inverse_transform(
            generated_scaled
        )
        .astype(
            np.float32
        )
    )

    # --------------------------------------------------------
    # Classification from generated GFVs
    # --------------------------------------------------------

    probability = (
        classifier
        .predict_proba(
            generated_scaled
        )[:, 1]
    )

    y = (
        test_meta[
            "diagnosis_binary"
        ]
        .astype(int)
        .to_numpy()
    )

    # --------------------------------------------------------
    # Evaluation
    # --------------------------------------------------------

    reconstruction = (
        CVAE.reconstruction_metrics(
            test_gfv_raw,
            generated_raw,
        )
    )

    classification = (
        CVAE.test_metrics(
            y,
            probability,
        )
    )

    ci = (
        CVAE.stratified_bootstrap_ci(
            y,
            probability,
            test_gfv_raw,
            generated_raw,
            bootstrap_replicates,
        )
    )

    # --------------------------------------------------------
    # Save participant-level output
    # --------------------------------------------------------

    rows = (
        test_meta[
            [
                "subject_id",
                "diagnosis_binary",
            ]
        ]
        .copy()
        .reset_index(
            drop=True
        )
    )

    rows[
        "probability_AD_from_generated_GFV"
    ] = probability

    observed_frame = (
        pd.DataFrame(
            test_gfv_raw,
            columns=[
                f"GFV{i+1:03d}"
                for i
                in range(
                    GFV_DIM
                )
            ],
        )
    )

    generated_frame = (
        pd.DataFrame(
            generated_raw,
            columns=[
                f"MLP_GFV{i+1:03d}"
                for i
                in range(
                    GFV_DIM
                )
            ],
        )
    )

    rows = pd.concat(
        [
            rows,
            observed_frame,
            generated_frame,
        ],
        axis=1,
    )

    rows.to_csv(
        d["Intermediate"]
        / "locked_test_generated_GFV128.csv",
        index=False,
    )

    # --------------------------------------------------------
    # Save final MLP model
    # --------------------------------------------------------

    torch.save(
        {
            "state_dict":
                cpu_state_dict(
                    model
                ),

            "final_epochs":
                int(
                    final_epochs
                ),

            "seed":
                int(
                    final_seed
                ),

            "config":
                asdict(
                    CFG
                ),

            "mfv_mean":
                mfv_scaler.mean_,

            "mfv_scale":
                mfv_scaler.scale_,

            "gfv_mean":
                gfv_scaler.mean_,

            "gfv_scale":
                gfv_scaler.scale_,

            "logistic_coef":
                classifier.coef_,

            "logistic_intercept":
                classifier.intercept_,

            "branch_provenance":
                dev_provenance,
        },

        d["Model"]
        / "final_development_fit.pt",
    )

    # --------------------------------------------------------
    # Save final report
    # --------------------------------------------------------

    report = {
        "cohort":
            cohort,

        "status":
            "locked_test_complete",

        "baseline":
            "deterministic_MLP_MFV128_to_GFV128",

        "run_tag":
            RUN_TAG,

        "n":
            int(
                len(y)
            ),

        "AD":
            int(
                (
                    y == 1
                ).sum()
            ),

        "CN":
            int(
                (
                    y == 0
                ).sum()
            ),

        "final_epochs_from_mean_CV_best_epoch":
            int(
                final_epochs
            ),

        "config":
            asdict(
                CFG
            ),

        "reconstruction":
            reconstruction,

        "classification":
            classification,

        "stratified_bootstrap_95ci":
            ci,

        "bootstrap_replicates":
            int(
                bootstrap_replicates
            ),

        "branch_provenance_development":
            dev_provenance,

        "branch_provenance_locked_test":
            test_provenance,
    }

    save_json(
        d["Table"]
        / "locked_test_results.json",
        report,
    )

    # Replace opening marker by completed report.
    save_json(
        marker,
        report,
    )

    logger.info(
        "%s TEST complete: "
        "MSE=%.4f | "
        "MAE=%.4f | "
        "Cos=%.4f | "
        "R2=%.4f | "
        "AUC=%.4f "
        "[%.4f, %.4f] | "
        "AP=%.4f "
        "[%.4f, %.4f]",
        cohort,
        reconstruction[
            "mse"
        ],
        reconstruction[
            "mae"
        ],
        reconstruction[
            "cosine_similarity"
        ],
        reconstruction[
            "r2"
        ],
        classification[
            "roc_auc"
        ],
        ci[
            "roc_auc"
        ][0],
        ci[
            "roc_auc"
        ][1],
        classification[
            "pr_auc"
        ],
        ci[
            "pr_auc"
        ][0],
        ci[
            "pr_auc"
        ][1],
    )


# ============================================================
# Main
# ============================================================

def main() -> None:

    parser = argparse.ArgumentParser(
        description=(
            "Deterministic MLP baseline: "
            "MRI MFV128 -> generated gene GFV128 "
            "using the same frozen branch representations "
            "and fixed folds as IG-CVAE."
        )
    )

    parser.add_argument(
        "command",
        choices=[
            "VALIDATE",
            "DEVELOP",
            "TEST",
        ],
    )

    parser.add_argument(
        "--cohort",
        choices=[
            "ADNI",
            "AddNeuroMed",
            "ALL",
        ],
        default="ALL",
    )

    parser.add_argument(
        "--device",
        choices=[
            "auto",
            "cpu",
            "cuda",
        ],
        default="auto",
    )

    parser.add_argument(
        "--allow-locked-test",
        action="store_true",
    )

    parser.add_argument(
        "--bootstrap-replicates",
        type=int,
        default=2000,
    )

    args = parser.parse_args()

    if args.device == "auto":

        device = torch.device(
            "cuda"
            if torch.cuda.is_available()
            else "cpu"
        )

    else:

        device = torch.device(
            args.device
        )

    if (
        device.type == "cuda"
        and not torch.cuda.is_available()
    ):

        raise RuntimeError(
            "CUDA requested "
            "but CUDA is unavailable."
        )

    cohorts = (
        [
            "ADNI",
            "AddNeuroMed",
        ]
        if args.cohort
        == "ALL"
        else [
            args.cohort
        ]
    )

    for cohort in cohorts:

        if (
            args.command
            == "VALIDATE"
        ):

            validate(
                cohort
            )

        elif (
            args.command
            == "DEVELOP"
        ):

            develop(
                cohort,
                device,
            )

        else:

            locked_test(
                cohort,
                device,
                args.bootstrap_replicates,
                args.allow_locked_test,
            )


if __name__ == "__main__":
    main()
