#!/usr/bin/env python3
"""
Covariate-adjusted IG-CVAE sensitivity analysis.

"""

from __future__ import annotations

import argparse
import importlib.util
import json
import logging
import os
import sys
from dataclasses import asdict
from pathlib import Path

import numpy as np
import pandas as pd
import torch
import warnings

from sklearn.exceptions import InconsistentVersionWarning
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler

# The pretrained branch StandardScaler objects were saved with an older
# scikit-learn version. Suppress only this known persistence-version warning.
warnings.filterwarnings(
    "ignore",
    category=InconsistentVersionWarning,
)


ROOT = Path(__file__).resolve().parents[1]
CVAE_SCRIPT = ROOT / "scripts" / "06_run_cvae.py"

RUN_TAG = "covadj_age_sex_apoe4_v1"
MODEL_GROUP = "CVAE_covadj"

ADNI_METADATA = (
    ROOT.parent
    / "data"
    / "ADNI"
    / "Paired_data_metadata.csv"
)

ANM_METADATA = (
    ROOT.parent
    / "data"
    / "AdNeuroMed"
    / "Paired_data_metadata.csv"
)


# ============================================================
# Load ORIGINAL current CVAE
# ============================================================

def load_current_cvae():
    name = "igcvae_primary_for_covadj"

    spec = importlib.util.spec_from_file_location(
        name,
        CVAE_SCRIPT,
    )

    if spec is None or spec.loader is None:
        raise RuntimeError(
            f"Cannot import {CVAE_SCRIPT}"
        )

    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)

    return mod


CV = load_current_cvae()


# ============================================================
# Separate covariate-adjusted output tree
# ============================================================

def outdir(kind: str, cohort: str) -> Path:
    p = (
        ROOT
        / "Output"
        / kind
        / MODEL_GROUP
        / RUN_TAG
        / cohort
    )

    p.mkdir(
        parents=True,
        exist_ok=True,
    )

    return p


def temp_mri_cache(cohort: str) -> Path:
    """
    Original MRI feature extraction expects a preprocessing cache.

    Keep this technical cache outside the project Output directory so
    no adjusted intermediate dataset is written into the project.
    """
    base = Path(
        os.environ.get(
            "TMPDIR",
            "/tmp",
        )
    )

    p = (
        base
        / f"igcvae_covadj_mri_cache_{os.getuid()}"
        / cohort
    )

    p.mkdir(
        parents=True,
        exist_ok=True,
    )

    return p


def _covadj_cv_outdir(
    kind: str,
    cohort: str,
) -> Path:

    if kind == "Intermediate":
        return temp_mri_cache(
            cohort
        )

    return outdir(
        kind,
        cohort,
    )


# The original fold/final MRI extraction uses CV.outdir only
# to choose its preprocessing cache. Redirect it away from
# the primary model output tree.
CV.outdir = _covadj_cv_outdir


# ============================================================
# General utilities
# ============================================================

def save_json(
    path: Path,
    obj,
) -> None:

    def convert(x):

        if isinstance(
            x,
            np.integer,
        ):
            return int(x)

        if isinstance(
            x,
            np.floating,
        ):
            return float(x)

        if isinstance(
            x,
            np.ndarray,
        ):
            return x.tolist()

        if isinstance(
            x,
            Path,
        ):
            return str(x)

        raise TypeError(
            type(x).__name__
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

    log = logging.getLogger(
        f"covadj_cvae_{cohort}_{phase}"
    )

    log.handlers.clear()
    log.setLevel(
        logging.INFO
    )

    fmt = logging.Formatter(
        "%(asctime)s | %(levelname)s | %(message)s"
    )

    sh = logging.StreamHandler(
        sys.stdout
    )

    sh.setFormatter(
        fmt
    )

    log.addHandler(
        sh
    )

    fh = logging.FileHandler(
        outdir(
            "Log",
            cohort,
        )
        / f"{phase}.log",
        mode="a",
    )

    fh.setFormatter(
        fmt
    )

    log.addHandler(
        fh
    )

    return log


def cpu_state_dict(
    model,
):

    return {
        k:
            v.detach()
            .cpu()
            .clone()

        for k, v
        in model.state_dict().items()
    }


def scale_transform(
    X,
    mean,
    scale,
):

    X = np.asarray(
        X,
        dtype=np.float64,
    )

    mean = np.asarray(
        mean,
        dtype=np.float64,
    )

    scale = np.asarray(
        scale,
        dtype=np.float64,
    )

    return (
        (
            X
            - mean[
                None,
                :
            ]
        )
        / scale[
            None,
            :
        ]
    ).astype(
        np.float32
    )


def scale_inverse(
    X,
    mean,
    scale,
):

    X = np.asarray(
        X,
        dtype=np.float64,
    )

    mean = np.asarray(
        mean,
        dtype=np.float64,
    )

    scale = np.asarray(
        scale,
        dtype=np.float64,
    )

    return (
        X
        * scale[
            None,
            :
        ]
        + mean[
            None,
            :
        ]
    ).astype(
        np.float32
    )


# ============================================================
# Covariates
# age + sex + APOE4 only
# NO batch covariate
# ============================================================

def normalize_sex(
    value,
):

    s = (
        str(value)
        .strip()
        .upper()
    )

    if s in {
        "M",
        "MALE",
    }:
        return 1.0

    if s in {
        "F",
        "FEMALE",
    }:
        return 0.0

    raise ValueError(
        f"Unsupported sex value: {value!r}"
    )


def normalize_apoe4(
    value,
):

    s = (
        str(value)
        .strip()
        .upper()
    )

    if s in {
        "0",
        "0.0",
    }:
        return 0.0

    if s in {
        "1",
        "1.0",
    }:
        return 1.0

    if s in {
        "2",
        "2.0",
    }:
        return 2.0

    s = "".join(
        ch
        for ch in s
        if ch.isalnum()
    )

    mapping = {
        "E2E2": 0.0,
        "E2E3": 0.0,
        "E3E2": 0.0,
        "E3E3": 0.0,

        "E2E4": 1.0,
        "E4E2": 1.0,
        "E3E4": 1.0,
        "E4E3": 1.0,

        "E4E4": 2.0,
    }

    if s not in mapping:
        raise ValueError(
            f"Unsupported APOE value: {value!r}"
        )

    return mapping[s]


def covariate_matrix(
    cohort: str,
    subject_ids,
):

    subject_ids = [
        str(x).strip()
        for x in subject_ids
    ]

    # --------------------------------------------------------
    # ADNI
    # --------------------------------------------------------

    if cohort == "ADNI":

        path = ADNI_METADATA

        if not path.is_file():
            raise FileNotFoundError(
                path
            )

        df = pd.read_csv(
            path,
            dtype=str,
        )

        required = {
            "subject_id",
            "age",
            "sex",
            "apoe4_allele_count",
        }

        missing_cols = sorted(
            required
            - set(
                df.columns
            )
        )

        if missing_cols:

            raise RuntimeError(
                "ADNI metadata missing columns: "
                f"{missing_cols}"
            )

        df["_sid"] = (
            df[
                "subject_id"
            ]
            .astype(str)
            .str.strip()
        )

        if (
            df[
                "_sid"
            ]
            .duplicated()
            .any()
        ):

            raise RuntimeError(
                "Duplicate ADNI paired metadata IDs"
            )

        df = df.set_index(
            "_sid"
        )

        missing = sorted(
            set(
                subject_ids
            )
            - set(
                df.index
            )
        )

        if missing:

            raise RuntimeError(
                "ADNI covariates missing for "
                f"{len(missing)} subjects: "
                f"{missing[:10]}"
            )

        q = df.loc[
            subject_ids
        ]

        age = pd.to_numeric(
            q[
                "age"
            ],
            errors="raise",
        ).to_numpy(
            float
        )

        sex = np.asarray(
            [
                normalize_sex(
                    x
                )
                for x in q[
                    "sex"
                ]
            ],
            dtype=float,
        )

        apoe = pd.to_numeric(
            q[
                "apoe4_allele_count"
            ],
            errors="raise",
        ).to_numpy(
            float
        )

    # --------------------------------------------------------
    # AddNeuroMed
    # --------------------------------------------------------

    elif cohort == "AddNeuroMed":

        path = ANM_METADATA

        if not path.is_file():
            raise FileNotFoundError(
                path
            )

        df = pd.read_csv(
            path,
            dtype=str,
        )

        required = {
            "Subject_ID",
            "Age",
            "Sex",
            "APOE",
        }

        missing_cols = sorted(
            required
            - set(
                df.columns
            )
        )

        if missing_cols:

            raise RuntimeError(
                "AddNeuroMed metadata missing columns: "
                f"{missing_cols}"
            )

        df["_sid"] = (
            df[
                "Subject_ID"
            ]
            .astype(str)
            .str.strip()
        )

        if (
            df[
                "_sid"
            ]
            .duplicated()
            .any()
        ):

            raise RuntimeError(
                "Duplicate AddNeuroMed paired metadata IDs"
            )

        df = df.set_index(
            "_sid"
        )

        missing = sorted(
            set(
                subject_ids
            )
            - set(
                df.index
            )
        )

        if missing:

            raise RuntimeError(
                "AddNeuroMed covariates missing for "
                f"{len(missing)} subjects: "
                f"{missing[:10]}"
            )

        q = df.loc[
            subject_ids
        ]

        age = pd.to_numeric(
            q[
                "Age"
            ],
            errors="raise",
        ).to_numpy(
            float
        )

        sex = np.asarray(
            [
                normalize_sex(
                    x
                )
                for x in q[
                    "Sex"
                ]
            ],
            dtype=float,
        )

        apoe = np.asarray(
            [
                normalize_apoe4(
                    x
                )
                for x in q[
                    "APOE"
                ]
            ],
            dtype=float,
        )

    else:

        raise ValueError(
            cohort
        )

    C = np.column_stack(
        [
            age,
            sex,
            apoe,
        ]
    ).astype(
        np.float64
    )

    names = [
        "age",
        "sex_male",
        "apoe4_count",
    ]

    if not np.isfinite(
        C
    ).all():

        raise RuntimeError(
            f"{cohort}: non-finite covariates"
        )

    return (
        C,
        names,
    )


# ============================================================
# Covariate residualization
# ============================================================

def fit_residualizer(
    Y,
    C,
    train_mask,
    covariate_names,
):

    Y = np.asarray(
        Y,
        dtype=np.float64,
    )

    C = np.asarray(
        C,
        dtype=np.float64,
    )

    train_mask = np.asarray(
        train_mask,
        dtype=bool,
    )

    if Y.ndim != 2:
        raise ValueError(
            "Y must be N x features"
        )

    if (
        C.ndim != 2
        or len(C) != len(Y)
    ):

        raise ValueError(
            "Covariate/representation size mismatch"
        )

    if (
        train_mask.shape
        != (
            len(Y),
        )
    ):

        raise ValueError(
            "Invalid train_mask"
        )

    # Standardize covariates using TRAINING subjects only.
    scaler = StandardScaler()

    C_train_scaled = (
        scaler
        .fit_transform(
            C[
                train_mask
            ]
        )
    )

    D_train = np.column_stack(
        [
            np.ones(
                train_mask.sum(),
                dtype=np.float64,
            ),

            C_train_scaled,
        ]
    )

    beta, _, _, _ = (
        np.linalg.lstsq(
            D_train,
            Y[
                train_mask
            ],
            rcond=None,
        )
    )

    reference_mean = (
        Y[
            train_mask
        ]
        .mean(
            axis=0
        )
    )

    return {

        "covariate_names":
            list(
                covariate_names
            ),

        "covariate_mean":
            scaler.mean_
            .astype(
                np.float64
            ),

        "covariate_scale":
            scaler.scale_
            .astype(
                np.float64
            ),

        "beta":
            beta.astype(
                np.float64
            ),

        "reference_mean":
            reference_mean.astype(
                np.float64
            ),

        "design_rank":
            int(
                np.linalg.matrix_rank(
                    D_train
                )
            ),

        "n_train":
            int(
                train_mask.sum()
            ),
    }


def apply_residualizer(
    Y,
    C,
    fit,
):

    Y = np.asarray(
        Y,
        dtype=np.float64,
    )

    C = np.asarray(
        C,
        dtype=np.float64,
    )

    mean = np.asarray(
        fit[
            "covariate_mean"
        ],
        dtype=np.float64,
    )

    scale = np.asarray(
        fit[
            "covariate_scale"
        ],
        dtype=np.float64,
    )

    beta = np.asarray(
        fit[
            "beta"
        ],
        dtype=np.float64,
    )

    reference_mean = np.asarray(
        fit[
            "reference_mean"
        ],
        dtype=np.float64,
    )

    C_scaled = (
        C
        - mean[
            None,
            :
        ]
    ) / scale[
        None,
        :
    ]

    D = np.column_stack(
        [
            np.ones(
                len(
                    C_scaled
                ),
                dtype=np.float64,
            ),

            C_scaled,
        ]
    )

    fitted = (
        D
        @ beta
    )

    adjusted = (
        Y
        - fitted
        + reference_mean[
            None,
            :
        ]
    )

    if not np.isfinite(
        adjusted
    ).all():

        raise RuntimeError(
            "Residualization produced non-finite values"
        )

    return adjusted.astype(
        np.float32
    )


def fit_and_apply_residualizer(
    Y,
    C,
    train_mask,
    covariate_names,
):

    fit = fit_residualizer(
        Y,
        C,
        train_mask,
        covariate_names,
    )

    adjusted = apply_residualizer(
        Y,
        C,
        fit,
    )

    return (
        adjusted,
        fit,
    )


# ============================================================
# DEVELOPMENT
# ============================================================

def develop(
    cohort: str,
    device: torch.device,
) -> None:

    model_dir = outdir(
        "Model",
        cohort,
    )

    table_dir = outdir(
        "Table",
        cohort,
    )

    qc_dir = outdir(
        "QC",
        cohort,
    )

    complete_marker = (
        qc_dir
        / "development_complete.json"
    )

    final_checkpoint = (
        model_dir
        / "final_development_fit.pt"
    )

    if (
        complete_marker.exists()
        or final_checkpoint.exists()
    ):

        raise RuntimeError(
            f"{cohort}: covadj development already exists; "
            "refusing to overwrite."
        )

    log = setup_logger(
        cohort,
        "development",
    )

    # --------------------------------------------------------
    # Same subjects/folds as primary CVAE
    # --------------------------------------------------------

    meta = CV.cvae_rows(
        cohort
    )

    CV.verify_fold_agreement(
        cohort,
        meta,
    )

    gene_matrix = (
        CV.load_gene_matrix(
            cohort,
            meta,
        )
    )

    C, covariate_names = (
        covariate_matrix(
            cohort,
            meta.subject_id.tolist(),
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

    log.info(
        "%s covadj CVAE development started on %s",
        cohort,
        device,
    )

    log.info(
        "Covariates: %s",
        covariate_names,
    )

    log.info(
        "Batch covariate: NOT INCLUDED"
    )

    log.info(
        "Locked test accessed: False"
    )

    fold_rows = []
    best_epochs = []

    # ========================================================
    # Five-fold development
    # ========================================================

    for fold in range(
        1,
        CV.N_FOLDS + 1,
    ):

        tr = (
            folds != fold
        )

        va = (
            folds == fold
        )

        # ----------------------------------------------------
        # Same frozen Gene/MRI representations as primary CVAE
        # ----------------------------------------------------

        (
            mfv_raw,
            gfv_raw,
            branch_info,
        ) = CV.fold_features(
            cohort,
            fold,
            meta,
            gene_matrix,
            device,
        )

        # ----------------------------------------------------
        # Covariate models fit TRAINING ONLY
        # ----------------------------------------------------

        (
            mfv_adj,
            mfv_residualizer,
        ) = fit_and_apply_residualizer(
            mfv_raw,
            C,
            tr,
            covariate_names,
        )

        (
            gfv_adj,
            gfv_residualizer,
        ) = fit_and_apply_residualizer(
            gfv_raw,
            C,
            tr,
            covariate_names,
        )

        # ----------------------------------------------------
        # Same train-only MFV/GFV scaling as primary CVAE
        # ----------------------------------------------------

        mfv_scaler = StandardScaler()
        gfv_scaler = StandardScaler()

        mfv_tr = (
            mfv_scaler
            .fit_transform(
                mfv_adj[
                    tr
                ]
            )
            .astype(
                np.float32
            )
        )

        gfv_tr = (
            gfv_scaler
            .fit_transform(
                gfv_adj[
                    tr
                ]
            )
            .astype(
                np.float32
            )
        )

        mfv_va = (
            mfv_scaler
            .transform(
                mfv_adj[
                    va
                ]
            )
            .astype(
                np.float32
            )
        )

        gfv_va = (
            gfv_scaler
            .transform(
                gfv_adj[
                    va
                ]
            )
            .astype(
                np.float32
            )
        )

        # ----------------------------------------------------
        # ORIGINAL CVAE training
        # ----------------------------------------------------

        (
            model,
            best_epoch,
        ) = CV.fit_fold(
            mfv_tr,
            gfv_tr,
            mfv_va,
            gfv_va,
            device,
            CV.CFG.seed
            + fold,
        )

        generated_scaled = (
            CV.predict_generated(
                model,
                mfv_va,
                device,
            )
        )

        generated = (
            gfv_scaler
            .inverse_transform(
                generated_scaled
            )
            .astype(
                np.float32
            )
        )

        reconstruction = (
            CV.reconstruction_metrics(
                gfv_adj[
                    va
                ],
                generated,
            )
        )

        # ----------------------------------------------------
        # Same downstream classifier as primary CVAE
        # ----------------------------------------------------

        classifier = (
            LogisticRegression(
                C=1.0,
                class_weight="balanced",
                max_iter=2000,
                random_state=(
                    CV.CFG.seed
                    + fold
                ),
            )
        )

        classifier.fit(
            gfv_tr,
            y_all[
                tr
            ],
        )

        probability = (
            classifier
            .predict_proba(
                generated_scaled
            )[
                :,
                1
            ]
        )

        classification = (
            CV.class_metrics(
                y_all[
                    va
                ],
                probability,
            )
        )

        fold_rows.append(
            {

                "fold":
                    int(
                        fold
                    ),

                "n_train":
                    int(
                        tr.sum()
                    ),

                "n_validation":
                    int(
                        va.sum()
                    ),

                "best_epoch":
                    int(
                        best_epoch
                    ),

                "mfv_covariate_design_rank":
                    int(
                        mfv_residualizer[
                            "design_rank"
                        ]
                    ),

                "gfv_covariate_design_rank":
                    int(
                        gfv_residualizer[
                            "design_rank"
                        ]
                    ),

                **reconstruction,
                **classification,
            }
        )

        best_epochs.append(
            int(
                best_epoch
            )
        )

        # ----------------------------------------------------
        # Separate covadj fold checkpoint
        # ----------------------------------------------------

        torch.save(
            {

                "analysis":
                    "covariate_adjusted_CVAE_sensitivity",

                "run_tag":
                    RUN_TAG,

                "cohort":
                    cohort,

                "fold":
                    int(
                        fold
                    ),

                "state_dict":
                    cpu_state_dict(
                        model
                    ),

                "best_epoch":
                    int(
                        best_epoch
                    ),

                "config":
                    asdict(
                        CV.CFG
                    ),

                "covariate_names":
                    covariate_names,

                "batch_covariate_included":
                    False,

                "mfv_residualizer":
                    mfv_residualizer,

                "gfv_residualizer":
                    gfv_residualizer,

                "mfv_mean":
                    mfv_scaler.mean_
                    .astype(
                        np.float64
                    ),

                "mfv_scale":
                    mfv_scaler.scale_
                    .astype(
                        np.float64
                    ),

                "gfv_mean":
                    gfv_scaler.mean_
                    .astype(
                        np.float64
                    ),

                "gfv_scale":
                    gfv_scaler.scale_
                    .astype(
                        np.float64
                    ),

                "logistic_coef":
                    classifier.coef_
                    .astype(
                        np.float64
                    ),

                "logistic_intercept":
                    classifier.intercept_
                    .astype(
                        np.float64
                    ),

                "branch_provenance":
                    branch_info,

                "primary_cvae_script":
                    str(
                        CVAE_SCRIPT
                        .relative_to(
                            ROOT
                        )
                    ),

                "primary_cvae_script_sha256":
                    CV.sha256(
                        CVAE_SCRIPT
                    ),
            },

            model_dir
            / f"fold_{fold}_best.pt",
        )

        log.info(
            "fold=%d epoch=%d "
            "MSE=%.4f MAE=%.4f "
            "Cos=%.4f R2=%.4f "
            "AUC=%.4f AP=%.4f",
            fold,
            best_epoch,
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
            classification[
                "pr_auc"
            ],
        )

        del model

        if torch.cuda.is_available():
            torch.cuda.empty_cache()

    # ========================================================
    # Development fold metrics
    # ========================================================

    metrics = pd.DataFrame(
        fold_rows
    )

    metrics.to_csv(
        table_dir
        / "fold_metrics.csv",
        index=False,
    )

    exclude = {
        "fold",
        "n_train",
        "n_validation",
        "best_epoch",
        "mfv_covariate_design_rank",
        "gfv_covariate_design_rank",
    }

    numeric = [
        c
        for c in metrics.columns
        if (
            c not in exclude
            and pd.api.types
            .is_numeric_dtype(
                metrics[
                    c
                ]
            )
        )
    ]

    summary = {

        "cohort":
            cohort,

        "analysis":
            "covariate_adjusted_CVAE_sensitivity",

        "run_tag":
            RUN_TAG,

        "n_development":
            int(
                len(
                    meta
                )
            ),

        "covariates":
            covariate_names,

        "batch_covariate_included":
            False,

        "best_epochs":
            best_epochs,

        "cvae_config":
            asdict(
                CV.CFG
            ),

        "locked_test_accessed":
            False,
    }

    for name in numeric:

        summary[
            f"{name}_mean"
        ] = float(
            metrics[
                name
            ].mean()
        )

        summary[
            f"{name}_sd"
        ] = float(
            metrics[
                name
            ].std(
                ddof=1
            )
        )

    save_json(
        table_dir
        / "fivefold_mean_sd.json",
        summary,
    )

    # ========================================================
    # Freeze final development-only covadj model
    #
    # TEST WILL NOT TRAIN OR REFIT.
    # ========================================================

    final_epochs = max(
        1,
        int(
            round(
                float(
                    metrics[
                        "best_epoch"
                    ]
                    .mean()
                )
            )
        ),
    )

    log.info(
        "Training final frozen development model "
        "for %d epochs.",
        final_epochs,
    )

    # Same original final Gene/MRI branch checkpoints
    # used by primary CVAE.
    (
        dev_mfv_raw,
        dev_gfv_raw,
        final_branch_info,
    ) = CV.final_features(
        cohort,
        meta,
        device,
    )

    all_development = np.ones(
        len(
            meta
        ),
        dtype=bool,
    )

    (
        dev_mfv_adj,
        final_mfv_residualizer,
    ) = fit_and_apply_residualizer(
        dev_mfv_raw,
        C,
        all_development,
        covariate_names,
    )

    (
        dev_gfv_adj,
        final_gfv_residualizer,
    ) = fit_and_apply_residualizer(
        dev_gfv_raw,
        C,
        all_development,
        covariate_names,
    )

    final_mfv_scaler = (
        StandardScaler()
    )

    final_gfv_scaler = (
        StandardScaler()
    )

    dev_mfv = (
        final_mfv_scaler
        .fit_transform(
            dev_mfv_adj
        )
        .astype(
            np.float32
        )
    )

    dev_gfv = (
        final_gfv_scaler
        .fit_transform(
            dev_gfv_adj
        )
        .astype(
            np.float32
        )
    )

    final_seed = (
        CV.CFG.seed
        + 1000
    )

    final_model = (
        CV.fit_final(
            dev_mfv,
            dev_gfv,
            device,
            final_epochs,
            final_seed,
        )
    )

    final_classifier = (
        LogisticRegression(
            C=1.0,
            class_weight="balanced",
            max_iter=2000,
            random_state=(
                final_seed
            ),
        )
    )

    final_classifier.fit(
        dev_gfv,
        y_all,
    )

    # --------------------------------------------------------
    # Final checkpoint contains everything TEST needs.
    # --------------------------------------------------------

    torch.save(
        {

            "analysis":
                "covariate_adjusted_CVAE_sensitivity",

            "status":
                "frozen_development_model",

            "run_tag":
                RUN_TAG,

            "cohort":
                cohort,

            "state_dict":
                cpu_state_dict(
                    final_model
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
                    CV.CFG
                ),

            "covariate_names":
                covariate_names,

            "batch_covariate_included":
                False,

            "mfv_residualizer":
                final_mfv_residualizer,

            "gfv_residualizer":
                final_gfv_residualizer,

            "mfv_mean":
                final_mfv_scaler.mean_
                .astype(
                    np.float64
                ),

            "mfv_scale":
                final_mfv_scaler.scale_
                .astype(
                    np.float64
                ),

            "gfv_mean":
                final_gfv_scaler.mean_
                .astype(
                    np.float64
                ),

            "gfv_scale":
                final_gfv_scaler.scale_
                .astype(
                    np.float64
                ),

            "logistic_coef":
                final_classifier.coef_
                .astype(
                    np.float64
                ),

            "logistic_intercept":
                final_classifier.intercept_
                .astype(
                    np.float64
                ),

            "branch_provenance_development":
                final_branch_info,

            "primary_cvae_script":
                str(
                    CVAE_SCRIPT
                    .relative_to(
                        ROOT
                    )
                ),

            "primary_cvae_script_sha256":
                CV.sha256(
                    CVAE_SCRIPT
                ),
        },

        final_checkpoint,
    )

    save_json(
        complete_marker,
        {

            "status":
                "development_complete_and_final_model_frozen",

            "cohort":
                cohort,

            "run_tag":
                RUN_TAG,

            "n_development":
                int(
                    len(
                        meta
                    )
                ),

            "covariates":
                covariate_names,

            "batch_covariate_included":
                False,

            "best_epochs":
                best_epochs,

            "final_epochs_from_mean_CV_best_epoch":
                int(
                    final_epochs
                ),

            "final_checkpoint":
                str(
                    final_checkpoint
                    .relative_to(
                        ROOT
                    )
                ),

            "locked_test_accessed":
                False,
        },
    )

    log.info(
        "%s development complete.",
        cohort,
    )

    log.info(
        "Final checkpoint: %s",
        final_checkpoint,
    )

    log.info(
        "Locked test remains unopened."
    )


# ============================================================
# LOCKED TEST
#
# Evaluation only. No fitting or refitting.
# ============================================================

def locked_test(
    cohort: str,
    device: torch.device,
    allow_locked_test: bool,
    bootstrap_replicates: int,
) -> None:

    if not allow_locked_test:

        raise PermissionError(
            "TEST requires --allow-locked-test"
        )

    model_dir = outdir(
        "Model",
        cohort,
    )

    table_dir = outdir(
        "Table",
        cohort,
    )

    qc_dir = outdir(
        "QC",
        cohort,
    )

    final_checkpoint = (
        model_dir
        / "final_development_fit.pt"
    )

    development_marker = (
        qc_dir
        / "development_complete.json"
    )

    test_marker = (
        qc_dir
        / "locked_test_opened_once.json"
    )

    result_path = (
        table_dir
        / "locked_test_metrics.json"
    )

    if not development_marker.exists():

        raise RuntimeError(
            f"{cohort}: run DEVELOP before TEST."
        )

    if not final_checkpoint.is_file():

        raise FileNotFoundError(
            final_checkpoint
        )

    if (
        test_marker.exists()
        or result_path.exists()
    ):

        raise RuntimeError(
            f"{cohort}: locked test already "
            "opened/completed; refusing second access."
        )

    log = setup_logger(
        cohort,
        "locked_test",
    )

    # --------------------------------------------------------
    # Load fully frozen development checkpoint.
    # --------------------------------------------------------

    checkpoint = torch.load(
        final_checkpoint,
        map_location="cpu",
        weights_only=False,
    )

    expected_covariates = [
        "age",
        "sex_male",
        "apoe4_count",
    ]

    if (
        checkpoint[
            "covariate_names"
        ]
        != expected_covariates
    ):

        raise RuntimeError(
            "Unexpected covariates in frozen checkpoint."
        )

    # --------------------------------------------------------
    # Mark test opened BEFORE loading test rows.
    # --------------------------------------------------------

    save_json(
        test_marker,
        {

            "cohort":
                cohort,

            "status":
                "opened_once_before_locked_data_load",

            "run_tag":
                RUN_TAG,

            "final_checkpoint":
                str(
                    final_checkpoint
                    .relative_to(
                        ROOT
                    )
                ),

            "final_epochs":
                int(
                    checkpoint[
                        "final_epochs"
                    ]
                ),

            "covariates":
                checkpoint[
                    "covariate_names"
                ],
        },
    )

    log.warning(
        "Opening %s covadj locked test exactly once.",
        cohort,
    )

    # ========================================================
    # Locked-test access begins here.
    # ========================================================

    test_meta = (
        CV.locked_rows(
            cohort
        )
    )

    (
        test_mfv_raw,
        test_gfv_raw,
        test_branch_info,
    ) = CV.final_features(
        cohort,
        test_meta,
        device,
    )

    (
        C_test,
        covariate_names,
    ) = covariate_matrix(
        cohort,
        test_meta
        .subject_id
        .tolist(),
    )

    if (
        covariate_names
        != checkpoint[
            "covariate_names"
        ]
    ):

        raise RuntimeError(
            "Development/test covariate definitions differ."
        )

    # --------------------------------------------------------
    # Apply development-fitted residualizers unchanged.
    # --------------------------------------------------------

    test_mfv_adj = (
        apply_residualizer(
            test_mfv_raw,
            C_test,
            checkpoint[
                "mfv_residualizer"
            ],
        )
    )

    test_gfv_adj = (
        apply_residualizer(
            test_gfv_raw,
            C_test,
            checkpoint[
                "gfv_residualizer"
            ],
        )
    )

    # --------------------------------------------------------
    # Development-fitted MFV scaling
    # --------------------------------------------------------

    test_mfv = (
        scale_transform(
            test_mfv_adj,
            checkpoint[
                "mfv_mean"
            ],
            checkpoint[
                "mfv_scale"
            ],
        )
    )

    # --------------------------------------------------------
    # Load frozen CVAE
    # --------------------------------------------------------

    model = (
        CV.ConditionalVAE()
        .to(
            device
        )
    )

    model.load_state_dict(
        checkpoint[
            "state_dict"
        ]
    )

    generated_scaled = (
        CV.predict_generated(
            model,
            test_mfv,
            device,
        )
    )

    generated = (
        scale_inverse(
            generated_scaled,
            checkpoint[
                "gfv_mean"
            ],
            checkpoint[
                "gfv_scale"
            ],
        )
    )

    # --------------------------------------------------------
    # Frozen development logistic-regression classifier.
    # Reproduce predict_proba analytically.
    # --------------------------------------------------------

    coef = np.asarray(
        checkpoint[
            "logistic_coef"
        ],
        dtype=np.float64,
    )

    intercept = np.asarray(
        checkpoint[
            "logistic_intercept"
        ],
        dtype=np.float64,
    )

    logits = (
        generated_scaled
        .astype(
            np.float64
        )
        @ coef[
            0
        ]
        + intercept[
            0
        ]
    )

    probability = (
        1.0
        / (
            1.0
            + np.exp(
                -np.clip(
                    logits,
                    -50.0,
                    50.0,
                )
            )
        )
    )

    y = (
        test_meta[
            "diagnosis_binary"
        ]
        .astype(int)
        .to_numpy()
    )

    # --------------------------------------------------------
    # Same metrics as primary CVAE.
    # Reconstruction is against ADJUSTED GFV target.
    # --------------------------------------------------------

    reconstruction = (
        CV.reconstruction_metrics(
            test_gfv_adj,
            generated,
        )
    )

    classification = (
        CV.test_metrics(
            y,
            probability,
        )
    )

    ci = (
        CV.stratified_bootstrap_ci(
            y,
            probability,
            test_gfv_adj,
            generated,
            bootstrap_replicates,
        )
    )

    result = {

        "cohort":
            cohort,

        "status":
            "locked_test_complete",

        "analysis":
            "covariate_adjusted_CVAE_sensitivity",

        "run_tag":
            RUN_TAG,

        "n":
            int(
                len(
                    y
                )
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

        "covariates":
            checkpoint[
                "covariate_names"
            ],

        "batch_covariate_included":
            False,

        "final_epochs_from_mean_CV_best_epoch":
            int(
                checkpoint[
                    "final_epochs"
                ]
            ),

        "bootstrap_replicates":
            int(
                bootstrap_replicates
            ),

        "reconstruction":
            reconstruction,

        "classification":
            classification,

        "stratified_bootstrap_95ci":
            ci,

        "branch_provenance_development":
            checkpoint[
                "branch_provenance_development"
            ],

        "branch_provenance_locked_test":
            test_branch_info,

        "final_checkpoint":
            str(
                final_checkpoint
                .relative_to(
                    ROOT
                )
            ),

        "primary_cvae_script":
            checkpoint[
                "primary_cvae_script"
            ],

        "primary_cvae_script_sha256":
            checkpoint[
                "primary_cvae_script_sha256"
            ],
    }

    save_json(
        result_path,
        result,
    )

    # Replace opening marker with completed result.
    save_json(
        test_marker,
        result,
    )

    log.info(
        "%s TEST complete: "
        "MSE=%.4f MAE=%.4f "
        "Cos=%.4f R2=%.4f "
        "AUC=%.4f AP=%.4f",
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
        classification[
            "pr_auc"
        ],
    )


# ============================================================
# Validation
# ============================================================

def validate(
    cohort: str,
) -> None:

    meta = (
        CV.cvae_rows(
            cohort
        )
    )

    CV.verify_fold_agreement(
        cohort,
        meta,
    )

    gene_matrix = (
        CV.load_gene_matrix(
            cohort,
            meta,
        )
    )

    (
        C,
        covariate_names,
    ) = covariate_matrix(
        cohort,
        meta
        .subject_id
        .tolist(),
    )

    assets = {}

    for fold in range(
        1,
        CV.N_FOLDS + 1,
    ):

        (
            gene_model,
            gene_scaler,
        ) = CV.gene_checkpoint(
            cohort,
            fold,
        )

        mri_model = (
            CV.mri_checkpoint(
                cohort,
                fold,
            )
        )

        assets[
            str(
                fold
            )
        ] = {

            "gene_model":
                str(
                    gene_model
                    .relative_to(
                        ROOT
                    )
                ),

            "gene_scaler":
                str(
                    gene_scaler
                    .relative_to(
                        ROOT
                    )
                ),

            "mri_model":
                str(
                    mri_model
                    .relative_to(
                        ROOT
                    )
                ),
        }

    report = {

        "cohort":
            cohort,

        "development_n":
            int(
                len(
                    meta
                )
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

        "fold_counts": {
            str(k):
                int(v)

            for k, v
            in meta[
                "cv_fold"
            ]
            .value_counts()
            .sort_index()
            .items()
        },

        "gene_matrix_shape":
            list(
                gene_matrix.shape
            ),

        "covariate_matrix_shape":
            list(
                C.shape
            ),

        "covariates":
            covariate_names,

        "batch_covariate_included":
            False,

        "fold_agreement_verified":
            True,

        "branch_assets":
            assets,

        "original_cvae_script":
            str(
                CVAE_SCRIPT
                .relative_to(
                    ROOT
                )
            ),

        "original_cvae_script_sha256":
            CV.sha256(
                CVAE_SCRIPT
            ),

        "model_output_root":
            str(
                (
                    ROOT
                    / "Output"
                    / "Model"
                    / MODEL_GROUP
                    / RUN_TAG
                )
                .relative_to(
                    ROOT
                )
            ),

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
# CLI
# ============================================================

def main() -> None:

    ap = argparse.ArgumentParser(
        description=(
            "Age/sex/APOE4-adjusted IG-CVAE "
            "sensitivity analysis using the original "
            "branch checkpoints and original CVAE."
        )
    )

    ap.add_argument(
        "command",
        choices=[
            "VALIDATE",
            "DEVELOP",
            "TEST",
        ],
    )

    ap.add_argument(
        "--cohort",
        choices=[
            "ADNI",
            "AddNeuroMed",
            "ALL",
        ],
        default="ALL",
    )

    ap.add_argument(
        "--device",
        choices=[
            "auto",
            "cpu",
            "cuda",
        ],
        default="auto",
    )

    ap.add_argument(
        "--allow-locked-test",
        action="store_true",
    )

    ap.add_argument(
        "--bootstrap-replicates",
        type=int,
        default=2000,
    )

    args = ap.parse_args()

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
        device.type
        == "cuda"
        and not torch.cuda.is_available()
    ):

        raise RuntimeError(
            "CUDA requested but unavailable."
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
                args.allow_locked_test,
                args.bootstrap_replicates,
            )


if __name__ == "__main__":
    main()
