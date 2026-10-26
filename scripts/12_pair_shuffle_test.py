#!/usr/bin/env python3
"""
Participant-pairing permutation control for the ORIGINAL IG-CVAE.

Test whether the original MRI->GFV mapping contains participant-specific
cross-modal information beyond diagnosis-level structure.

For each development CV fold:
    1. Use the SAME frozen Gene and MRI branch representations as the
       original scripts/06_run_cvae.py.
    2. Preserve MRI order.
    3. Shuffle GFV training targets ONLY within diagnosis (AD within AD,
       CN within CN), with no participant allowed to remain paired with
       themselves.
    4. Retrain the ORIGINAL CVAE using the shuffled training pairs.
    5. Evaluate on the ORIGINAL, correctly paired validation fold.

Thus:
    - diagnosis distribution is preserved;
    - fixed folds are preserved;
    - branch checkpoints are preserved;
    - CVAE architecture/loss/hyperparameters are preserved;
    - participant MRI-GFV correspondence is destroyed only in training.

"""

from __future__ import annotations

import argparse
import importlib.util
import json
import logging
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import torch

from sklearn.exceptions import InconsistentVersionWarning
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler


# ============================================================
# Project
# ============================================================

ROOT = Path(__file__).resolve().parents[1]

CVAE_SCRIPT = (
    ROOT
    / "scripts"
    / "06_run_cvae.py"
)

MODEL_GROUP = "CVAE_pairshuffle"
RUN_TAG = "same_dx_training_target_v1"

BASE_SHUFFLE_SEED = 20261010
DEFAULT_N_SHUFFLES = 200


# Suppress only the known sklearn persistence-version warning
# from the already-frozen branch scalers.
warnings.filterwarnings(
    "ignore",
    category=InconsistentVersionWarning,
)


# ============================================================
# Load ORIGINAL CVAE
# ============================================================

def load_original_cvae():

    name = "igcvae_original_for_pairshuffle"

    spec = importlib.util.spec_from_file_location(
        name,
        CVAE_SCRIPT,
    )

    if (
        spec is None
        or spec.loader is None
    ):
        raise RuntimeError(
            f"Cannot import {CVAE_SCRIPT}"
        )

    mod = importlib.util.module_from_spec(
        spec
    )

    sys.modules[name] = mod

    spec.loader.exec_module(
        mod
    )

    return mod


CV = load_original_cvae()


# ============================================================
# Separate output tree
# ============================================================

def outdir(
    kind: str,
    cohort: str,
) -> Path:

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


# The original CVAE's MRI feature extractor uses CV.outdir()
# only to choose its preprocessing cache. Redirect that cache
# to this shuffle-control output tree so primary outputs are
# never modified.
_ORIGINAL_CV_OUTDIR = CV.outdir


def _pairshuffle_cv_outdir(
    kind: str,
    cohort: str,
) -> Path:

    if kind == "Intermediate":
        return outdir(
            "Intermediate",
            cohort,
        )

    return _ORIGINAL_CV_OUTDIR(
        kind,
        cohort,
    )


CV.outdir = _pairshuffle_cv_outdir


# ============================================================
# Utilities
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
) -> logging.Logger:

    log = logging.getLogger(
        f"pairshuffle_{cohort}"
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
        / "pair_shuffle.log",
        mode="a",
    )

    fh.setFormatter(
        fmt
    )

    log.addHandler(
        fh
    )

    return log


# ============================================================
# Same-diagnosis derangement
# ============================================================

def same_dx_permutation(
    y: np.ndarray,
    seed: int,
) -> np.ndarray:

    """
    Diagnosis-preserving derangement.

    Every training MRI is assigned the GFV of a DIFFERENT
    participant with the SAME diagnosis.
    """

    y = np.asarray(
        y,
        dtype=int,
    )

    rng = np.random.default_rng(
        seed
    )

    perm = np.arange(
        len(y)
    )

    for diagnosis in np.unique(
        y
    ):

        idx = np.flatnonzero(
            y == diagnosis
        )

        if len(idx) < 2:

            raise RuntimeError(
                f"Cannot derange diagnosis group "
                f"{diagnosis}: n={len(idx)}"
            )

        # Rejection sampling gives a true no-fixed-point
        # permutation within the diagnosis group.
        while True:

            shuffled = rng.permutation(
                idx
            )

            if np.all(
                shuffled != idx
            ):
                break

        perm[
            idx
        ] = shuffled

    # --------------------------------------------------------
    # Safety checks
    # --------------------------------------------------------

    original = np.arange(
        len(y)
    )

    if np.any(
        perm == original
    ):

        raise RuntimeError(
            "Derangement failure: "
            "at least one true pair remains."
        )

    if not np.array_equal(
        y[
            perm
        ],
        y,
    ):

        raise RuntimeError(
            "Diagnosis-preserving shuffle failure."
        )

    if not np.array_equal(
        np.sort(
            perm
        ),
        original,
    ):

        raise RuntimeError(
            "Permutation is not one-to-one."
        )

    return perm


# ============================================================
# Original primary CVAE metrics
# ============================================================

def primary_fold_metrics_path(
    cohort: str,
) -> Path:

    return (
        ROOT
        / "Output"
        / "Table"
        / "CVAE"
        / CV.RUN_TAG
        / cohort
        / "fold_metrics.csv"
    )


def primary_metrics(
    cohort: str,
) -> dict[str, float]:

    path = primary_fold_metrics_path(
        cohort
    )

    if not path.is_file():

        raise FileNotFoundError(
            path
        )

    df = pd.read_csv(
        path
    )

    required = [
        "mse",
        "mae",
        "cosine_similarity",
        "r2",
        "roc_auc",
        "pr_auc",
    ]

    missing = [
        name
        for name in required
        if name not in df.columns
    ]

    if missing:

        raise RuntimeError(
            f"{cohort}: primary CVAE fold metrics "
            f"missing {missing}"
        )

    if len(
        df
    ) != CV.N_FOLDS:

        raise RuntimeError(
            f"{cohort}: expected {CV.N_FOLDS} "
            f"primary CVAE folds, found {len(df)}"
        )

    return {

        metric:
            float(
                df[
                    metric
                ].mean()
            )

        for metric in required
    }


# ============================================================
# Prepare original fold representations once
# ============================================================

def prepare_folds(
    cohort: str,
    device: torch.device,
):

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

    folds = (
        meta[
            "cv_fold"
        ]
        .astype(int)
        .to_numpy()
    )

    y = (
        meta[
            "diagnosis_binary"
        ]
        .astype(int)
        .to_numpy()
    )

    prepared = {}
    provenance = {}

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
        # EXACT SAME original branch features as primary CVAE
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
        # EXACT SAME train-only scaling as primary CVAE
        # ----------------------------------------------------

        mfv_scaler = (
            StandardScaler()
        )

        gfv_scaler = (
            StandardScaler()
        )

        mfv_tr = (
            mfv_scaler
            .fit_transform(
                mfv_raw[
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
                gfv_raw[
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
                mfv_raw[
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
                gfv_raw[
                    va
                ]
            )
            .astype(
                np.float32
            )
        )

        prepared[
            fold
        ] = {

            "mfv_tr":
                mfv_tr,

            "gfv_tr":
                gfv_tr,

            "mfv_va":
                mfv_va,

            "gfv_va":
                gfv_va,

            "gfv_raw_va":
                gfv_raw[
                    va
                ].astype(
                    np.float32
                ),

            "gfv_scaler":
                gfv_scaler,

            "y_tr":
                y[
                    tr
                ],

            "y_va":
                y[
                    va
                ],
        }

        provenance[
            str(
                fold
            )
        ] = branch_info

    return (
        meta,
        prepared,
        provenance,
    )


# ============================================================
# One shuffle replicate
# ============================================================

def evaluate_shuffle(
    prepared: dict,
    shuffle_number: int,
    device: torch.device,
) -> dict:

    fold_rows = []

    for fold in range(
        1,
        CV.N_FOLDS + 1,
    ):

        d = prepared[
            fold
        ]

        # ----------------------------------------------------
        # Destroy MRI-GFV participant pairing ONLY within DX
        # ----------------------------------------------------

        permutation = (
            same_dx_permutation(

                d[
                    "y_tr"
                ],

                BASE_SHUFFLE_SEED
                + shuffle_number
                * 100
                + fold,
            )
        )

        shuffled_gfv_tr = (
            d[
                "gfv_tr"
            ][
                permutation
            ]
            .copy()
        )

        # ----------------------------------------------------
        # Train EXACT ORIGINAL CVAE.
        #
        # Same initialization seed as corresponding true fold.
        # Validation remains correctly paired.
        # ----------------------------------------------------

        (
            model,
            best_epoch,
        ) = CV.fit_fold(

            d[
                "mfv_tr"
            ],

            shuffled_gfv_tr,

            d[
                "mfv_va"
            ],

            d[
                "gfv_va"
            ],

            device,

            CV.CFG.seed
            + fold,
        )

        # ----------------------------------------------------
        # MRI-only generation on correctly paired validation
        # subjects
        # ----------------------------------------------------

        generated_scaled = (
            CV.predict_generated(
                model,
                d[
                    "mfv_va"
                ],
                device,
            )
        )

        generated_raw = (
            d[
                "gfv_scaler"
            ]
            .inverse_transform(
                generated_scaled
            )
            .astype(
                np.float32
            )
        )

        reconstruction = (
            CV.reconstruction_metrics(

                d[
                    "gfv_raw_va"
                ],

                generated_raw,
            )
        )

        # ----------------------------------------------------
        # Same diagnostic evaluator as primary CVAE.
        #
        # Classifier is trained on TRUE GFV representations.
        # Since shuffle is within diagnosis, diagnosis-level
        # information is intentionally preserved in the null.
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

            d[
                "gfv_tr"
            ],

            d[
                "y_tr"
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

                d[
                    "y_va"
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

                "best_epoch":
                    int(
                        best_epoch
                    ),

                "mse":
                    reconstruction[
                        "mse"
                    ],

                "mae":
                    reconstruction[
                        "mae"
                    ],

                "cosine_similarity":
                    reconstruction[
                        "cosine_similarity"
                    ],

                "r2":
                    reconstruction[
                        "r2"
                    ],

                "roc_auc":
                    classification[
                        "roc_auc"
                    ],

                "pr_auc":
                    classification[
                        "pr_auc"
                    ],
            }
        )

        del model

        if (
            device.type
            == "cuda"
        ):
            torch.cuda.empty_cache()

    fold_df = pd.DataFrame(
        fold_rows
    )

    return {

        "shuffle":
            int(
                shuffle_number
            ),

        "mse":
            float(
                fold_df[
                    "mse"
                ].mean()
            ),

        "mae":
            float(
                fold_df[
                    "mae"
                ].mean()
            ),

        "cosine_similarity":
            float(
                fold_df[
                    "cosine_similarity"
                ].mean()
            ),

        "r2":
            float(
                fold_df[
                    "r2"
                ].mean()
            ),

        "roc_auc":
            float(
                fold_df[
                    "roc_auc"
                ].mean()
            ),

        "pr_auc":
            float(
                fold_df[
                    "pr_auc"
                ].mean()
            ),
    }


# ============================================================
# Empirical permutation P value
# ============================================================

def empirical_p(
    true_value: float,
    null_values: np.ndarray,
    higher_is_better: bool,
) -> float:

    null_values = np.asarray(
        null_values,
        dtype=float,
    )

    if higher_is_better:

        # shuffled performs at least as well as true
        count = int(
            np.sum(
                null_values
                >= true_value
            )
        )

    else:

        # lower is better:
        # shuffled error <= true error
        count = int(
            np.sum(
                null_values
                <= true_value
            )
        )

    return float(
        (
            count
            + 1
        )
        /
        (
            len(
                null_values
            )
            + 1
        )
    )


# ============================================================
# Validation
# ============================================================

def validate(
    cohort: str,
) -> None:

    meta = CV.cvae_rows(
        cohort
    )

    CV.verify_fold_agreement(
        cohort,
        meta,
    )

    primary_path = (
        primary_fold_metrics_path(
            cohort
        )
    )

    true_metrics = (
        primary_metrics(
            cohort
        )
    )

    branch_assets = {}

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

        branch_assets[
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

        "test":
            "same_diagnosis_training_target_pair_shuffle",

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

        "fold_agreement_verified":
            True,

        "shuffle_rule":
            (
                "training GFV targets deranged "
                "within diagnosis; MRI order fixed"
            ),

        "validation_pairing":
            "true participant pairing",

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

        "primary_fold_metrics":
            str(
                primary_path
                .relative_to(
                    ROOT
                )
            ),

        "primary_fold_metrics_sha256":
            CV.sha256(
                primary_path
            ),

        "true_primary_metrics":
            true_metrics,

        "branch_assets":
            branch_assets,

        "output_root":
            str(
                (
                    ROOT
                    / "Output"
                    / "Table"
                    / MODEL_GROUP
                    / RUN_TAG
                    / cohort
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
# Run permutation control
# ============================================================

def run_cohort(
    cohort: str,
    device: torch.device,
    n_shuffles: int,
    restart: bool,
) -> None:

    log = setup_logger(
        cohort
    )

    table_dir = outdir(
        "Table",
        cohort,
    )

    qc_dir = outdir(
        "QC",
        cohort,
    )

    replicate_path = (
        table_dir
        / "shuffle_replicates.csv"
    )

    summary_path = (
        table_dir
        / "summary.json"
    )

    provenance_path = (
        qc_dir
        / "provenance.json"
    )

    # --------------------------------------------------------
    # Optional fresh restart
    # --------------------------------------------------------

    if restart:

        for path in (
            replicate_path,
            summary_path,
            provenance_path,
        ):

            if path.exists():
                path.unlink()

    # --------------------------------------------------------
    # Resume support
    # --------------------------------------------------------

    if replicate_path.exists():

        existing = pd.read_csv(
            replicate_path
        )

        completed = set(
            existing[
                "shuffle"
            ]
            .astype(int)
            .tolist()
        )

    else:

        completed = set()

    log.info(
        "%s pair-shuffle test starting.",
        cohort,
    )

    log.info(
        "Original CVAE: %s",
        CVAE_SCRIPT,
    )

    log.info(
        "Same-diagnosis derangement; "
        "true participant pairs destroyed in training."
    )

    log.info(
        "Locked test accessed: False"
    )

    log.info(
        "%d/%d shuffle replicates already complete.",
        len(
            completed
        ),
        n_shuffles,
    )

    # --------------------------------------------------------
    # Extract original branch representations ONCE
    # --------------------------------------------------------

    (
        meta,
        prepared,
        branch_provenance,
    ) = prepare_folds(
        cohort,
        device,
    )

    primary_path = (
        primary_fold_metrics_path(
            cohort
        )
    )

    true = primary_metrics(
        cohort
    )

    save_json(
        provenance_path,
        {

            "cohort":
                cohort,

            "test":
                "same_diagnosis_training_target_pair_shuffle",

            "run_tag":
                RUN_TAG,

            "n_requested_shuffles":
                int(
                    n_shuffles
                ),

            "base_shuffle_seed":
                int(
                    BASE_SHUFFLE_SEED
                ),

            "shuffle_rule":
                (
                    "Within each CV training partition, "
                    "GFV targets are assigned to a different "
                    "participant of the same diagnosis."
                ),

            "validation_rule":
                (
                    "Validation subjects retain their true "
                    "MRI-GFV correspondence."
                ),

            "model_rule":
                (
                    "Original 06_run_cvae.py fit_fold, "
                    "architecture, loss, hyperparameters, "
                    "folds, branch checkpoints, and fold model "
                    "seeds are reused unchanged."
                ),

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

            "primary_fold_metrics":
                str(
                    primary_path
                    .relative_to(
                        ROOT
                    )
                ),

            "primary_fold_metrics_sha256":
                CV.sha256(
                    primary_path
                ),

            "n_development":
                int(
                    len(
                        meta
                    )
                ),

            "branch_provenance":
                branch_provenance,

            "locked_test_accessed":
                False,
        },
    )

    # ========================================================
    # Null distribution
    # ========================================================

    for shuffle_number in range(
        1,
        n_shuffles + 1,
    ):

        if shuffle_number in completed:
            continue

        row = evaluate_shuffle(
            prepared,
            shuffle_number,
            device,
        )

        pd.DataFrame(
            [
                row
            ]
        ).to_csv(

            replicate_path,

            mode="a",

            header=(
                not replicate_path.exists()
            ),

            index=False,
        )

        completed.add(
            shuffle_number
        )

        if (
            shuffle_number == 1
            or shuffle_number % 25 == 0
            or shuffle_number == n_shuffles
        ):

            log.info(
                "shuffle %d/%d complete",
                shuffle_number,
                n_shuffles,
            )

    # ========================================================
    # Final null distribution
    # ========================================================

    null = pd.read_csv(
        replicate_path
    )

    null = (
        null[
            null[
                "shuffle"
            ]
            <= n_shuffles
        ]
        .sort_values(
            "shuffle"
        )
        .drop_duplicates(
            "shuffle",
            keep="last",
        )
        .reset_index(
            drop=True
        )
    )

    if (
        len(
            null
        )
        != n_shuffles
    ):

        raise RuntimeError(
            f"{cohort}: expected "
            f"{n_shuffles} shuffles, "
            f"found {len(null)}"
        )

    # --------------------------------------------------------
    # Metric directions
    # --------------------------------------------------------

    directions = {

        "mse":
            False,

        "mae":
            False,

        "cosine_similarity":
            True,

        "r2":
            True,

        "roc_auc":
            True,

        "pr_auc":
            True,
    }

    summary = {

        "cohort":
            cohort,

        "test":
            "same_diagnosis_training_target_pair_shuffle",

        "run_tag":
            RUN_TAG,

        "n_shuffles":
            int(
                n_shuffles
            ),

        "shuffle_rule":
            (
                "training-fold GFV targets "
                "deranged within diagnosis only"
            ),

        "validation_pairing":
            "true participant pairing",

        "locked_test_accessed":
            False,

        "metrics":
            {},
    }

    for (
        metric,
        higher_is_better,
    ) in directions.items():

        values = (
            null[
                metric
            ]
            .to_numpy(
                dtype=float
            )
        )

        true_value = float(
            true[
                metric
            ]
        )

        summary[
            "metrics"
        ][
            metric
        ] = {

            "direction":
                (
                    "higher_is_better"
                    if higher_is_better
                    else "lower_is_better"
                ),

            "true":
                true_value,

            "shuffle_mean":
                float(
                    values.mean()
                ),

            "shuffle_sd":
                float(
                    values.std(
                        ddof=1
                    )
                ),

            "shuffle_2p5":
                float(
                    np.quantile(
                        values,
                        0.025,
                    )
                ),

            "shuffle_97p5":
                float(
                    np.quantile(
                        values,
                        0.975,
                    )
                ),

            "true_minus_shuffle_mean":
                float(
                    true_value
                    - values.mean()
                ),

            "p_value_one_sided":
                empirical_p(
                    true_value,
                    values,
                    higher_is_better,
                ),
        }

    save_json(
        summary_path,
        summary,
    )

    m = summary[
        "metrics"
    ]

    log.info(
        "%s FINAL | "
        "MSE true=%.4f shuffle=%.4f p=%.4g | "
        "MAE true=%.4f shuffle=%.4f p=%.4g | "
        "Cos true=%.4f shuffle=%.4f p=%.4g | "
        "R2 true=%.4f shuffle=%.4f p=%.4g | "
        "AUC true=%.4f shuffle=%.4f p=%.4g | "
        "AP true=%.4f shuffle=%.4f p=%.4g",
        cohort,

        m[
            "mse"
        ][
            "true"
        ],
        m[
            "mse"
        ][
            "shuffle_mean"
        ],
        m[
            "mse"
        ][
            "p_value_one_sided"
        ],

        m[
            "mae"
        ][
            "true"
        ],
        m[
            "mae"
        ][
            "shuffle_mean"
        ],
        m[
            "mae"
        ][
            "p_value_one_sided"
        ],

        m[
            "cosine_similarity"
        ][
            "true"
        ],
        m[
            "cosine_similarity"
        ][
            "shuffle_mean"
        ],
        m[
            "cosine_similarity"
        ][
            "p_value_one_sided"
        ],

        m[
            "r2"
        ][
            "true"
        ],
        m[
            "r2"
        ][
            "shuffle_mean"
        ],
        m[
            "r2"
        ][
            "p_value_one_sided"
        ],

        m[
            "roc_auc"
        ][
            "true"
        ],
        m[
            "roc_auc"
        ][
            "shuffle_mean"
        ],
        m[
            "roc_auc"
        ][
            "p_value_one_sided"
        ],

        m[
            "pr_auc"
        ][
            "true"
        ],
        m[
            "pr_auc"
        ][
            "shuffle_mean"
        ],
        m[
            "pr_auc"
        ][
            "p_value_one_sided"
        ],
    )


# ============================================================
# CLI
# ============================================================

def main() -> None:

    ap = argparse.ArgumentParser(
        description=(
            "Same-diagnosis participant-pair "
            "shuffle control for original IG-CVAE."
        )
    )

    ap.add_argument(
        "command",
        choices=[
            "VALIDATE",
            "RUN",
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
        "--n-shuffles",
        type=int,
        default=DEFAULT_N_SHUFFLES,
    )

    ap.add_argument(
        "--restart",
        action="store_true",
        help=(
            "Delete existing partial shuffle "
            "results and start from shuffle 1."
        ),
    )

    args = ap.parse_args()

    if args.n_shuffles < 1:

        raise ValueError(
            "--n-shuffles must be >= 1"
        )

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

    if args.command == "VALIDATE":

        for cohort in cohorts:

            validate(
                cohort
            )

        return

    print(
        "PAIR SHUFFLE TEST — "
        "same-diagnosis training-target derangement"
    )

    print(
        f"Shuffles per cohort: {args.n_shuffles}"
    )

    print(
        "Locked test will NOT be accessed."
    )

    for cohort in cohorts:

        run_cohort(
            cohort,
            device,
            args.n_shuffles,
            args.restart,
        )


if __name__ == "__main__":
    main()
