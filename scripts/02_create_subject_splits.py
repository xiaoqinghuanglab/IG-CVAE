#!/usr/bin/env python3

"""
Create deterministic subject-level splits for the dual-cohort XIG-CVAE study.

Primary design
--------------
1. ADNI and AddNeuroMed are analyzed as independent parallel cohorts.
2. Each cohort receives model-specific development and locked-test sets.
3. Approximately 20% of each eligible model population is locked for testing.
4. Within each cohort, the CVAE test is a common paired core included in the
   corresponding gene- and MRI-branch locked tests.
5. Gene and MRI locked tests contain the paired core plus modality-specific
   subjects when necessary to reach approximately 20%.
6. Five-fold assignments are created independently within each cohort.
7. Paired development subjects retain the same fold across all three models.
8. No preprocessing fitting, DEG selection, model fitting, or hyperparameter
   selection may use locked-test subjects.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd


SEED = 20260716
SPLIT_VERSION = "dual_cohort_independent_v1"
TEST_FRACTION = 0.20
N_FOLDS = 5

EXPECTED_COUNTS = {
    ("ADNI", "cvae", "development"): {"AD": 60, "CN": 171},
    ("ADNI", "cvae", "locked_test"): {"AD": 25, "CN": 33},
    ("ADNI", "gene_branch", "development"): {"AD": 64, "CN": 178},
    ("ADNI", "gene_branch", "locked_test"): {"AD": 26, "CN": 34},
    ("ADNI", "mri_branch", "development"): {"AD": 502, "CN": 694},
    ("ADNI", "mri_branch", "locked_test"): {"AD": 129, "CN": 170},
    ("AddNeuroMed", "cvae", "development"): {"AD": 37, "CN": 41},
    ("AddNeuroMed", "cvae", "locked_test"): {"AD": 9, "CN": 11},
    ("AddNeuroMed", "gene_branch", "development"): {
        "AD": 145,
        "CN": 135,
    },
    ("AddNeuroMed", "gene_branch", "locked_test"): {
        "AD": 36,
        "CN": 34,
    },
    ("AddNeuroMed", "mri_branch", "development"): {
        "AD": 45,
        "CN": 44,
    },
    ("AddNeuroMed", "mri_branch", "locked_test"): {
        "AD": 11,
        "CN": 11,
    },
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create deterministic independent ADNI and AddNeuroMed splits."
    )
    parser.add_argument(
        "--root",
        type=Path,
        default=Path.cwd(),
        help="Project root. Default: current directory.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite an existing split manifest.",
    )
    return parser.parse_args()


def clean_text(series: pd.Series) -> pd.Series:
    return series.astype("string").str.strip()


def normalize_diagnosis(series: pd.Series) -> pd.Series:
    value = clean_text(series).str.upper()

    mapping = {
        "AD": "AD",
        "ALZHEIMER'S DISEASE": "AD",
        "ALZHEIMERS DISEASE": "AD",
        "CN": "CN",
        "CONTROL": "CN",
        "CTL": "CN",
        "COGNITIVELY NORMAL": "CN",
    }

    result = value.map(mapping)

    if result.isna().any():
        bad = sorted(value[result.isna()].dropna().unique().tolist())
        raise RuntimeError(
            "Unrecognized diagnosis values: "
            + ", ".join(map(str, bad))
        )

    return result


def normalize_sex(series: pd.Series) -> pd.Series:
    value = clean_text(series).str.upper()

    mapping = {
        "F": "Female",
        "FEMALE": "Female",
        "M": "Male",
        "MALE": "Male",
    }

    result = value.map(mapping).fillna("Missing")
    return result


def normalize_bids_id(value: object, subject_id: str) -> str:
    if pd.isna(value) or not str(value).strip():
        compact = subject_id.replace("_", "")
        return f"sub-{compact}"

    text = str(value).strip()
    return text if text.startswith("sub-") else f"sub-{text}"


def proportional_allocation(
    capacities: pd.Series,
    requested: int,
    rng: np.random.Generator,
) -> pd.Series:
    capacities = capacities.astype(int)

    if requested < 0 or requested > int(capacities.sum()):
        raise RuntimeError(
            f"Cannot select {requested} subjects from "
            f"{int(capacities.sum())} available subjects."
        )

    if requested == 0:
        return pd.Series(0, index=capacities.index, dtype=int)

    raw = capacities * requested / capacities.sum()
    allocation = np.floor(raw).astype(int)
    remaining = requested - int(allocation.sum())

    fractional = raw - allocation
    tie_break = pd.Series(
        rng.random(len(capacities)),
        index=capacities.index,
    )

    while remaining > 0:
        eligible = allocation < capacities

        if not eligible.any():
            raise RuntimeError("Unable to complete proportional allocation.")

        ordering = pd.DataFrame(
            {
                "fractional": fractional[eligible],
                "tie_break": tie_break[eligible],
            }
        ).sort_values(
            ["fractional", "tie_break"],
            ascending=[False, True],
        )

        chosen = ordering.index[0]
        allocation.loc[chosen] += 1
        fractional.loc[chosen] = -1
        remaining -= 1

    return allocation.astype(int)


def select_exact_test_subjects(
    frame: pd.DataFrame,
    diagnosis_targets: dict[str, int],
    rng: np.random.Generator,
) -> set[str]:
    """
    Select exact diagnosis counts while approximately preserving sex and age.

    Within each diagnosis, subjects are stratified by sex and age quintile.
    The requested diagnosis count is allocated proportionally across those
    demographic strata.
    """

    selected: set[str] = set()

    for diagnosis, target in diagnosis_targets.items():
        pool = frame.loc[
            frame["diagnosis"].eq(diagnosis)
        ].copy()

        if len(pool) < target:
            raise RuntimeError(
                f"Requested {target} {diagnosis} subjects, "
                f"but only {len(pool)} are available."
            )

        age_for_bins = pool["age"].copy()

        if age_for_bins.isna().any():
            age_for_bins = age_for_bins.fillna(age_for_bins.median())

        number_of_bins = min(5, max(1, len(pool)))

        age_rank = age_for_bins.rank(method="first")

        pool["_age_bin"] = pd.qcut(
            age_rank,
            q=number_of_bins,
            labels=False,
            duplicates="drop",
        ).astype(str)

        pool["_stratum"] = (
            pool["sex"].fillna("Missing").astype(str)
            + "|age"
            + pool["_age_bin"]
        )

        capacities = pool.groupby("_stratum").size()
        allocation = proportional_allocation(
            capacities,
            target,
            rng,
        )

        for stratum, number in allocation.items():
            if number == 0:
                continue

            candidates = pool.loc[
                pool["_stratum"].eq(stratum),
                "subject_id",
            ].to_numpy()

            chosen = rng.choice(
                candidates,
                size=int(number),
                replace=False,
            )

            selected.update(map(str, chosen))

    expected = int(sum(diagnosis_targets.values()))

    if len(selected) != expected:
        raise RuntimeError(
            f"Expected {expected} selected subjects; "
            f"obtained {len(selected)}."
        )

    return selected


def assign_development_folds(
    frame: pd.DataFrame,
    rng: np.random.Generator,
) -> dict[str, int]:
    """
    Assign five development folds using diagnosis/sex stratification.

    Fold balancing is performed without looking at gene expression or MRI
    values. The same mapping is reused for paired subjects in every model.
    """

    working = frame[
        ["subject_id", "diagnosis", "sex"]
    ].copy()

    working["_stratum"] = (
        working["diagnosis"].astype(str)
        + "|"
        + working["sex"].astype(str)
    )

    fold_map: dict[str, int] = {}
    total_counts = np.zeros(N_FOLDS, dtype=int)
    diagnosis_counts = {
        "AD": np.zeros(N_FOLDS, dtype=int),
        "CN": np.zeros(N_FOLDS, dtype=int),
    }

    for _, group in working.groupby("_stratum", sort=True):
        subject_ids = group["subject_id"].to_numpy().copy()
        rng.shuffle(subject_ids)

        diagnosis = str(group["diagnosis"].iloc[0])

        for subject_id in subject_ids:
            random_ties = rng.random(N_FOLDS)

            ordering = sorted(
                range(N_FOLDS),
                key=lambda fold: (
                    diagnosis_counts[diagnosis][fold],
                    total_counts[fold],
                    random_ties[fold],
                ),
            )

            chosen_fold = ordering[0]

            fold_map[str(subject_id)] = chosen_fold + 1
            total_counts[chosen_fold] += 1
            diagnosis_counts[diagnosis][chosen_fold] += 1

    if len(fold_map) != len(frame):
        raise RuntimeError("Development fold assignment is incomplete.")

    return fold_map


def prepare_adni_paired(root: Path) -> pd.DataFrame:
    path = root / "data/ADNI/Paired_data_metadata.csv"
    raw = pd.read_csv(path, low_memory=False)

    frame = pd.DataFrame(
        {
            "subject_id": clean_text(raw["subject_id"]),
            "diagnosis": normalize_diagnosis(raw["diagnosis"]),
            "sex": normalize_sex(raw["sex"]),
            "age": pd.to_numeric(raw["age"], errors="coerce"),
        }
    )

    frame["bids_id"] = [
        normalize_bids_id(value, subject)
        for value, subject in zip(
            raw["bids_id"],
            frame["subject_id"],
        )
    ]

    frame["mri_path"] = frame["bids_id"].map(
        lambda bids: str(
            (
                root
                / "data/ADNI/paired_MRI"
                / f"{bids}.nii.gz"
            ).resolve()
        )
    )

    frame["mri_available"] = frame["mri_path"].map(
        lambda path_value: Path(path_value).is_file()
    )

    frame["gene_available"] = True
    frame["paired_subject"] = frame["mri_available"]
    frame["source_subset"] = np.where(
        frame["mri_available"],
        "paired_MRI_GE",
        "gene_only",
    )

    return frame


def prepare_adni_mri_only(root: Path) -> pd.DataFrame:
    path = root / "data/ADNI/mri_only_metadata.csv"
    raw = pd.read_csv(path, low_memory=False)

    frame = pd.DataFrame(
        {
            "subject_id": clean_text(raw["subject_id"]),
            "diagnosis": normalize_diagnosis(raw["diagnosis"]),
            "sex": normalize_sex(raw["sex"]),
            "age": pd.to_numeric(raw["age"], errors="coerce"),
        }
    )

    frame["bids_id"] = [
        normalize_bids_id(value, subject)
        for value, subject in zip(
            raw["bids_id"],
            frame["subject_id"],
        )
    ]

    frame["mri_path"] = frame["bids_id"].map(
        lambda bids: str(
            (
                root
                / "data/ADNI/only_MRI"
                / f"{bids}.nii.gz"
            ).resolve()
        )
    )

    frame["mri_available"] = frame["mri_path"].map(
        lambda path_value: Path(path_value).is_file()
    )

    frame["gene_available"] = False
    frame["paired_subject"] = False
    frame["source_subset"] = "MRI_only"

    return frame


def prepare_anm_paired(root: Path) -> pd.DataFrame:
    path = root / "data/AdNeuroMed/Paired_data_metadata.csv"
    raw = pd.read_csv(path, low_memory=False)

    frame = pd.DataFrame(
        {
            "subject_id": clean_text(raw["Subject_ID"]),
            "diagnosis": normalize_diagnosis(raw["Diagnosis"]),
            "sex": normalize_sex(raw["Sex"]),
            "age": pd.to_numeric(raw["Age"], errors="coerce"),
        }
    )

    frame["bids_id"] = [
        normalize_bids_id(value, subject)
        for value, subject in zip(
            raw["participant_id"],
            frame["subject_id"],
        )
    ]

    frame["mri_path"] = frame["bids_id"].map(
        lambda bids: str(
            (
                root
                / "data/AdNeuroMed/paired_MRI"
                / f"{bids}.nii.gz"
            ).resolve()
        )
    )

    frame["mri_available"] = frame["mri_path"].map(
        lambda path_value: Path(path_value).is_file()
    )

    frame["gene_available"] = True
    frame["paired_subject"] = True
    frame["source_subset"] = "paired_MRI_GE"

    return frame


def prepare_anm_gene(root: Path, paired: pd.DataFrame) -> pd.DataFrame:
    path = root / "data/AdNeuroMed/gene_expression_metadata.csv"
    raw = pd.read_csv(path, low_memory=False)

    frame = pd.DataFrame(
        {
            "subject_id": clean_text(raw["Subject_ID"]),
            "diagnosis": normalize_diagnosis(raw["Diagnosis"]),
            "sex": normalize_sex(raw["Sex"]),
            "age": pd.to_numeric(raw["Age"], errors="coerce"),
        }
    )

    paired_path_map = paired.set_index("subject_id")["mri_path"]
    paired_ids = set(paired["subject_id"])

    frame["bids_id"] = frame["subject_id"].map(
        lambda subject: f"sub-{subject}"
    )

    frame["mri_path"] = (
        frame["subject_id"].map(paired_path_map).fillna("")
    )

    frame["mri_available"] = frame["subject_id"].isin(paired_ids)
    frame["gene_available"] = True
    frame["paired_subject"] = frame["subject_id"].isin(paired_ids)
    frame["source_subset"] = np.where(
        frame["paired_subject"],
        "paired_MRI_GE",
        "gene_only",
    )

    return frame


def prepare_anm_mri_only(root: Path) -> pd.DataFrame:
    path = root / "data/AdNeuroMed/mri_only_metadata.csv"
    raw = pd.read_csv(path, low_memory=False)

    frame = pd.DataFrame(
        {
            "subject_id": clean_text(raw["Subject_ID"]),
            "diagnosis": normalize_diagnosis(raw["Diagnosis"]),
            "sex": normalize_sex(raw["Sex"]),
            "age": pd.to_numeric(raw["Age"], errors="coerce"),
        }
    )

    frame["bids_id"] = [
        normalize_bids_id(value, subject)
        for value, subject in zip(
            raw["participant_id"],
            frame["subject_id"],
        )
    ]

    frame["mri_path"] = frame["bids_id"].map(
        lambda bids: str(
            (
                root
                / "data/AdNeuroMed/only_MRI"
                / f"{bids}.nii.gz"
            ).resolve()
        )
    )

    frame["mri_available"] = frame["mri_path"].map(
        lambda path_value: Path(path_value).is_file()
    )

    frame["gene_available"] = False
    frame["paired_subject"] = False
    frame["source_subset"] = "MRI_only"

    return frame


def manifest_rows(
    frame: pd.DataFrame,
    cohort: str,
    model: str,
    role: str,
    component: str,
    fold_map: dict[str, int] | None = None,
) -> pd.DataFrame:
    result = frame.copy()

    result.insert(0, "split_version", SPLIT_VERSION)
    result.insert(1, "seed", SEED)
    result.insert(2, "cohort", cohort)
    result.insert(3, "model", model)
    result.insert(4, "role", role)
    result.insert(5, "split_component", component)

    if fold_map is None:
        result["cv_fold"] = pd.Series(
            pd.NA,
            index=result.index,
            dtype="Int64",
        )
    else:
        result["cv_fold"] = (
            result["subject_id"]
            .map(fold_map)
            .astype("Int64")
        )

        if result["cv_fold"].isna().any():
            missing = result.loc[
                result["cv_fold"].isna(),
                "subject_id",
            ].tolist()

            raise RuntimeError(
                "Missing development fold assignments: "
                + ", ".join(missing[:10])
            )

    result["diagnosis_binary"] = (
        result["diagnosis"].eq("AD").astype(int)
    )

    columns = [
        "split_version",
        "seed",
        "cohort",
        "model",
        "role",
        "split_component",
        "cv_fold",
        "subject_id",
        "bids_id",
        "diagnosis",
        "diagnosis_binary",
        "sex",
        "age",
        "source_subset",
        "paired_subject",
        "gene_available",
        "mri_available",
        "mri_path",
    ]

    return result[columns]


def count_diagnosis(frame: pd.DataFrame) -> dict[str, int]:
    counts = frame["diagnosis"].value_counts().to_dict()
    return {
        "AD": int(counts.get("AD", 0)),
        "CN": int(counts.get("CN", 0)),
    }


def standardized_mean_difference(
    development: pd.Series,
    test: pd.Series,
) -> float:
    development = pd.to_numeric(
        development,
        errors="coerce",
    ).dropna()

    test = pd.to_numeric(
        test,
        errors="coerce",
    ).dropna()

    if len(development) < 2 or len(test) < 2:
        return float("nan")

    pooled_sd = np.sqrt(
        (
            development.var(ddof=1)
            + test.var(ddof=1)
        )
        / 2
    )

    if pooled_sd == 0:
        return 0.0

    return float(
        (development.mean() - test.mean())
        / pooled_sd
    )


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()

    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)

    return digest.hexdigest()


def validate_source_data(
    adni_paired: pd.DataFrame,
    adni_mri_only: pd.DataFrame,
    anm_paired: pd.DataFrame,
    anm_gene: pd.DataFrame,
    anm_mri_only: pd.DataFrame,
) -> None:
    checks = {
        "ADNI gene subjects": (len(adni_paired), 302),
        "ADNI paired MRI-GE subjects": (
            int(adni_paired["mri_available"].sum()),
            289,
        ),
        "ADNI gene-only subjects": (
            int((~adni_paired["mri_available"]).sum()),
            13,
        ),
        "ADNI MRI-only subjects": (len(adni_mri_only), 1206),
        "AddNeuroMed paired subjects": (len(anm_paired), 98),
        "AddNeuroMed gene subjects": (len(anm_gene), 350),
        "AddNeuroMed MRI-only subjects": (len(anm_mri_only), 13),
    }

    failures = []

    for label, (observed, expected) in checks.items():
        if observed != expected:
            failures.append(
                f"{label}: expected {expected}, observed {observed}"
            )

    if failures:
        raise RuntimeError(
            "Source-data counts changed:\n  "
            + "\n  ".join(failures)
        )

    for label, frame in {
        "ADNI MRI-only": adni_mri_only,
        "AddNeuroMed paired": anm_paired,
        "AddNeuroMed MRI-only": anm_mri_only,
    }.items():
        missing = frame.loc[
            ~frame["mri_available"],
            ["subject_id", "mri_path"],
        ]

        if not missing.empty:
            first = missing.iloc[0]
            raise RuntimeError(
                f"{label} MRI missing for "
                f"{first['subject_id']}: {first['mri_path']}"
            )


def validate_manifest(
    manifest: pd.DataFrame,
    paired_test_ids: set[str],
    paired_fold_map: dict[str, int],
    anm_paired_test_ids: set[str],
    anm_paired_fold_map: dict[str, int],
) -> list[str]:
    validation_lines: list[str] = []

    duplicates = manifest.duplicated(
        ["cohort", "model", "subject_id"]
    ).sum()

    if duplicates:
        raise RuntimeError(
            f"Duplicate model-subject rows: {duplicates}"
        )

    validation_lines.append(
        f"Duplicate model-subject rows: {duplicates}"
    )

    for key, expected in EXPECTED_COUNTS.items():
        cohort, model, role = key

        subset = manifest.loc[
            manifest["cohort"].eq(cohort)
            & manifest["model"].eq(model)
            & manifest["role"].eq(role)
        ]

        observed = count_diagnosis(subset)

        if observed != expected:
            raise RuntimeError(
                f"Count mismatch for {cohort}/{model}/{role}: "
                f"expected {expected}, observed {observed}"
            )

    validation_lines.append(
        "All prespecified model/role/diagnosis counts: PASS"
    )

    if set(manifest["role"]) != {
        "development",
        "locked_test",
    }:
        raise RuntimeError(
            "Manifest contains a role other than development "
            "or locked_test."
        )

    validation_lines.append(
        "Only development and locked_test roles are present: PASS"
    )

    cohort_cores = {
        "ADNI": {
            "test_ids": paired_test_ids,
            "fold_map": paired_fold_map,
            "expected_test_size": 58,
        },
        "AddNeuroMed": {
            "test_ids": anm_paired_test_ids,
            "fold_map": anm_paired_fold_map,
            "expected_test_size": 20,
        },
    }

    for cohort, specification in cohort_cores.items():
        cohort_manifest = manifest.loc[
            manifest["cohort"].eq(cohort)
        ]

        for model in [
            "cvae",
            "gene_branch",
            "mri_branch",
        ]:
            development_ids = set(
                cohort_manifest.loc[
                    cohort_manifest["model"].eq(model)
                    & cohort_manifest["role"].eq("development"),
                    "subject_id",
                ]
            )

            test_ids = set(
                cohort_manifest.loc[
                    cohort_manifest["model"].eq(model)
                    & cohort_manifest["role"].eq("locked_test"),
                    "subject_id",
                ]
            )

            overlap = development_ids & test_ids

            if overlap:
                raise RuntimeError(
                    f"{cohort}/{model} development-test "
                    f"overlap: {len(overlap)}"
                )

            if not specification["test_ids"].issubset(
                test_ids
            ):
                raise RuntimeError(
                    f"{cohort} paired test core is incomplete "
                    f"in {model}."
                )

        validation_lines.append(
            f"{cohort} within-model development/test overlap: 0"
        )

        validation_lines.append(
            f"{cohort} common "
            f"{specification['expected_test_size']}-subject "
            "paired test core across all models: PASS"
        )

        paired_development_ids = set(
            specification["fold_map"]
        )

        for model in [
            "cvae",
            "gene_branch",
            "mri_branch",
        ]:
            subset = cohort_manifest.loc[
                cohort_manifest["model"].eq(model)
                & cohort_manifest["role"].eq("development")
                & cohort_manifest["subject_id"].isin(
                    paired_development_ids
                ),
                ["subject_id", "cv_fold"],
            ]

            observed = dict(
                zip(
                    subset["subject_id"],
                    subset["cv_fold"].astype(int),
                )
            )

            if observed != specification["fold_map"]:
                raise RuntimeError(
                    f"{cohort} paired fold assignments "
                    f"differ in {model}."
                )

        validation_lines.append(
            f"{cohort} paired development folds are "
            "aligned across all models: PASS"
        )

    development = manifest.loc[
        manifest["role"].eq("development")
    ]

    if development["cv_fold"].isna().any():
        raise RuntimeError(
            "At least one development subject lacks a fold."
        )

    if not development["cv_fold"].between(
        1,
        N_FOLDS,
    ).all():
        raise RuntimeError(
            "Invalid development fold number."
        )

    validation_lines.append(
        "All development folds are between 1 and 5: PASS"
    )

    locked_test = manifest.loc[
        manifest["role"].eq("locked_test")
    ]

    if locked_test["cv_fold"].notna().any():
        raise RuntimeError(
            "A locked-test subject received a CV fold."
        )

    validation_lines.append(
        "Locked-test subjects have no CV fold: PASS"
    )

    mri_rows = manifest.loc[
        manifest["mri_available"]
    ]

    missing_paths = [
        path_value
        for path_value in mri_rows["mri_path"]
        if not Path(path_value).is_file()
    ]

    if missing_paths:
        raise RuntimeError(
            f"{len(missing_paths)} manifest MRI paths "
            "do not exist."
        )

    validation_lines.append(
        "MRI path existence: PASS"
    )

    return validation_lines


def main() -> int:
    args = parse_args()
    root = args.root.resolve()

    output_dir = root / "Output/Intermediate/Splits"
    output_dir.mkdir(parents=True, exist_ok=True)

    manifest_path = output_dir / "model_split_manifest.csv"
    summary_path = output_dir / "split_summary.txt"
    config_path = output_dir / "split_config.json"
    hashes_path = output_dir / "split_hashes.sha256"
    log_path = output_dir / "split_run.log"

    managed_outputs = [
        manifest_path,
        summary_path,
        config_path,
        hashes_path,
        log_path,
    ]

    existing = [path for path in managed_outputs if path.exists()]

    if existing and not args.force:
        names = "\n  ".join(str(path) for path in existing)
        raise RuntimeError(
            "Split outputs already exist. Use --force to replace them:\n  "
            + names
        )

    rng = np.random.default_rng(SEED)

    adni_paired_all = prepare_adni_paired(root)
    adni_paired = adni_paired_all.loc[
        adni_paired_all["mri_available"]
    ].copy()

    adni_gene_only = adni_paired_all.loc[
        ~adni_paired_all["mri_available"]
    ].copy()

    adni_mri_only = prepare_adni_mri_only(root)

    anm_paired = prepare_anm_paired(root)
    anm_gene = prepare_anm_gene(root, anm_paired)
    anm_mri_only = prepare_anm_mri_only(root)

    validate_source_data(
        adni_paired_all,
        adni_mri_only,
        anm_paired,
        anm_gene,
        anm_mri_only,
    )

    # Common paired locked-test core: 58/289 subjects.
    paired_test_ids = select_exact_test_subjects(
        adni_paired,
        {"AD": 25, "CN": 33},
        rng,
    )

    paired_test = adni_paired.loc[
        adni_paired["subject_id"].isin(paired_test_ids)
    ].copy()

    paired_development = adni_paired.loc[
        ~adni_paired["subject_id"].isin(paired_test_ids)
    ].copy()

    # Gene-only extension: 2 subjects, yielding 60/302 gene test subjects.
    gene_extra_test_ids = select_exact_test_subjects(
        adni_gene_only,
        {"AD": 1, "CN": 1},
        rng,
    )

    gene_extra_test = adni_gene_only.loc[
        adni_gene_only["subject_id"].isin(gene_extra_test_ids)
    ].copy()

    gene_extra_development = adni_gene_only.loc[
        ~adni_gene_only["subject_id"].isin(gene_extra_test_ids)
    ].copy()

    # MRI-only extension: 241 subjects, yielding 299/1495 MRI test subjects.
    mri_extra_test_ids = select_exact_test_subjects(
        adni_mri_only,
        {"AD": 104, "CN": 137},
        rng,
    )

    mri_extra_test = adni_mri_only.loc[
        adni_mri_only["subject_id"].isin(mri_extra_test_ids)
    ].copy()

    mri_extra_development = adni_mri_only.loc[
        ~adni_mri_only["subject_id"].isin(mri_extra_test_ids)
    ].copy()

    paired_fold_map = assign_development_folds(
        paired_development,
        rng,
    )

    gene_extra_fold_map = assign_development_folds(
        gene_extra_development,
        rng,
    )

    mri_extra_fold_map = assign_development_folds(
        mri_extra_development,
        rng,
    )

    gene_fold_map = {
        **paired_fold_map,
        **gene_extra_fold_map,
    }

    mri_fold_map = {
        **paired_fold_map,
        **mri_extra_fold_map,
    }

    rows = []

    # ADNI CVAE.
    rows.append(
        manifest_rows(
            paired_development,
            "ADNI",
            "cvae",
            "development",
            "common_paired_development",
            paired_fold_map,
        )
    )

    rows.append(
        manifest_rows(
            paired_test,
            "ADNI",
            "cvae",
            "locked_test",
            "common_paired_test_core",
        )
    )

    # ADNI gene branch.
    rows.append(
        manifest_rows(
            paired_development,
            "ADNI",
            "gene_branch",
            "development",
            "common_paired_development",
            gene_fold_map,
        )
    )

    rows.append(
        manifest_rows(
            gene_extra_development,
            "ADNI",
            "gene_branch",
            "development",
            "gene_only_development_extension",
            gene_fold_map,
        )
    )

    rows.append(
        manifest_rows(
            paired_test,
            "ADNI",
            "gene_branch",
            "locked_test",
            "common_paired_test_core",
        )
    )

    rows.append(
        manifest_rows(
            gene_extra_test,
            "ADNI",
            "gene_branch",
            "locked_test",
            "gene_only_test_extension",
        )
    )

    # ADNI MRI branch.
    rows.append(
        manifest_rows(
            paired_development,
            "ADNI",
            "mri_branch",
            "development",
            "common_paired_development",
            mri_fold_map,
        )
    )

    rows.append(
        manifest_rows(
            mri_extra_development,
            "ADNI",
            "mri_branch",
            "development",
            "MRI_only_development_extension",
            mri_fold_map,
        )
    )

    rows.append(
        manifest_rows(
            paired_test,
            "ADNI",
            "mri_branch",
            "locked_test",
            "common_paired_test_core",
        )
    )

    rows.append(
        manifest_rows(
            mri_extra_test,
            "ADNI",
            "mri_branch",
            "locked_test",
            "MRI_only_test_extension",
        )
    )

    # AddNeuroMed common paired locked-test core: 20/98 subjects.
    anm_paired_test_ids = select_exact_test_subjects(
        anm_paired,
        {"AD": 9, "CN": 11},
        rng,
    )

    anm_paired_test = anm_paired.loc[
        anm_paired["subject_id"].isin(
            anm_paired_test_ids
        )
    ].copy()

    anm_paired_development = anm_paired.loc[
        ~anm_paired["subject_id"].isin(
            anm_paired_test_ids
        )
    ].copy()

    # AddNeuroMed gene-only population excludes the 98 paired subjects.
    anm_gene_only = anm_gene.loc[
        ~anm_gene["paired_subject"]
    ].copy()

    if len(anm_gene_only) != 252:
        raise RuntimeError(
            "Expected 252 AddNeuroMed gene-only subjects; "
            f"observed {len(anm_gene_only)}."
        )

    # Gene-only extension: 50 subjects, producing 70/350 gene tests.
    anm_gene_extra_test_ids = select_exact_test_subjects(
        anm_gene_only,
        {"AD": 27, "CN": 23},
        rng,
    )

    anm_gene_extra_test = anm_gene_only.loc[
        anm_gene_only["subject_id"].isin(
            anm_gene_extra_test_ids
        )
    ].copy()

    anm_gene_extra_development = anm_gene_only.loc[
        ~anm_gene_only["subject_id"].isin(
            anm_gene_extra_test_ids
        )
    ].copy()

    # MRI-only extension: two AD subjects, producing 22/111 MRI tests.
    anm_mri_extra_test_ids = select_exact_test_subjects(
        anm_mri_only,
        {"AD": 2, "CN": 0},
        rng,
    )

    anm_mri_extra_test = anm_mri_only.loc[
        anm_mri_only["subject_id"].isin(
            anm_mri_extra_test_ids
        )
    ].copy()

    anm_mri_extra_development = anm_mri_only.loc[
        ~anm_mri_only["subject_id"].isin(
            anm_mri_extra_test_ids
        )
    ].copy()

    anm_paired_fold_map = assign_development_folds(
        anm_paired_development,
        rng,
    )

    anm_gene_extra_fold_map = assign_development_folds(
        anm_gene_extra_development,
        rng,
    )

    anm_mri_extra_fold_map = assign_development_folds(
        anm_mri_extra_development,
        rng,
    )

    anm_gene_fold_map = {
        **anm_paired_fold_map,
        **anm_gene_extra_fold_map,
    }

    anm_mri_fold_map = {
        **anm_paired_fold_map,
        **anm_mri_extra_fold_map,
    }

    # AddNeuroMed CVAE.
    rows.append(
        manifest_rows(
            anm_paired_development,
            "AddNeuroMed",
            "cvae",
            "development",
            "common_paired_development",
            anm_paired_fold_map,
        )
    )

    rows.append(
        manifest_rows(
            anm_paired_test,
            "AddNeuroMed",
            "cvae",
            "locked_test",
            "common_paired_test_core",
        )
    )

    # AddNeuroMed gene branch.
    rows.append(
        manifest_rows(
            anm_paired_development,
            "AddNeuroMed",
            "gene_branch",
            "development",
            "common_paired_development",
            anm_gene_fold_map,
        )
    )

    rows.append(
        manifest_rows(
            anm_gene_extra_development,
            "AddNeuroMed",
            "gene_branch",
            "development",
            "gene_only_development_extension",
            anm_gene_fold_map,
        )
    )

    rows.append(
        manifest_rows(
            anm_paired_test,
            "AddNeuroMed",
            "gene_branch",
            "locked_test",
            "common_paired_test_core",
        )
    )

    rows.append(
        manifest_rows(
            anm_gene_extra_test,
            "AddNeuroMed",
            "gene_branch",
            "locked_test",
            "gene_only_test_extension",
        )
    )

    # AddNeuroMed MRI branch.
    rows.append(
        manifest_rows(
            anm_paired_development,
            "AddNeuroMed",
            "mri_branch",
            "development",
            "common_paired_development",
            anm_mri_fold_map,
        )
    )

    rows.append(
        manifest_rows(
            anm_mri_extra_development,
            "AddNeuroMed",
            "mri_branch",
            "development",
            "MRI_only_development_extension",
            anm_mri_fold_map,
        )
    )

    rows.append(
        manifest_rows(
            anm_paired_test,
            "AddNeuroMed",
            "mri_branch",
            "locked_test",
            "common_paired_test_core",
        )
    )

    rows.append(
        manifest_rows(
            anm_mri_extra_test,
            "AddNeuroMed",
            "mri_branch",
            "locked_test",
            "MRI_only_test_extension",
        )
    )

    manifest = pd.concat(rows, ignore_index=True)

    cohort_order = pd.Categorical(
        manifest["cohort"],
        categories=["ADNI", "AddNeuroMed"],
        ordered=True,
    )

    model_order = pd.Categorical(
        manifest["model"],
        categories=["cvae", "gene_branch", "mri_branch"],
        ordered=True,
    )

    role_order = pd.Categorical(
        manifest["role"],
        categories=["development", "locked_test"],
        ordered=True,
    )

    manifest = (
        manifest.assign(
            _cohort_order=cohort_order,
            _model_order=model_order,
            _role_order=role_order,
        )
        .sort_values(
            [
                "_cohort_order",
                "_model_order",
                "_role_order",
                "diagnosis",
                "subject_id",
            ]
        )
        .drop(
            columns=[
                "_cohort_order",
                "_model_order",
                "_role_order",
            ]
        )
        .reset_index(drop=True)
    )

    validation_lines = validate_manifest(
        manifest,
        paired_test_ids,
        paired_fold_map,
        anm_paired_test_ids,
        anm_paired_fold_map,
    )

    config = {
        "split_version": SPLIT_VERSION,
        "seed": SEED,
        "locked_test_fraction_target": TEST_FRACTION,
        "development_folds": N_FOLDS,
        "cohort_policy": (
            "ADNI and AddNeuroMed are independent parallel "
            "development and locked-test cohorts."
        ),
        "test_design": (
            "Model-specific approximately 20% locked tests "
            "containing a common paired core within each cohort."
        ),
        "cohort_test_designs": {
            "ADNI": {
                "common_paired_test_core": {
                    "total": 58,
                    "AD": 25,
                    "CN": 33,
                },
                "models": {
                    "cvae": {
                        "population": 289,
                        "test_total": 58,
                        "AD": 25,
                        "CN": 33,
                    },
                    "gene_branch": {
                        "population": 302,
                        "test_total": 60,
                        "AD": 26,
                        "CN": 34,
                    },
                    "mri_branch": {
                        "population": 1495,
                        "test_total": 299,
                        "AD": 129,
                        "CN": 170,
                    },
                },
            },
            "AddNeuroMed": {
                "common_paired_test_core": {
                    "total": 20,
                    "AD": 9,
                    "CN": 11,
                },
                "models": {
                    "cvae": {
                        "population": 98,
                        "test_total": 20,
                        "AD": 9,
                        "CN": 11,
                    },
                    "gene_branch": {
                        "population": 350,
                        "test_total": 70,
                        "AD": 36,
                        "CN": 34,
                    },
                    "mri_branch": {
                        "population": 111,
                        "test_total": 22,
                        "AD": 11,
                        "CN": 11,
                    },
                },
            },
        },
        "selection_information_used": [
            "diagnosis",
            "sex",
            "age",
        ],
        "selection_information_not_used": [
            "gene-expression values",
            "MRI voxel values",
            "model predictions",
            "model performance",
        ],
        "feature_selection_scope": (
            "Development subjects within each cohort independently"
        ),
        "preprocessing_fit_scope": (
            "Each cohort and training fold independently; "
            "locked-test transformations use parameters learned "
            "from the corresponding development cohort only."
        ),
        "status": "locked_before_feature_selection_and_modeling",
    }

    summary_lines = [
        "DUAL-COHORT LOCKED SUBJECT SPLIT SUMMARY",
        "=" * 72,
        f"Split version: {SPLIT_VERSION}",
        f"Seed: {SEED}",
        f"Target locked-test fraction: {TEST_FRACTION:.0%}",
        f"Development folds: {N_FOLDS}",
        "",
        "DESIGN",
        (
            "  ADNI and AddNeuroMed are independent parallel "
            "development/test cohorts."
        ),
        (
            "  ADNI paired test core: "
            "58 subjects (25 AD, 33 CN)"
        ),
        (
            "  AddNeuroMed paired test core: "
            "20 subjects (9 AD, 11 CN)"
        ),
        (
            "  Paired development folds are aligned across gene, "
            "MRI, and CVAE models within each cohort."
        ),
        "",
        "MODEL POPULATIONS",
    ]

    reporting_order = []

    for cohort in ["ADNI", "AddNeuroMed"]:
        for model in [
            "cvae",
            "gene_branch",
            "mri_branch",
        ]:
            for role in [
                "development",
                "locked_test",
            ]:
                reporting_order.append(
                    (cohort, model, role)
                )

    summary_lines.append(
        f"{'Cohort':<14} {'Model':<13} {'Role':<13} "
        f"{'Total':>6} {'AD':>5} {'CN':>5}"
    )

    for cohort, model, role in reporting_order:
        subset = manifest.loc[
            manifest["cohort"].eq(cohort)
            & manifest["model"].eq(model)
            & manifest["role"].eq(role)
        ]

        counts = count_diagnosis(subset)

        summary_lines.append(
            f"{cohort:<14} {model:<13} {role:<13} "
            f"{len(subset):>6} "
            f"{counts['AD']:>5} {counts['CN']:>5}"
        )

    summary_lines.extend(
        [
            "",
            "MODEL-SPECIFIC TEST QC",
        ]
    )

    for cohort in ["ADNI", "AddNeuroMed"]:
        summary_lines.append(f"{cohort}:")

        for model in [
            "cvae",
            "gene_branch",
            "mri_branch",
        ]:
            development = manifest.loc[
                manifest["cohort"].eq(cohort)
                & manifest["model"].eq(model)
                & manifest["role"].eq("development")
            ]

            locked_test = manifest.loc[
                manifest["cohort"].eq(cohort)
                & manifest["model"].eq(model)
                & manifest["role"].eq("locked_test")
            ]

            full_size = len(development) + len(locked_test)
            fraction = len(locked_test) / full_size

            age_smd = standardized_mean_difference(
                development["age"],
                locked_test["age"],
            )

            development_sex = (
                development["sex"].value_counts().to_dict()
            )

            test_sex = (
                locked_test["sex"].value_counts().to_dict()
            )

            summary_lines.extend(
                [
                    f"  {model}:",
                    (
                        f"    Test fraction: "
                        f"{len(locked_test)}/{full_size} "
                        f"({fraction:.2%})"
                    ),
                    (
                        "    Development diagnosis: "
                        + str(count_diagnosis(development))
                    ),
                    (
                        "    Locked-test diagnosis: "
                        + str(count_diagnosis(locked_test))
                    ),
                    (
                        f"    Development sex: "
                        f"{development_sex}"
                    ),
                    f"    Locked-test sex: {test_sex}",
                    (
                        "    Age standardized difference "
                        f"(development - test): {age_smd:.4f}"
                    ),
                ]
            )

    summary_lines.extend(
        [
            "",
            "FIVE-FOLD DEVELOPMENT COUNTS",
        ]
    )

    for cohort in ["ADNI", "AddNeuroMed"]:
        summary_lines.append(f"{cohort}:")

        for model in [
            "cvae",
            "gene_branch",
            "mri_branch",
        ]:
            development = manifest.loc[
                manifest["cohort"].eq(cohort)
                & manifest["model"].eq(model)
                & manifest["role"].eq("development")
            ]

            fold_table = (
                development.groupby(
                    ["cv_fold", "diagnosis"],
                    observed=True,
                )
                .size()
                .unstack(fill_value=0)
            )

            summary_lines.append(f"  {model}:")

            for fold in range(1, N_FOLDS + 1):
                ad_count = int(
                    fold_table.loc[fold, "AD"]
                    if fold in fold_table.index
                    and "AD" in fold_table.columns
                    else 0
                )

                cn_count = int(
                    fold_table.loc[fold, "CN"]
                    if fold in fold_table.index
                    and "CN" in fold_table.columns
                    else 0
                )

                summary_lines.append(
                    f"    Fold {fold}: "
                    f"AD={ad_count}, CN={cn_count}, "
                    f"Total={ad_count + cn_count}"
                )

    summary_lines.extend(
        [
            "",
            "VALIDATION CHECKS",
            *[
                f"  {line}"
                for line in validation_lines
            ],
            "",
            "LOCKING RULE",
            (
                "  Assignments must not be changed after "
                "inspecting gene values, MRI values, "
                "predictions, or performance."
            ),
            (
                "  DEG selection and all fitted preprocessing "
                "must use development subjects from the "
                "corresponding cohort only."
            ),
            "=" * 72,
        ]
    )

    summary_text = "\n".join(summary_lines) + "\n"

    config_text = json.dumps(
        config,
        indent=2,
        sort_keys=True,
    ) + "\n"

    manifest.to_csv(
        manifest_path,
        index=False,
        na_rep="",
    )

    config_path.write_text(config_text)
    summary_path.write_text(summary_text)

    log_text = (
        "SCRIPT: scripts/02_create_subject_splits.py\n"
        f"PROJECT_ROOT: {root}\n"
        "STATUS: PASS\n\n"
        + summary_text
    )

    log_path.write_text(log_text)

    hash_targets = [
        manifest_path,
        config_path,
        summary_path,
        log_path,
    ]

    hash_lines = [
        f"{sha256_file(path)}  {path.relative_to(root)}"
        for path in hash_targets
    ]

    hashes_path.write_text("\n".join(hash_lines) + "\n")

    print(summary_text, end="")

    print("\nSPLIT OUTPUTS")

    for path in [
        manifest_path,
        config_path,
        summary_path,
        hashes_path,
        log_path,
    ]:
        print(
            f"{path.relative_to(root)}\t"
            f"{path.stat().st_size} bytes"
        )

    print("\nSUBJECT SPLIT CREATION: PASS")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as error:
        print(
            f"\nSUBJECT SPLIT CREATION: FAIL\n{error}",
            file=sys.stderr,
        )
        sys.exit(1)
