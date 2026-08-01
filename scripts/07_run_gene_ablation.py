#!/usr/bin/env python3
"""Development-only permutation importance for the frozen top-500 gene branch."""
from __future__ import annotations
import argparse, hashlib, importlib.util, json, sys
from pathlib import Path

import joblib
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import torch
from sklearn.metrics import average_precision_score, roc_auc_score

ROOT = Path(__file__).resolve().parents[1]
TAG = "dual_cohort_top500_v1"
N_REPEATS = 5
BASE_SEED = 20260718

def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod

GENE = load_module("xig_gene_branch", ROOT / "scripts/04_run_gene_branch.py")

def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda: f.read(1024 * 1024), b""):
            h.update(b)
    return h.hexdigest()

def outdir(kind, cohort):
    p = ROOT / "Output" / kind / "Gene_ablation" / TAG / cohort
    p.mkdir(parents=True, exist_ok=True)
    return p

def device_from_arg(name):
    if name == "auto":
        return torch.device("cuda" if torch.cuda.is_available() else "cpu")
    if name == "cuda" and not torch.cuda.is_available():
        raise RuntimeError("CUDA requested but unavailable")
    return torch.device(name)

def load_fold_models(cohort, fold, device):
    base = ROOT / "Output/Model/Gene_branch" / TAG / cohort / f"fold_{fold}"
    scaler_path = base / "standard_scaler.joblib"
    scaler = joblib.load(scaler_path)
    models, provenance = [], []
    for seed in (1, 2, 3):
        p = base / f"seed_{seed}_gene_encoder_classifier.pt"
        ckpt = torch.load(p, map_location="cpu", weights_only=False)
        model = GENE.GeneBranch(GENE.TrainConfig(**ckpt["config"])).to(device)
        model.load_state_dict(ckpt["state_dict"])
        model.eval()
        models.append(model)
        provenance.append({"path": str(p.relative_to(ROOT)), "sha256": sha256(p),
                           "seed": int(ckpt["seed"])})
    return scaler, models, {"scaler": str(scaler_path.relative_to(ROOT)),
                             "scaler_sha256": sha256(scaler_path),
                             "models": provenance}

def ensemble_probability(models, X, device):
    return np.mean([GENE.predict(m, X.astype(np.float32), device)[0]
                    for m in models], axis=0)

def all_gene_permutations(models, Xval, device, rng):
    """One shuffled validation copy per gene; returns [500, n_validation] probabilities."""
    n, p = Xval.shape
    block = np.broadcast_to(Xval[None, :, :], (p, n, p)).copy()
    for j in range(p):
        block[j, :, j] = rng.permutation(Xval[:, j])
    probability = ensemble_probability(models, block.reshape(p * n, p), device)
    return probability.reshape(p, n)

def save_figure(summary, cohort):
    top = summary.head(20).iloc[::-1]
    fig, ax = plt.subplots(figsize=(7.2, 6.2))
    ax.barh(top["gene"], top["delta_roc_auc_mean"],
            xerr=top["delta_roc_auc_sd"], color="#2A6FBB",
            edgecolor="black", linewidth=0.6, capsize=2)
    ax.axvline(0, color="black", linewidth=0.8)
    ax.set_xlabel("Decrease in ROC-AUC after validation-fold permutation")
    ax.set_ylabel("Gene")
    ax.set_title(f"{cohort}: top gene-branch permutation importance")
    ax.grid(axis="x", alpha=0.25)
    fig.tight_layout()
    d = ROOT / "Output/Figure/Main/Figure4_components"
    d.mkdir(parents=True, exist_ok=True)
    stem = d / f"Figure4_component_{'A' if cohort == 'ADNI' else 'B'}_{cohort}_gene_permutation_importance_development"
    fig.savefig(str(stem) + ".png", dpi=600, bbox_inches="tight")
    fig.savefig(str(stem) + ".pdf", bbox_inches="tight")
    plt.close(fig)

def run(cohort, device):
    records = GENE.development_rows(cohort)
    matrix_path, order_path = GENE.paths_for(cohort)
    genes = GENE.read_gene_order(order_path)
    X, y, ids, folds = GENE.read_matrix_for_ids(matrix_path, genes, records)

    all_rows, fold_baselines, provenance = [], [], {}
    for fold in range(1, 6):
        va = folds == fold
        scaler, models, prov = load_fold_models(cohort, fold, device)
        Xval = scaler.transform(X[va]).astype(np.float32)
        yval = y[va]
        baseline = ensemble_probability(models, Xval, device)
        baseline_auc = float(roc_auc_score(yval, baseline))
        baseline_pr = float(average_precision_score(yval, baseline))
        fold_baselines.append({"fold": fold, "n_validation": int(va.sum()),
                               "roc_auc": baseline_auc, "pr_auc": baseline_pr})
        provenance[str(fold)] = prov

        auc_drops, pr_drops = [], []
        for repeat in range(N_REPEATS):
            rng = np.random.default_rng(BASE_SEED + 10000 * fold + repeat)
            permuted = all_gene_permutations(models, Xval, device, rng)
            auc_drops.append([
                baseline_auc - roc_auc_score(yval, permuted[j])
                for j in range(len(genes))])
            pr_drops.append([
                baseline_pr - average_precision_score(yval, permuted[j])
                for j in range(len(genes))])

        fold_df = pd.DataFrame({
            "cohort": cohort, "fold": fold, "gene": genes,
            "baseline_roc_auc": baseline_auc, "baseline_pr_auc": baseline_pr,
            "delta_roc_auc": np.mean(auc_drops, axis=0),
            "delta_pr_auc": np.mean(pr_drops, axis=0),
        })
        all_rows.append(fold_df)
        print(f"{cohort} fold={fold} baseline_AUC={baseline_auc:.3f} baseline_PR={baseline_pr:.3f}")

    detailed = pd.concat(all_rows, ignore_index=True)
    summary = detailed.groupby("gene", as_index=False).agg(
        delta_roc_auc_mean=("delta_roc_auc", "mean"),
        delta_roc_auc_sd=("delta_roc_auc", "std"),
        delta_pr_auc_mean=("delta_pr_auc", "mean"),
        delta_pr_auc_sd=("delta_pr_auc", "std"),
        positive_auc_folds=("delta_roc_auc", lambda x: int((x > 0).sum())),
    ).sort_values(["delta_roc_auc_mean", "delta_pr_auc_mean"],
                  ascending=False).reset_index(drop=True)

    table = outdir("Table", cohort)
    qc = outdir("QC", cohort)
    detailed.to_csv(table / "fold_permutation_importance.csv", index=False)
    summary.to_csv(table / "gene_permutation_importance_summary.csv", index=False)
    (table / "top20_genes.txt").write_text(
        summary.head(20).to_csv(index=False), encoding="utf-8")
    with (qc / "provenance.json").open("w") as f:
        json.dump({
            "cohort": cohort, "method": "validation_fold_permutation_importance",
            "permutation_repeats_per_gene_per_fold": N_REPEATS,
            "n_genes": len(genes), "fold_baselines": fold_baselines,
            "frozen_model_provenance": provenance,
            "locked_test_accessed": False,
        }, f, indent=2, sort_keys=True)

    save_figure(summary, cohort)
    print("\nTOP 10 GENES")
    print(summary.head(10).to_string(index=False, float_format=lambda x: f"{x:.5f}"))
    print(f"\nSaved: {table / 'gene_permutation_importance_summary.csv'}")

def validate(cohort):
    records = GENE.development_rows(cohort)
    matrix_path, order_path = GENE.paths_for(cohort)
    X, y, _, folds = GENE.read_matrix_for_ids(
        matrix_path, GENE.read_gene_order(order_path), records)
    print(json.dumps({
        "cohort": cohort, "development_n": int(len(y)),
        "AD": int(y.sum()), "CN": int((y == 0).sum()),
        "matrix_shape": list(X.shape),
        "fold_counts": {str(k): int((folds == k).sum()) for k in range(1, 6)},
        "method": "five_repeat_validation_fold_permutation_importance",
        "locked_test_accessed": False,
    }, indent=2, sort_keys=True))

def main():
    ap = argparse.ArgumentParser(description="Development-only gene permutation ablation.")
    ap.add_argument("command", choices=["VALIDATE", "RUN"])
    ap.add_argument("--cohort", choices=["ADNI", "AddNeuroMed", "ALL"], default="ALL")
    ap.add_argument("--device", choices=["auto", "cpu", "cuda"], default="auto")
    args = ap.parse_args()
    device = device_from_arg(args.device)
    cohorts = ["ADNI", "AddNeuroMed"] if args.cohort == "ALL" else [args.cohort]
    for cohort in cohorts:
        validate(cohort) if args.command == "VALIDATE" else run(cohort, device)

if __name__ == "__main__":
    main()
