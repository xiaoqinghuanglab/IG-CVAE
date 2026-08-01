#!/usr/bin/env python3
"""Development-only CVAE regional ablation and projected gene-feature analysis."""
import argparse
import hashlib
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import torch

ROOT = Path(__file__).resolve().parents[1]
ATLAS = ROOT / "Resources/Atlas/DesikanKilliany/desikanKilliany_on_modelgrid_48x56x48.npy"
ATLAS_DIR = ROOT / "Resources/Atlas/DesikanKilliany"
CVAE_ROOT = ROOT / "Output/Model/CVAE/dual_cohort_independent_v1"
MRI_ROOT = ROOT / "Output/Model/MRI_branch/dual_cohort_independent_v1"
GENE_ABLATION_ROOT = ROOT / "Output/Table/Gene_ablation/dual_cohort_top500_v1"
N_FOLDS = 5
TOP_GENES = 20

def load_local_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module

CVAE = load_local_module("cvae_stage09", ROOT / "scripts/06_run_cvae.py")
MRI = load_local_module("mri_stage09", ROOT / "scripts/05_run_mri_branch.py")

def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()

def outdir(kind, cohort):
    p = ROOT / "Output" / kind / "CVAE_ablation" / "dual_cohort_independent_v1" / cohort
    p.mkdir(parents=True, exist_ok=True)
    return p

def atlas_regions():
    labels = np.load(ATLAS, allow_pickle=False).astype(np.int16)
    ids = [int(x.strip()) for x in
           (ATLAS_DIR / "desikanKillianyNodeIndex.1D").read_text().splitlines()
           if x.strip()]
    names = [x.strip() for x in
             (ATLAS_DIR / "desikanKillianyNodeNames.txt").read_text().splitlines()
             if x.strip()]
    mapping = dict(zip(ids, names))
    present = sorted(int(x) for x in np.unique(labels) if x > 0)
    if len(ids) != len(names) or any(x not in mapping for x in present):
        raise RuntimeError("Invalid atlas label-to-name mapping")
    return labels, present, mapping

def top_gene_indices(cohort):
    table = pd.read_csv(
        GENE_ABLATION_ROOT / cohort / "gene_permutation_importance_summary.csv")
    genes = table.head(TOP_GENES)["gene"].astype(str).tolist()
    _, order_path = CVAE.GENE.paths_for(cohort)
    order = CVAE.GENE.read_gene_order(order_path)
    position = {gene: i for i, gene in enumerate(order)}
    missing = [gene for gene in genes if gene not in position]
    if missing:
        raise RuntimeError(f"{cohort}: ablation genes absent from top-500 order: {missing}")
    return np.asarray([position[g] for g in genes], dtype=int), genes

def paired_images(rows, cohort, device):
    cache = outdir("Intermediate", cohort) / "mri_preprocessing_cache"
    cache.mkdir(parents=True, exist_ok=True)
    by_id = {}
    for xb, _, subject_ids in MRI.loader(rows, cache, False, device):
        for i, sid in enumerate(subject_ids):
            by_id[str(sid)] = xb[i].cpu()
    expected = rows.subject_id.astype(str).tolist()
    if set(by_id) != set(expected):
        raise RuntimeError(f"{cohort}: paired MRI/manifest IDs do not match")
    return torch.stack([by_id[sid] for sid in expected], dim=0)

def load_mri_fold(cohort, fold, device):
    path = MRI_ROOT / cohort / f"fold_{fold}_best.pt"
    ckpt = torch.load(path, map_location="cpu", weights_only=False)
    model = MRI.MRIBranch(float(ckpt["config"]["dropout"])).to(device)
    model.load_state_dict(ckpt["state_dict"])
    model.eval()
    return model, ckpt, path

def load_cvae_fold(cohort, fold, device):
    path = CVAE_ROOT / cohort / f"fold_{fold}_best.pt"
    ckpt = torch.load(path, map_location="cpu", weights_only=False)
    model = CVAE.ConditionalVAE().to(device)
    model.load_state_dict(ckpt["state_dict"])
    model.eval()
    return model, ckpt, path

def mfv_from_images(model, images, device, batch_size=8):
    out = []
    with torch.no_grad():
        for start in range(0, len(images), batch_size):
            _, mfv = model(images[start:start + batch_size].to(device, non_blocking=True))
            out.append(mfv.cpu().numpy())
    return np.concatenate(out).astype(np.float32)

def generated_gfv(model, mfv_raw, ckpt, device):
    mfv = (mfv_raw - ckpt["mfv_mean"]) / ckpt["mfv_scale"]
    generated_scaled = CVAE.predict_generated(
        model, mfv.astype(np.float32), device)
    generated_raw = generated_scaled * ckpt["gfv_scale"] + ckpt["gfv_mean"]
    return generated_scaled.astype(np.float32), generated_raw.astype(np.float32)

def occluded_mfvs(model, images, atlas_labels, label_ids, device,
                  region_chunk=4, batch_size=8):
    all_regions = []
    model.eval()
    with torch.no_grad():
        for r0 in range(0, len(label_ids), region_chunk):
            codes = label_ids[r0:r0 + region_chunk]
            masks = torch.from_numpy(np.stack(
                [(atlas_labels == code) for code in codes])).bool().to(device)
            batches = []
            for b0 in range(0, len(images), batch_size):
                xb = images[b0:b0 + batch_size].to(device, non_blocking=True)
                keep = (~masks).unsqueeze(1).unsqueeze(1).to(dtype=xb.dtype)
                occluded = (xb.unsqueeze(0) * keep).reshape(
                    len(codes) * len(xb), *xb.shape[1:])
                _, mfv = model(occluded)
                batches.append(mfv.reshape(len(codes), len(xb), -1).cpu().numpy())
            all_regions.append(np.concatenate(batches, axis=1))
    return np.concatenate(all_regions, axis=0).astype(np.float32)

def validate(cohort):
    meta = CVAE.cvae_rows(cohort)
    X = CVAE.load_gene_matrix(cohort, meta)
    labels, region_ids, _ = atlas_regions()
    gene_indices, genes = top_gene_indices(cohort)
    print(json.dumps({
        "cohort": cohort,
        "development_n": int(len(meta)),
        "AD": int(meta.diagnosis_binary.sum()),
        "CN": int((meta.diagnosis_binary == 0).sum()),
        "fold_counts": {str(k): int(v) for k, v in
                        meta.cv_fold.value_counts().sort_index().items()},
        "gene_matrix_shape": list(X.shape),
        "top_gene_projection_n": int(len(gene_indices)),
        "top_gene_projection_genes": genes,
        "atlas_model_grid_shape": list(labels.shape),
        "atlas_regions_tested": int(len(region_ids)),
        "method": "held_out_fold_CVAE_region_occlusion_with_train_fold_GFV_to_gene_projection",
        "locked_test_accessed": False,
    }, indent=2, sort_keys=True))


def make_figures(region_summary, effect_mean, gene_names, region_ids, region_names, cohort):
    import matplotlib.pyplot as plt

    panel_region = "A" if cohort == "ADNI" else "B"
    panel_matrix = "C" if cohort == "ADNI" else "D"
    figure_dir = ROOT / "Output/Figure/Main/Figure5_components"
    figure_dir.mkdir(parents=True, exist_ok=True)

    top = region_summary.nlargest(20, "generated_gfv_shift_mean").sort_values(
        "generated_gfv_shift_mean", ascending=True)
    fig, ax = plt.subplots(figsize=(8.2, 6.8))
    ax.barh(top.region_name, top.generated_gfv_shift_mean,
            xerr=top.generated_gfv_shift_sd, color="#762A83",
            edgecolor="black", linewidth=0.45, capsize=2.5)
    ax.set_xlabel("Mean shift in generated GFV (scaled-feature L2 distance)")
    ax.set_ylabel("")
    ax.set_title(f"{cohort}: CVAE regional occlusion (development only)")
    ax.tick_params(axis="y", labelsize=8)
    ax.grid(axis="x", alpha=0.22, linewidth=0.6)
    fig.tight_layout()
    stem = figure_dir / f"Figure5_component_{panel_region}_{cohort}_CVAE_regional_GFV_shift_development"
    fig.savefig(str(stem) + ".png", dpi=400, bbox_inches="tight")
    fig.savefig(str(stem) + ".pdf", bbox_inches="tight")
    plt.close(fig)

    selected_regions = np.argsort(region_summary.generated_gfv_shift_mean.to_numpy())[::-1][:12]
    selected_genes = np.argsort(effect_mean[selected_regions].mean(axis=0))[::-1][:12]
    heat = effect_mean[np.ix_(selected_regions, selected_genes)]

    fig, ax = plt.subplots(figsize=(9.0, 7.0))
    image = ax.imshow(heat, aspect="auto", cmap="magma")
    ax.set_xticks(range(len(selected_genes)))
    ax.set_xticklabels([gene_names[i] for i in selected_genes], rotation=55,
                       ha="right", fontsize=8)
    ax.set_yticks(range(len(selected_regions)))
    ax.set_yticklabels([region_names[region_ids[i]] for i in selected_regions],
                       fontsize=8)
    ax.set_title(f"{cohort}: projected gene-feature shifts after regional occlusion")
    cb = fig.colorbar(image, ax=ax, fraction=0.035, pad=0.03)
    cb.set_label("Mean absolute projected gene-feature shift (z units)", fontsize=9)
    fig.tight_layout()
    stem2 = figure_dir / f"Figure5_component_{panel_matrix}_{cohort}_CVAE_gene_region_projection_development"
    fig.savefig(str(stem2) + ".png", dpi=400, bbox_inches="tight")
    fig.savefig(str(stem2) + ".pdf", bbox_inches="tight")
    plt.close(fig)
    return [str(stem.relative_to(ROOT)) + ".png",
            str(stem.relative_to(ROOT)) + ".pdf",
            str(stem2.relative_to(ROOT)) + ".png",
            str(stem2.relative_to(ROOT)) + ".pdf"]

def run(cohort, device):
    from sklearn.linear_model import Ridge
    from sklearn.preprocessing import StandardScaler

    meta = CVAE.cvae_rows(cohort)
    X = CVAE.load_gene_matrix(cohort, meta)
    atlas_labels, region_ids, region_names = atlas_regions()
    gene_indices, gene_names = top_gene_indices(cohort)

    fold_effects, fold_signed, fold_gfv_shifts, provenance = [], [], [], []
    for fold in range(1, N_FOLDS + 1):
        train = meta.cv_fold.to_numpy() != fold
        valid = meta.cv_fold.to_numpy() == fold
        val_meta = meta.loc[valid].reset_index(drop=True)

        gfv_all, gene_info = CVAE.extract_gene_features(cohort, fold, X, device)
        gfv_scaler = StandardScaler().fit(gfv_all[train])
        gene_scaler = StandardScaler().fit(X[train][:, gene_indices])
        projection = Ridge(alpha=1.0).fit(
            gfv_scaler.transform(gfv_all[train]),
            gene_scaler.transform(X[train][:, gene_indices]))

        images = paired_images(val_meta, cohort, device)
        if tuple(images.shape[-3:]) != tuple(atlas_labels.shape):
            raise RuntimeError(f"{cohort} fold {fold}: MRI/atlas grid mismatch")

        mri_model, _, mri_path = load_mri_fold(cohort, fold, device)
        cvae_model, cvae_ckpt, cvae_path = load_cvae_fold(cohort, fold, device)

        baseline_mfv = mfv_from_images(mri_model, images, device)
        baseline_scaled, baseline_raw = generated_gfv(
            cvae_model, baseline_mfv, cvae_ckpt, device)

        occluded_mfv = occluded_mfvs(
            mri_model, images, atlas_labels, region_ids, device)
        n_regions, n_subjects, _ = occluded_mfv.shape
        occluded_scaled, occluded_raw = generated_gfv(
            cvae_model, occluded_mfv.reshape(-1, 128), cvae_ckpt, device)
        occluded_scaled = occluded_scaled.reshape(n_regions, n_subjects, 128)
        occluded_raw = occluded_raw.reshape(n_regions, n_subjects, 128)

        baseline_gene = projection.predict(gfv_scaler.transform(baseline_raw))
        occluded_gene = projection.predict(gfv_scaler.transform(
            occluded_raw.reshape(-1, 128))).reshape(n_regions, n_subjects, -1)
        delta_gene = occluded_gene - baseline_gene[None, :, :]

        fold_effects.append(np.abs(delta_gene).mean(axis=1))
        fold_signed.append(delta_gene.mean(axis=1))
        fold_gfv_shifts.append(np.linalg.norm(
            occluded_scaled - baseline_scaled[None, :, :], axis=2).mean(axis=1))
        provenance.append({
            "fold": fold,
            "cvae_checkpoint": str(cvae_path.relative_to(ROOT)),
            "cvae_checkpoint_sha256": sha256(cvae_path),
            "mri_checkpoint": str(mri_path.relative_to(ROOT)),
            "mri_checkpoint_sha256": sha256(mri_path),
            "projection": "Ridge(alpha=1.0), fitted on CVAE-fold training subjects only",
            **gene_info,
        })
        print(f"{cohort} fold={fold} paired_validation_n={len(val_meta)} "
              f"mean_GFV_shift={fold_gfv_shifts[-1].mean():.5f}")

    effects = np.stack(fold_effects)
    signed = np.stack(fold_signed)
    gfv_shifts = np.stack(fold_gfv_shifts)
    effect_mean, effect_sd = effects.mean(axis=0), effects.std(axis=0, ddof=1)
    signed_mean = signed.mean(axis=0)
    sign_consistency = np.maximum((signed > 0).sum(axis=0), (signed < 0).sum(axis=0))

    region_summary = pd.DataFrame({
        "label_id": region_ids,
        "region_name": [region_names[x] for x in region_ids],
        "n_voxels_model_grid": [int((atlas_labels == x).sum()) for x in region_ids],
        "generated_gfv_shift_mean": gfv_shifts.mean(axis=0),
        "generated_gfv_shift_sd": gfv_shifts.std(axis=0, ddof=1),
        "mean_projected_gene_shift_z": effect_mean.mean(axis=1),
    }).sort_values("generated_gfv_shift_mean", ascending=False).reset_index(drop=True)

    edges = []
    for r, region_id in enumerate(region_ids):
        for g, gene in enumerate(gene_names):
            edges.append({
                "label_id": region_id,
                "region_name": region_names[region_id],
                "gene": gene,
                "mean_absolute_projected_gene_shift_z": effect_mean[r, g],
                "sd_absolute_projected_gene_shift_z": effect_sd[r, g],
                "mean_signed_projected_gene_shift_z": signed_mean[r, g],
                "direction_consistent_folds": int(sign_consistency[r, g]),
            })
    edge_table = pd.DataFrame(edges).sort_values(
        "mean_absolute_projected_gene_shift_z", ascending=False).reset_index(drop=True)

    matrix = pd.DataFrame(effect_mean, columns=gene_names)
    matrix.insert(0, "region_name", [region_names[x] for x in region_ids])
    matrix.insert(0, "label_id", region_ids)

    table_dir = outdir("Table", cohort)
    region_summary.to_csv(table_dir / "region_generated_GFV_shift_summary.csv", index=False)
    matrix.to_csv(table_dir / "region_gene_projected_shift_matrix.csv", index=False)
    edge_table.to_csv(table_dir / "gene_region_projected_effects.csv", index=False)
    edge_table.head(30).to_csv(
        table_dir / "top30_gene_region_projected_effects.csv", index=False)
    np.savez_compressed(
        outdir("Intermediate", cohort) / "fold_gene_region_projected_effects.npz",
        absolute_effect=effects, signed_effect=signed, generated_gfv_shift=gfv_shifts,
        label_ids=np.asarray(region_ids), gene_names=np.asarray(gene_names))

    components = make_figures(
        region_summary, effect_mean, gene_names, region_ids, region_names, cohort)
    outdir("QC", cohort).joinpath("provenance.json").write_text(json.dumps({
        "cohort": cohort,
        "method": "held_out_fold_CVAE_region_occlusion_with_train_fold_GFV_to_gene_projection",
        "locked_test_accessed": False,
        "paired_development_n": int(len(meta)),
        "atlas_regions_tested": int(len(region_ids)),
        "genes_projected": gene_names,
        "fold_provenance": provenance,
        "figure_components": components,
    }, indent=2))

    print("\nTOP 10 PROJECTED GENE-REGION EFFECTS")
    print(edge_table.head(10).to_string(index=False, float_format=lambda x: f"{x:.5f}"))
    print(f"\nSaved: {table_dir / 'gene_region_projected_effects.csv'}")

def main():
    ap = argparse.ArgumentParser(
        description="Development-only CVAE regional occlusion and projected gene-feature analysis.")
    ap.add_argument("command", choices=["VALIDATE", "RUN"])
    ap.add_argument("--cohort", choices=["ADNI", "AddNeuroMed", "ALL"], default="ALL")
    ap.add_argument("--device", choices=["auto", "cpu", "cuda"], default="auto")
    args = ap.parse_args()

    device = torch.device("cuda" if args.device == "auto" and torch.cuda.is_available()
                          else ("cpu" if args.device == "auto" else args.device))
    if device.type == "cuda" and not torch.cuda.is_available():
        raise RuntimeError("CUDA requested but unavailable.")
    cohorts = ["ADNI", "AddNeuroMed"] if args.cohort == "ALL" else [args.cohort]
    for cohort in cohorts:
        if args.command == "VALIDATE":
            validate(cohort)
        else:
            run(cohort, device)

if __name__ == "__main__":
    main()
