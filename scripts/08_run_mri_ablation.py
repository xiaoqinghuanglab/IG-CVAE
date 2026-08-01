#!/usr/bin/env python3
"""MRI regional-ablation utilities; atlas preparation is affine-aware."""
from __future__ import annotations
import argparse, hashlib, json
from pathlib import Path

import nibabel as nib
from nibabel.processing import resample_from_to
import numpy as np
import pandas as pd
import torch
import torch.nn.functional as F

ROOT = Path(__file__).resolve().parents[1]
SPLIT = ROOT / "Output/Intermediate/Splits/model_split_manifest.csv"
ATLAS_DIR = ROOT / "Resources/Atlas/DesikanKilliany"
ATLAS = ATLAS_DIR / "desikanKillianyMNI152_2mm.nii.gz"
RAW_ALIGNED = ATLAS_DIR / "desikanKilliany_on_MRIgrid_97x115x97.nii.gz"
MODEL_LABELS = ATLAS_DIR / "desikanKilliany_on_modelgrid_48x56x48.npy"
TARGET_SHAPE = (48, 56, 48)

def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda: f.read(1024 * 1024), b""):
            h.update(b)
    return h.hexdigest()

def prepare_atlas():
    if not ATLAS.is_file():
        raise FileNotFoundError(ATLAS)
    d = pd.read_csv(SPLIT)
    dev = d[(d.model == "mri_branch") & (d.role == "development")].copy()
    ref_path = Path(dev[dev.cohort == "ADNI"].iloc[0].mri_path)
    ref = nib.load(str(ref_path))
    atlas = nib.load(str(ATLAS))

    aligned = resample_from_to(atlas, (ref.shape, ref.affine), order=0)
    labels_raw = np.rint(aligned.get_fdata(dtype=np.float32)).astype(np.int16)
    nib.save(nib.Nifti1Image(labels_raw, ref.affine, ref.header), str(RAW_ALIGNED))

    labels_model = F.interpolate(
        torch.from_numpy(labels_raw.astype(np.float32))[None, None],
        size=TARGET_SHAPE, mode="nearest")[0, 0].numpy().astype(np.int16)
    np.save(MODEL_LABELS, labels_model, allow_pickle=False)

    failures = {}
    for cohort in ["ADNI", "AddNeuroMed"]:
        paths = dev.loc[dev.cohort == cohort, "mri_path"].astype(str).unique()
        checked = paths[:3]
        bad = []
        for image_path in checked:
            img = nib.load(image_path)
            if img.shape != ref.shape or not np.allclose(img.affine, ref.affine, atol=1e-4):
                bad.append(image_path)
        failures[cohort] = {"n_development_images": int(len(paths)),
                            "n_geometry_checked": int(len(checked)),
                            "incompatible_images": int(len(bad)),
                            "examples": bad[:3]}

    raw_labels = np.unique(labels_raw)
    model_labels = np.unique(labels_model)
    report = {
        "atlas_source": str(ATLAS.relative_to(ROOT)),
        "atlas_source_sha256": sha256(ATLAS),
        "reference_image": str(ref_path),
        "reference_shape": list(ref.shape),
        "reference_affine": ref.affine.tolist(),
        "atlas_original_shape": list(atlas.shape),
        "atlas_original_affine": atlas.affine.tolist(),
        "aligned_raw_atlas": str(RAW_ALIGNED.relative_to(ROOT)),
        "aligned_raw_atlas_sha256": sha256(RAW_ALIGNED),
        "model_label_array": str(MODEL_LABELS.relative_to(ROOT)),
        "model_shape": list(labels_model.shape),
        "nonzero_labels_raw_grid": int((raw_labels > 0).sum()),
        "nonzero_labels_model_grid": int((model_labels > 0).sum()),
        "development_geometry_qc": failures,
        "compatible_for_shared_ADNI_AddNeuroMed_region_occlusion":
            all(v["incompatible_images"] == 0 for v in failures.values()),
    }
    qc = ROOT / "Output/QC/MRI_ablation"
    qc.mkdir(parents=True, exist_ok=True)
    (qc / "atlas_preparation.json").write_text(json.dumps(report, indent=2))
    print(json.dumps(report, indent=2))


def ablation_metrics(y, probability):
    from sklearn.metrics import average_precision_score, roc_auc_score
    return {
        "roc_auc": float(roc_auc_score(y, probability)),
        "pr_auc": float(average_precision_score(y, probability)),
    }

def validate(cohort):
    atlas_labels, regions = region_mapping()
    rows = development_rows(cohort)
    print(json.dumps({
        "cohort": cohort,
        "development_n": int(len(rows)),
        "AD": int(rows.diagnosis_binary.sum()),
        "CN": int((rows.diagnosis_binary == 0).sum()),
        "fold_counts": {str(k): int(v) for k, v in
                        rows.cv_fold.value_counts().sort_index().items()},
        "atlas_model_grid_shape": list(atlas_labels.shape),
        "atlas_regions_tested": int(len(regions)),
        "method": "held_out_fold_regional_zero_occlusion",
        "locked_test_accessed": False,
    }, indent=2, sort_keys=True))

def make_region_figure(summary, cohort):
    import matplotlib.pyplot as plt

    top = summary.nlargest(20, "delta_roc_auc_mean").sort_values(
        "delta_roc_auc_mean", ascending=True)
    fig, ax = plt.subplots(figsize=(8.2, 6.8))
    ax.barh(top["region_name"], top["delta_roc_auc_mean"],
            xerr=top["delta_roc_auc_sd"], color="#2166AC",
            edgecolor="black", linewidth=0.45, capsize=2.5)
    ax.axvline(0, color="black", linewidth=0.8)
    ax.set_xlabel("Decrease in held-out ROC-AUC after regional occlusion")
    ax.set_ylabel("")
    ax.set_title(f"{cohort}: MRI regional importance (development only)")
    ax.tick_params(axis="y", labelsize=8)
    ax.grid(axis="x", alpha=0.22, linewidth=0.6)
    fig.tight_layout()

    panel = "C" if cohort == "ADNI" else "D"
    fig_dir = ROOT / "Output/Figure/Main/Figure4_components"
    fig_dir.mkdir(parents=True, exist_ok=True)
    stem = fig_dir / f"Figure4_component_{panel}_{cohort}_MRI_region_occlusion_development"
    fig.savefig(str(stem) + ".png", dpi=400, bbox_inches="tight")
    fig.savefig(str(stem) + ".pdf", bbox_inches="tight")
    plt.close(fig)
    return [str(stem.relative_to(ROOT)) + ".png",
            str(stem.relative_to(ROOT)) + ".pdf"]

def run(cohort, device):
    mri = load_mri_module()
    atlas_labels, regions = region_mapping()
    label_codes = [code for code, _ in regions]
    region_names = dict(regions)

    rows = development_rows(cohort)
    cache = ablation_outdir("Intermediate", cohort) / "mri_preprocessing_cache"
    cache.mkdir(parents=True, exist_ok=True)

    all_rows, baseline_rows, provenance = [], [], []
    for fold in range(1, 6):
        validation = rows[rows.cv_fold == fold].copy().reset_index(drop=True)
        images, y, ids = load_images(mri, validation, cache, device)

        if tuple(images.shape[-3:]) != tuple(atlas_labels.shape):
            raise RuntimeError(
                f"{cohort} fold {fold}: image grid {tuple(images.shape[-3:])} "
                f"does not match atlas grid {tuple(atlas_labels.shape)}")
        if set(ids) != set(validation.subject_id.astype(str)):
            raise RuntimeError(f"{cohort} fold {fold}: image/manifest subject mismatch")

        model, model_info = load_fold_model(mri, cohort, fold, device)
        baseline_probability = predict_probabilities(model, images, device)
        baseline = ablation_metrics(y, baseline_probability)
        occluded_probability = regional_occlusion_probabilities(
            model, images, atlas_labels, label_codes, device)

        baseline_rows.append({
            "fold": fold, "n_validation": int(len(y)),
            "AD": int(y.sum()), "CN": int((y == 0).sum()), **baseline})
        provenance.append({"fold": fold, **model_info})

        for i, code in enumerate(label_codes):
            occluded = ablation_metrics(y, occluded_probability[i])
            all_rows.append({
                "fold": fold,
                "label_id": int(code),
                "region_name": region_names[code],
                "n_voxels_model_grid": int((atlas_labels == code).sum()),
                "baseline_roc_auc": baseline["roc_auc"],
                "baseline_pr_auc": baseline["pr_auc"],
                "occluded_roc_auc": occluded["roc_auc"],
                "occluded_pr_auc": occluded["pr_auc"],
                "delta_roc_auc": baseline["roc_auc"] - occluded["roc_auc"],
                "delta_pr_auc": baseline["pr_auc"] - occluded["pr_auc"],
            })
        print(f"{cohort} fold={fold} baseline_AUC={baseline['roc_auc']:.3f} "
              f"baseline_PR={baseline['pr_auc']:.3f}")

    detail = pd.DataFrame(all_rows)
    summary = detail.groupby(["label_id", "region_name"], as_index=False).agg(
        n_voxels_model_grid=("n_voxels_model_grid", "first"),
        delta_roc_auc_mean=("delta_roc_auc", "mean"),
        delta_roc_auc_sd=("delta_roc_auc", "std"),
        delta_pr_auc_mean=("delta_pr_auc", "mean"),
        delta_pr_auc_sd=("delta_pr_auc", "std"),
    )
    positive = detail.assign(positive=detail.delta_roc_auc > 0).groupby(
        ["label_id", "region_name"], as_index=False)["positive"].sum()
    summary = summary.merge(positive, on=["label_id", "region_name"])
    summary = summary.rename(columns={"positive": "positive_auc_folds"}).sort_values(
        ["delta_roc_auc_mean", "delta_pr_auc_mean"], ascending=False).reset_index(drop=True)

    table_dir = ablation_outdir("Table", cohort)
    detail.to_csv(table_dir / "fold_region_occlusion.csv", index=False)
    summary.to_csv(table_dir / "region_occlusion_summary.csv", index=False)
    pd.DataFrame(baseline_rows).to_csv(
        table_dir / "fold_baseline_metrics.csv", index=False)
    summary.head(20).to_csv(table_dir / "top20_regions.csv", index=False)
    (table_dir / "top20_regions.txt").write_text(
        summary.head(20).to_string(index=False, float_format=lambda x: f"{x:.5f}") + "\n")

    figure_components = make_region_figure(summary, cohort)
    qc_dir = ablation_outdir("QC", cohort)
    (qc_dir / "provenance.json").write_text(json.dumps({
        "cohort": cohort,
        "method": "held_out_fold_regional_zero_occlusion",
        "locked_test_accessed": False,
        "atlas_regions_tested": len(regions),
        "fold_model_provenance": provenance,
        "figure_components": figure_components,
    }, indent=2))

    print("\nTOP 10 REGIONS")
    print(summary.head(10).to_string(index=False, float_format=lambda x: f"{x:.5f}"))
    print(f"\nSaved: {table_dir / 'region_occlusion_summary.csv'}")

def main():
    ap = argparse.ArgumentParser(
        description="Development-only MRI regional occlusion analysis.")
    ap.add_argument("command", choices=["PREPARE_ATLAS", "VALIDATE", "RUN"])
    ap.add_argument("--cohort", choices=["ADNI", "AddNeuroMed", "ALL"],
                    default="ALL")
    ap.add_argument("--device", choices=["auto", "cpu", "cuda"], default="auto")
    args = ap.parse_args()

    if args.command == "PREPARE_ATLAS":
        prepare_atlas()
        return

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
        else:
            run(cohort, device)


def load_mri_module():
    import importlib.util, sys
    p = ROOT / "scripts/05_run_mri_branch.py"
    spec = importlib.util.spec_from_file_location("mri_ablation_source", p)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod

def ablation_outdir(kind, cohort):
    p = ROOT / "Output" / kind / "MRI_ablation" / "dual_cohort_independent_v1" / cohort
    p.mkdir(parents=True, exist_ok=True)
    return p

def region_mapping():
    ids = [int(x.strip()) for x in
           (ATLAS_DIR / "desikanKillianyNodeIndex.1D").read_text().splitlines() if x.strip()]
    names = [x.strip() for x in
             (ATLAS_DIR / "desikanKillianyNodeNames.txt").read_text().splitlines() if x.strip()]
    if len(ids) != len(names):
        raise RuntimeError("Atlas index/name mapping length mismatch")
    labels = np.load(MODEL_LABELS, allow_pickle=False).astype(np.int16)
    present = sorted(int(v) for v in np.unique(labels) if v > 0)
    mapping = dict(zip(ids, names))
    missing = [v for v in present if v not in mapping]
    if missing:
        raise RuntimeError(f"Atlas labels missing names: {missing}")
    return labels, [(v, mapping[v]) for v in present]

def development_rows(cohort):
    d = pd.read_csv(SPLIT)
    x = d[(d.cohort == cohort) & (d.model == "mri_branch") &
          (d.role == "development")].copy()
    x["diagnosis_binary"] = pd.to_numeric(x.diagnosis_binary).astype(int)
    x["cv_fold"] = pd.to_numeric(x.cv_fold).astype(int)
    if x.empty or x.subject_id.duplicated().any():
        raise RuntimeError(f"{cohort}: invalid MRI development rows")
    return x.sort_values("subject_id").reset_index(drop=True)

def load_images(mri, rows, cache, device):
    dl = mri.loader(rows, cache, False, device)
    images, labels, ids = [], [], []
    for x, y, sid in dl:
        images.append(x)
        labels.extend(y.numpy().astype(int).tolist())
        ids.extend(sid)
    return torch.cat(images, dim=0), np.asarray(labels), list(ids)

def load_fold_model(mri, cohort, fold, device):
    p = ROOT / "Output/Model/MRI_branch/dual_cohort_independent_v1" / cohort / f"fold_{fold}_best.pt"
    ckpt = torch.load(p, map_location="cpu", weights_only=False)
    model = mri.MRIBranch(float(ckpt["config"]["dropout"])).to(device)
    model.load_state_dict(ckpt["state_dict"])
    model.eval()
    return model, {"path": str(p.relative_to(ROOT)), "sha256": sha256(p),
                   "best_epoch": int(ckpt["best_epoch"])}

def predict_probabilities(model, images, device, batch_size=8):
    out = []
    model.eval()
    with torch.no_grad():
        for start in range(0, len(images), batch_size):
            logits, _ = model(images[start:start+batch_size].to(device, non_blocking=True))
            out.append(torch.sigmoid(logits).cpu().numpy())
    return np.concatenate(out)

def regional_occlusion_probabilities(
        model, images, atlas_labels, label_codes, device, region_chunk=4, batch_size=8):
    """Return AD probabilities after zeroing one atlas region at a time.

    Images are per-volume standardized; zero is therefore the training-independent
    mean-intensity value. This function uses development validation images only.
    """
    all_probabilities = []
    model.eval()
    with torch.no_grad():
        for r0 in range(0, len(label_codes), region_chunk):
            codes = label_codes[r0:r0 + region_chunk]
            masks = torch.from_numpy(
                np.stack([(atlas_labels == code) for code in codes])
            ).bool().to(device)
            chunk_probabilities = []
            for b0 in range(0, len(images), batch_size):
                xb = images[b0:b0 + batch_size].to(device, non_blocking=True)
                # [regions, subjects, channels, D, H, W]
                keep = (~masks).unsqueeze(1).unsqueeze(1).to(dtype=xb.dtype)
                occluded = (xb.unsqueeze(0) * keep).reshape(
                    len(codes) * len(xb), *xb.shape[1:])
                logits, _ = model(occluded)
                chunk_probabilities.append(
                    torch.sigmoid(logits).reshape(len(codes), len(xb)).cpu().numpy())
            all_probabilities.append(np.concatenate(chunk_probabilities, axis=1))
    return np.concatenate(all_probabilities, axis=0)
if __name__ == "__main__":
    main()
