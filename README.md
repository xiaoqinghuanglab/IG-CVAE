# IG-CVAE

IG-CVAE is an MRI-conditioned conditional variational autoencoder for learning gene-feature representations associated with Alzheimer's disease. ADNI and AddNeuroMed are analyzed independently. The gene and MRI branches provide representations to the CVAE; the remaining scripts evaluate the models and run sensitivity and regional analyses.

## Repository structure

```text
IG-CVAE/
├── README.md
├── .gitignore
├── environment/requirements.txt         Python dependencies for Scripts 01–12
├── data/
│   ├── ADNI/                       Placeholder for private inputs
│   └── AddNeuroMed/                Placeholder for private inputs
├── Resources/Atlas/DesikanKilliany/  Atlas images, labels, and region names
├── scripts/
│   ├── 00_install_r_deg_packages.R
│   ├── 01_generate_table1.py
│   ├── 02_create_subject_splits.py
│   ├── 03_run_deg_pipeline.R
│   ├── 04_run_gene_branch.py
│   ├── 05_run_mri_branch.py
│   ├── 06_run_cvae.py
│   ├── 07_run_gene_ablation.py
│   ├── 08_run_mri_ablation.py
│   ├── 09_run_cvae_ablation.py
│   ├── 10_mlp.py
│   ├── 11_covadj_cvae.py
│   └── 12_pair_shuffle_test.py
└── Output/Table/                  
```

`Output/Table/` contains results for the gene and MRI branches, CVAE, regional and gene ablations, MLP, covariate-adjusted CVAE, pair-shuffle control, and supplementary tables. Participant-level expression and DEG design matrices in that tree are explicitly excluded by `.gitignore`. Other local output directories contain checkpoints, figures, logs, and intermediate representations and are not published here.

## Setup and inputs

Run commands from the repository root. Obtain ADNI and AddNeuroMed data independently under their access requirements. The tracked `data/` directories are placeholders. Some scripts use the path `data/AdNeuroMed/` while the tracked placeholder is `data/AddNeuroMed/`; ensure the local input directory matches the paths in the scripts before running.

Use Python 3.12 for the versions listed in `environment/requirements.txt`. The file covers the direct third-party Python imports in Scripts 01–12, including PyTorch (`torch`) for training and NiBabel (`nibabel`) for MRI handling. Install a CPU or CUDA build of PyTorch suitable for your machine using the [official installation selector](https://pytorch.org/get-started/locally/), then install the remaining Python packages:

```bash
python -m pip install -r environment/requirements.txt
```

The requirement versions specify the package environment; they do not include the private input data or pretrained model checkpoints. Script 03 uses R; install its DEG packages into the project library with:

```bash
Rscript scripts/00_install_r_deg_packages.R
```

The following preparatory steps require private participant inputs. They create Table 1, the subject split manifest, and cohort-specific top-500 gene inputs:

```bash
python scripts/01_generate_table1.py
python scripts/02_create_subject_splits.py --root .
Rscript scripts/03_run_deg_pipeline.R ALL
```

The split manifest is written to `Output/Intermediate/Splits/model_split_manifest.csv`; selected gene matrices and orders are written to `Output/Model_input/Gene/`. Script 02 refuses to overwrite existing split outputs unless `--force` is supplied. Preserve the original split when reproducing downstream results.

## Model and analysis commands

Run `VALIDATE` before each stage. The examples process both cohorts independently with `--cohort ALL`; use `ADNI` or `AddNeuroMed` for one cohort. `DEVELOP` fits development models, and `RUN` performs development-only analyses. Downstream stages require upstream locally generated inputs and checkpoints.

### 1. Gene-branch training

```bash
python scripts/04_run_gene_branch.py VALIDATE --cohort ALL
python scripts/04_run_gene_branch.py DEVELOP --cohort ALL --device auto
python scripts/04_run_gene_branch.py SUMMARY
```

The gene branch writes checkpoints under `Output/Model/Gene_branch/dual_cohort_top500_v1/`. After development review, locked-test evaluation is available with `python scripts/04_run_gene_branch.py TEST --cohort ALL --allow-locked-test`.

### 2. MRI-branch training

```bash
python scripts/05_run_mri_branch.py VALIDATE --cohort ALL
python scripts/05_run_mri_branch.py DEVELOP --cohort ALL --device auto
```

The MRI branch writes checkpoints under `Output/Model/MRI_branch/dual_cohort_independent_v1/`. After development review, use `python scripts/05_run_mri_branch.py TEST --cohort ALL --allow-locked-test` for the locked test.

### 3. MRI-conditioned CVAE training

Requires the gene and MRI branch checkpoints and paired development subjects.

```bash
python scripts/06_run_cvae.py VALIDATE --cohort ALL
python scripts/06_run_cvae.py DEVELOP --cohort ALL --device auto
```

CVAE checkpoints are written under `Output/Model/CVAE/dual_cohort_independent_v1/`. After development review, use `python scripts/06_run_cvae.py TEST --cohort ALL --allow-locked-test` for the locked test.

### 4. Gene permutation analysis

Uses frozen gene-branch fold models and development validation folds.

```bash
python scripts/07_run_gene_ablation.py VALIDATE --cohort ALL
python scripts/07_run_gene_ablation.py RUN --cohort ALL --device auto
```

### 5. MRI regional zero-occlusion analysis

Uses the atlas and frozen MRI-branch fold models. `PREPARE_ATLAS` is available for a new local MRI grid; it rewrites aligned atlas files and is unnecessary merely to use the tracked atlas files.

```bash
python scripts/08_run_mri_ablation.py VALIDATE --cohort ALL
python scripts/08_run_mri_ablation.py RUN --cohort ALL --device auto
```

### 6. Region–gene-feature analysis

Uses frozen MRI and CVAE fold models, atlas labels, and the gene permutation ranking. This version of Script 09 writes region and gene-effect tables and does not generate a plot.

```bash
python scripts/09_run_cvae_ablation.py VALIDATE --cohort ALL
python scripts/09_run_cvae_ablation.py RUN --cohort ALL --device auto
```

### 7. MLP baseline

Uses paired participants and frozen branch representations for comparison with the CVAE.

```bash
python scripts/10_mlp.py VALIDATE --cohort ALL
python scripts/10_mlp.py DEVELOP --cohort ALL --device auto
```

After development review, use `python scripts/10_mlp.py TEST --cohort ALL --allow-locked-test` for the locked test.

### 8. Covariate-adjusted CVAE

Uses the original branch checkpoints and adds age, sex, and APOE4 covariates.

```bash
python scripts/11_covadj_cvae.py VALIDATE --cohort ALL
python scripts/11_covadj_cvae.py DEVELOP --cohort ALL --device auto
```

After development review, use `python scripts/11_covadj_cvae.py TEST --cohort ALL --allow-locked-test` for the locked test.

### 9. Pair-shuffle control

Shuffles training targets among subjects with the same diagnosis. The default is 200 shuffles per cohort; `RUN` does not access the locked test.

```bash
python scripts/12_pair_shuffle_test.py VALIDATE --cohort ALL
python scripts/12_pair_shuffle_test.py RUN --cohort ALL --device auto --n-shuffles 200
```

## Local path requirements for Scripts 07–09

Scripts 07–09 currently refer to branch checkpoints under `Output/main_model/Model/`, although Scripts 04–06 write under `Output/Model/`. Their `RUN` commands therefore require the expected checkpoint paths to be present or adjusted locally. Script 09 additionally sets an absolute `ROOT` at the top of its source; change that setting to the local project root if necessary. These paths are not supplied by cloning the repository.

## Data availability and interpretation

ADNI and AddNeuroMed participant-level data, trained checkpoints, and intermediate representations are not distributed here. The published `Output/Table/` files contain analysis results, but private inputs and locally fitted models are required for a complete rerun.

IG-CVAE generates MRI-conditioned gene-feature representations. It does not reconstruct raw transcript abundance, measure regional brain gene expression, or establish causal MRI–gene relationships.


