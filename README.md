# XIG-CVAE

XIG-CVAE is an MRI-conditioned conditional variational autoencoder framework for learning gene-feature representations associated with Alzheimer's disease.

The pipeline independently analyzes ADNI and AddNeuroMed and contains:

1. Cohort summary generation
2. Deterministic train-validation and test splitting
3. Differential-expression analysis
4. Gene-branch training
5. MRI-branch training
6. MRI-conditioned CVAE training
7. Gene permutation analysis
8. MRI regional zero-occlusion analysis
9. Region–gene-feature analysis
10. Supplementary-table export

## Repository structure

```text
scripts/       Executable analysis scripts
Resources/     Desikan–Killiany atlas resources used for regional analyses
environment/   Software dependency specifications
```

Participant-level data, trained models, generated outputs, logs, and local software environments are not included.

## Analysis workflow

The principal scripts are intended to be run in the following order:

1. `01_generate_table1.py` — generate cohort summary tables
2. `02_create_subject_splits.py` — create deterministic train-validation and test assignments
3. `03_run_deg_pipeline.R` — perform cohort-specific differential-expression analysis
4. `04_run_gene_branch.py` — train and evaluate the Gene branch
5. `05_run_mri_branch.py` — train and evaluate the MRI branch
6. `06_run_cvae.py` — train and evaluate the MRI-conditioned CVAE
7. `07_run_gene_ablation.py` — perform gene permutation analysis
8. `08_run_mri_ablation.py` — perform regional MRI zero-occlusion analysis
9. `09_run_cvae_ablation.py` — evaluate region-linked changes in generated gene-feature representations
10. `17_export_full_ablation_supplementary_tables.py` — export supplementary ablation tables

The R package installation helper is:

```text
scripts/00_install_r_deg_packages.R
```

## Atlas resources

The repository includes Desikan–Killiany atlas files used by the MRI and CVAE regional analyses:

```text
Resources/Atlas/DesikanKilliany/
```

These atlas files are analysis resources and do not contain participant data.

## Software environment

The currently available Python dependency specification is:

```text
environment/requirements_gene_branch.txt
```

This file was originally prepared for the Gene branch and should not be interpreted as a complete environment specification for every pipeline component.

The differential-expression scripts require R and the relevant Bioconductor packages, including `limma`.

## Data availability

ADNI and AddNeuroMed participant data are governed by their respective access and data-use requirements and are therefore not distributed through this repository.

Users must obtain authorized access independently and organize the required input files locally.

## Generated outputs

Generated tables, trained models, intermediate representations, figures, and logs are written under local output directories such as:

```text
Output/
```

These files are excluded from version control.

## Interpretation

XIG-CVAE generates MRI-conditioned gene-feature representations. It does not reconstruct raw transcript abundance, measure regional brain gene expression, or establish causal MRI–gene relationships.
