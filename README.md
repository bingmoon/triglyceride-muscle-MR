# Code for "Triglycerides and skeletal muscle: resolving the body mass confound through relative muscle indices and multivariable Mendelian randomization"

This repository contains the analytical R scripts used in the three-stage study by Chen, Tang, and Zhang.

## Overview

- **Stage 1** (`stage1_metabolomics.R`): Untargeted metabolomics analysis of murine disuse atrophy (Metabolomics Workbench ST003803).
- **Stage 2** (`stage2_nhanes.R`): Population‑based cohort analysis using NHANES 2011–2018.
- **Stage 3** (`stage3_mr.R`): Univariable and multivariable Mendelian randomization.

## Requirements

R ≥ 4.2.2 with packages: `nhanesA`, `TwoSampleMR`, `ieugwasr`, `MVMR`, `ropls`, `impute`, `ggplot2`, `ggpubr`, `broom`, `dplyr`, `tidyr`, `ggrepel`, `openxlsx`.

## Usage

1. Set the working directory to the repository root.
2. Run each script in order (`stage1_metabolomics.R` → `stage2_nhanes.R` → `stage3_mr.R`).
3. Output figures and tables are saved to the same folder.

## Data Availability

All datasets are publicly available:
- Metabolomics: [Metabolomics Workbench ST003803](https://doi.org/10.21228/M8S3803)
- NHANES: [CDC NHANES](https://wwwn.cdc.gov/nchs/nhanes/)
- GWAS: [IEU OpenGWAS](https://gwas.mrcieu.ac.uk/)

## Current Submission

This manuscript is currently under review at **PLOS ONE**.

## License

MIT License. See [LICENSE](LICENSE) for details.
