
# Code for "Triglycerides and Skeletal Muscle: Resolving the Body Mass Confound through Relative Muscle Indices and Mendelian Randomization"

This repository contains the complete analytical R scripts for the three‑stage study by Chen, Tang, and Zhang (2026).

## Overview

- **Stage 1** (`stage1_metabolomics.R`): Untargeted metabolomics of murine disuse atrophy (Metabolomics Workbench ST003803).
- **Stage 2** (`stage2_nhanes.R`): Population‑based cross‑sectional analysis of NHANES 2011–2018.
- **Stage 3** (`stage3_mr.R`): Univariable and multivariable Mendelian randomization of triglycerides on grip strength, including confounder screening and sensitivity analyses.

## Key Findings

- The inverse observational association between triglycerides and muscle outcomes is almost entirely driven by positive confounding from body mass (BMI).
- Univariable MR using 112 well‑powered genetic instruments finds no evidence of a causal effect of TG on grip strength (IVW β = 0.001, P = 0.929).
- After excluding 210 confounder‑associated SNPs, the null result remains robust (β = −0.007, P = 0.780), with no directional pleiotropy.
- Multivariable MR conditioning on BMI yields similarly null results, and bidirectional MR shows no reverse causation.

## Requirements

R ≥ 4.2.2 with packages: `nhanesA`, `TwoSampleMR`, `ieugwasr`, `MVMR`, `ropls`, `impute`, `ggplot2`, `ggpubr`, `broom`, `dplyr`, `tidyr`, `ggrepel`.

Install missing packages:

```r
install.packages(c("nhanesA", "ggplot2", "ggpubr", "broom", "dplyr", "tidyr", "ggrepel"))
# For Bioconductor / GitHub packages:
# remotes::install_github("MRCIEU/TwoSampleMR")
# remotes::install_github("MRCIEU/ieugwasr")
```

## Usage

1. Clone or download the repository.
2. Set the working directory to the repository root.
3. Run the scripts in order:

```r
source("stage1_metabolomics.R")
source("stage2_nhanes.R")
source("stage3_mr.R")
```

Figures and tables are saved directly in the working directory.

## Data Availability

All datasets are publicly available:

- Metabolomics: [Metabolomics Workbench ST003803](https://doi.org/10.21228/M8S3803)  
- NHANES 2011–2018: [CDC NHANES](https://wwwn.cdc.gov/nchs/nhanes/)  
- GWAS summary statistics: [IEU OpenGWAS](https://gwas.mrcieu.ac.uk/)  
  - Triglycerides: `ieu-b-111`  
  - Hand grip strength: `ukb-b-10215`  
  - Body mass index: `ieu-b-40`

## License

MIT License. See `LICENSE` for details.

## Citation

If you use this code or data, please cite:

> Chen Y, Tang B, Zhang J. Triglycerides and Skeletal Muscle: Resolving the Body Mass Confound through Relative Muscle Indices and Mendelian Randomization. *Manuscript under review*, 2026.

```bibtex
@article{chen2026tgmuscle,
  title   = {Triglycerides and Skeletal Muscle: Resolving the Body Mass Confound through Relative Muscle Indices and Mendelian Randomization},
  author  = {Chen, Yanru and Tang, Bing and Zhang, Jie},
  year    = {2026},
  note    = {Manuscript under review}
}
```
