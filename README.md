# Triglycerides and Skeletal Muscle: Resolving the Body Mass Confound through Relative Muscle Indices and Mendelian Randomization

This repository contains the complete analytical R scripts for the three-stage study by Chen, Tang, and Zhang.

---

## Overview

This study investigates whether circulating triglycerides (TG) independently impair skeletal muscle function and mass, after accounting for the mechanical confounding effect of body mass. A triangulation framework was employed:

- **Stage 1** (`stage1_metabolomics.R`): Untargeted metabolomics of murine disuse atrophy (Metabolomics Workbench study ST003803). Characterizes intramuscular lipid mobilization during muscle loss.
- **Stage 2** (`stage2_nhanes.R`): Population-based cross-sectional analysis of NHANES 2011--2018 adults. Examines associations of TG with relative grip strength and DXA-derived relative appendicular skeletal muscle mass.
- **Stage 3** (`stage3_mr.R`): Primary univariable and multivariable Mendelian randomization of TG on grip strength (GWAS: TG ieu-b-111; outcome: ukb-b-10215; BMI: ieu-b-40).
- **Post-rejection sensitivity analyses** (`stage3_post_rejection_analysis.R`): PhenoScanner-based confounder screening, re-analysis after excluding pleiotropic SNPs, biological annotation of top instruments.

---

## Key Findings

- Observational inverse associations between TG and muscle outcomes are predominantly driven by positive confounding from body mass (BMI).
- Univariable MR using 112 well-powered genetic instruments reveals no evidence of a causal effect of TG on grip strength (IVW β = 0.001, P = 0.929).
- After excluding 210 confounder-associated SNPs, the null result persists (β = −0.007, P = 0.780), with no evidence of directional pleiotropy.
- Multivariable MR conditioning on BMI yields similarly null findings.
- Bidirectional MR provides no evidence of reverse causation.

---

## Repository Structure

```
triglyceride-muscle-MR/
├── README.md
├── LICENSE
├── stage1_metabolomics.R
├── stage2_nhanes.R
├── stage3_mr.R
├── stage3_post_rejection_analysis.R
├── references.bib
└── output/                     # Generated figures and tables
    ├── Figure_*.pdf
    └── Table_*.csv
```

---

## Requirements

- **R** ≥ 4.2.2
- **R packages**: `nhanesA`, `TwoSampleMR`, `ieugwasr`, `MVMR`, `ropls`, `impute`, `ggplot2`, `ggpubr`, `broom`, `dplyr`, `tidyr`, `ggrepel`, `writexl`, `readr`

Install missing packages:
```r
install.packages(c("nhanesA", "ggplot2", "ggpubr", "broom", "dplyr", "tidyr", "ggrepel", "writexl", "readr"))
install.packages(c("TwoSampleMR", "ieugwasr", "MVMR", "ropls", "impute"))
```

---

## Usage

1. Clone the repository or download all scripts to a local working directory.
2. Set the working directory in R:
   ```r
   setwd("/path/to/triglyceride-muscle-MR")
   ```
3. Run scripts in order:
   ```r
   source("stage1_metabolomics.R")
   source("stage2_nhanes.R")
   source("stage3_mr.R")
   source("stage3_post_rejection_analysis.R")
   ```
4. Output figures (PDF) and tables (CSV) are saved to the `output/` subdirectory.
5. Combine CSV tables into a single Excel file:
   ```r
   source("merge_tables_to_excel.R")
   ```

---

## Data Availability

All datasets analyzed in this study are publicly accessible:

| Data Source | Access |
|-------------|--------|
| Metabolomics | [Metabolomics Workbench ST003803](https://doi.org/10.21228/M8S3803) |
| NHANES 2011--2018 | [CDC NHANES](https://wwwn.cdc.gov/nchs/nhanes/) |
| TG GWAS (ieu-b-111) | [IEU OpenGWAS](https://gwas.mrcieu.ac.uk/datasets/ieu-b-111/) |
| BMI GWAS (ieu-b-40) | [IEU OpenGWAS](https://gwas.mrcieu.ac.uk/datasets/ieu-b-40/) |
| Grip strength GWAS (ukb-b-10215) | [IEU OpenGWAS](https://gwas.mrcieu.ac.uk/datasets/ukb-b-10215/) |

---

## Reproducibility

All analyses are fully reproducible from the raw data sources listed above. The scripts include:
- Complete data extraction and cleaning pipelines
- Instrument selection and validation procedures
- PhenoScanner-based confounder screening with network-dependent API queries
- All sensitivity analyses (MR-Egger, weighted median, leave-one-out, bidirectional MR)
- Biological annotation of key instrumental SNPs

---

## Current Status

This manuscript has undergone major revision addressing reviewer concerns regarding:
- Power calculation and instrument strength justification
- Detailed rationale for instrument selection thresholds
- Comprehensive confounder screening via PhenoScanner
- Literature-based biological validation of genetic instruments
- Restrained conclusions consistent with the evidence

---

## License

MIT License. See [LICENSE](LICENSE) for details.

---

## Citation

If you use these scripts or data in your own work, please cite:

Chen Y, Tang B, Zhang J. Triglycerides and Skeletal Muscle: Resolving the Body Mass Confound through Relative Muscle Indices and Mendelian Randomization. *(manuscript under review)*.

```bibtex
@article{chen2025tgmuscle,
  title = {Triglycerides and Skeletal Muscle: Resolving the Body Mass Confound through Relative Muscle Indices and Mendelian Randomization},
  author = {Chen, Yanru and Tang, Bing and Zhang, Jie},
  year = {2026}
}
```
