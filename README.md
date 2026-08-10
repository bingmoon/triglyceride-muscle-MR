# Code for "Triglycerides and Skeletal Muscle: Resolving the Body Mass Confound through Relative Muscle Indices and Mendelian Randomization"

This repository contains the complete analytical R scripts and the metabolomics dataset for the three‑stage study by Chen, Tang, and Zhang.

## Repository structure


├── README.md
├── LICENSE
├── stage1_metabolomics.R
├── stage2_nhanes.R
├── stage3_mr.R
├── stage_supplementary_analyses.R
└── MSdata_ST003803_1.txt
```

## Stage overview

- **Stage 1** (`stage1_metabolomics.R`): Untargeted metabolomics of murine disuse atrophy (Metabolomics Workbench ST003803). The script reads `MSdata_ST003803_1.txt` as input.
- **Stage 2** (`stage2_nhanes.R`): Population‑based cross‑sectional analysis of NHANES 2011–2018.
- **Stage 3** (`stage3_mr.R`): Univariable and multivariable Mendelian randomization of triglycerides on grip strength, including instrument selection, PhenoScanner‑based confounder screening, sensitivity analyses, and biological annotation.
- **Supplementary analyses** (`stage_supplementary_analyses.R`): Additional analyses addressing methodological reviewer comments, including:
  - PLS-DA model fit indices (R²Y, Q²)
  - NHANES survey-weighted regression with complex sampling design
  - MR-PRESSO pleiotropy assessment
  - Conditional F-statistic documentation for MVMR

## Key findings

- The observational inverse association between triglycerides and muscle outcomes is largely driven by positive confounding from body mass (BMI).
- In survey-weighted NHANES analyses, TG was inversely associated with relative grip strength (β = −0.060, P = 3.93×10⁻¹⁰) and relative ASM (β = −0.012, P = 2.19×10⁻²¹).
- Univariable MR with 112 genetic instruments shows no evidence of a causal effect of TG on grip strength (IVW β = 0.001, P = 0.929).
- After excluding 210 confounder‑associated SNPs, the null result remains robust (β = −0.007, P = 0.780).
- MR-PRESSO outlier correction does not alter the null finding (β = 0.0013, P = 0.85).
- Multivariable MR and bidirectional MR yield similarly null findings.

## Requirements

- R ≥ 4.2.2
- Packages: `nhanesA`, `TwoSampleMR`, `ieugwasr`, `MVMR`, `ropls`, `impute`, `ggplot2`, `ggpubr`, `broom`, `dplyr`, `tidyr`, `ggrepel`, `survey`, `MRPRESSO`

Install missing packages:
```r
install.packages(c("nhanesA", "ggplot2", "ggpubr", "broom", "dplyr", "tidyr", "ggrepel", "survey"))
# For packages not on CRAN, use:
# remotes::install_github("MRCIEU/TwoSampleMR")
# remotes::install_github("MRCIEU/ieugwasr")
# remotes::install_github("rondolab/MR-PRESSO")
```

## Usage

1. Clone or download the repository.
2. Set the working directory to the repository root.
3. Ensure `MSdata_ST003803_1.txt` is in the root folder.
4. Run the main scripts in order:
   ```r
   source("stage1_metabolomics.R")
   source("stage2_nhanes.R")
   source("stage3_mr.R")
   ```
5. Run the supplementary analyses (requires the main scripts to have been executed first):
   ```r
   source("stage_supplementary_analyses.R")
   ```
Output figures (PDF) and tables (CSV) are saved in the working directory.

## Data availability

| Data | Source |
|------|--------|
| Metabolomics dataset | `MSdata_ST003803_1.txt` (included in this repository; originally from [Metabolomics Workbench ST003803](https://doi.org/10.21228/M8S3803)) |
| NHANES 2011–2018 | [CDC NHANES](https://wwwn.cdc.gov/nchs/nhanes/) |
| GWAS summary statistics | [IEU OpenGWAS](https://gwas.mrcieu.ac.uk/) |
| · Triglycerides | `ieu-b-111` |
| · Hand grip strength | `ukb-b-10215` |
| · Body mass index | `ieu-b-40` |

## License

MIT License. See `LICENSE` for details.

## Citation

If you use this code or data, please cite:

> Chen Y, Tang B, Zhang J. Triglycerides and Skeletal Muscle: Resolving the Body Mass Confound through Relative Muscle Indices and Mendelian Randomization. *Manuscript under review*.

```bibtex
@article{chen2026tgmuscle,
  title   = {Triglycerides and Skeletal Muscle: Resolving the Body Mass Confound through Relative Muscle Indices and Mendelian Randomization},
  author  = {Chen, Yanru and Tang, Bing and Zhang, Jie},
  year    = {2026},
  note    = {Manuscript under review}
}
```
