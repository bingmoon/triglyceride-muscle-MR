# 代码归档溯源（Code Archive Provenance）

本仓库的 `code/*.R` **不是手写存放的散脚本**，而是从项目唯一代码载体 md 中**按代码块抽取**得到的。
本文件记录每个文件的来源与正文指纹，便于审计与回溯。

## 源载体

- 源文件：`4.MSM（ST003803） 1.md`
- 路径（原始项目）：`/Users/bing/Documents/Thinking Library/Paper_Project/4.MSM（ST003803） 1.md`
- 源文件 md5：`b8d02fcfa17cbeeed39539ee4e6c3cf0`
- 🔒 **该 md 未纳入本仓库**：它同时承载内部工作笔记（多轮勘误索引、待办、稿件口径讨论），不适合公开发布。本仓库只取其中**代码块**。

## 代码文件 ← 块映射

| 本仓库文件 | 源 md 块 | 源 md 行范围 | 正文行数 | 正文字节 | 正文 md5（不含本工具添加的 3 行头注释）|
|---|---|---|---|---|---|
| `code/01_stage1_metabolomics.R` | 1 | `L30-280` | 249 | 9531 | `cbd3820621894d235fbef2a3ef7466be` |
| `code/02_stage2_nhanes.R` | 2 | `L283-724` | 440 | 18996 | `a0d5e255a764533d43924eb4e67df30b` |
| `code/03_stage3_mr.R` | 3 | `L726-1279` | 552 | 19852 | `f639b74280b49017cc8482f286fa8a02` |
| `code/04_stage3_mr_ukb_grip.R` | 4 | `L1282-1557` | 274 | 12152 | `4b91075e6aa8b360fbe0be1949c7dcb6` |
| `code/05_supplementary_analyses.R` | 5 | `L1559-1655` | 95 | 3955 | `5397d7e07d29e8b3d2eed6361becd4af` |
| `code/06_stage2_nhanes_weighted_FIXED.R` | 6 | `L1665-1768` | 102 | 5143 | `acbfe5ea8f07d443dde657cd21b296ec` |
| `code/07_power_analysis.R` | 7 | `L1772-1842` | 69 | 3306 | `191d382e680bd67f083eec45e29e4cb7` |
| `code/08_mr_presso_bidirectional.R` | 8 | `L1846-1949` | 102 | 4618 | `28e25760139af53ac59f4ba3a9d28c37` |
| `code/09_abs_grip_bmi_weighted.R` | 9 | `L1967-2146` | 178 | 6110 | `8cf2f6253500a195a01d7511a20a82cc` |
| `code/10_wtsaf2yr_sensitivity.R` | 10 | `L2170-2303` | 132 | 5418 | `045a06985014fa38b61c94d45f183a96` |
| `code/11_wtsaf2yr_sensitivity_trueN.R` | standalone | `/Users/bing/MS/results/WTSAF2YR_sensitivity_trueN.R (verbatim copy)` | 136 | 5593 | `a591abe2ee04cbab3af5b6b04ab0bc6c` |

> 每个 `.R` 文件开头的 3 行 `# ---` 注释由归档工具添加（含来源块与正文 md5），**代码正文逐字未改**。

## 表格来源

- 候选 31 项，按 md5 去重后写入（重复副本只保留一份，来源见下）。

| 仓库内文件名 | 处置 | 详情 |
|---|---|---|
| `Stage1_Table_VIP_Scores.csv` | COPIED | 1296 B |
| `Stage2_Results_Summary.txt` | COPIED | 2277 B |
| `Table_2.8_Adult_DXA_Regression.csv` | COPIED | 945 B |
| `Table_2.9_Adult_Grip_Regression.csv` | COPIED | 675 B |
| `Table_S1_UVMR_Raw_Results.csv` | COPIED | 636 B |
| `Table_S2_Harmonised_SNPs_UVMR.csv` | COPIED | 47365 B |
| `Table_S3_MVMR_Instrumental_Variables.csv` | COPIED | 87756 B |
| `Table_S4_Sensitivity_Analysis.csv` | COPIED | 418 B |
| `Table_S5_Bidirectional_MR.csv` | COPIED | 634 B |
| `Table_S6_Excluded_Confounder_SNPs.csv` | COPIED | 20812 B |
| `Table_S7_UVMR_Clean_Results.csv` | COPIED | 632 B |
| `Table_S8_MVMR_Clean_Results.csv` | COPIED | 399 B |
| `Table_S9_Heterogeneity_Clean.csv` | COPIED | 408 B |
| `Table_S10_Pleiotropy_Clean.csv` | COPIED | 241 B |
| `Table_S11_Top10_SNPs_By_Fstat.csv` | COPIED | 923 B |
| `Table_Power_Analysis_UVMR.csv` | COPIED | 207 B |
| `Table_Power_Curve.csv` | COPIED | 487 B |
| `Table_Weighted_DXA_Regression.csv` | COPIED | 979 B |
| `Table_Weighted_Grip_Regression.csv` | COPIED | 710 B |
| `Supplementary_Tables.xlsx` | COPIED | 106871 B |
| `Table_S12_WTSAF2YR_Sensitivity.csv` | COPIED | 529 B |
| `Table_S4_Sensitivity_Analysis_weighted.csv` | COPIED | 384 B |
| `Table_S4_Sensitivity_Analysis_weighted_full.csv` | COPIED | 950 B |
| `WTSAF2YR_sensitivity.csv` | COPIED | 564 B |
| `weighted_abs_grip_BMI.csv` | COPIED | 794 B |
| `MRPRESSO_forward_main.csv` | COPIED | 278 B |
| `MRPRESSO_reverse_main.csv` | COPIED | 277 B |
| `MRPRESSO_forward_112SNP_main.csv` | COPIED | 277 B |
| `MRPRESSO_forward_112SNP_outliers.csv` | COPIED | 447 B |
| `MS_taskY_steiger_rerun.csv` | COPIED | 510 B |
| `MS_taskY_steiger_heterogeneity.csv` | COPIED | 115 B |

## 图来源（md5 去重）

| 文件名 | 处置 | 详情 |
|---|---|---|
| `Figure_2.6_Adult_RelativeGrip_TG.pdf` | COPIED | 245614 B |
| `Figure_2.8_Adult_RelativeASM_TG.pdf` | COPIED | 274044 B |
| `Figure_2.9_Regression_ForestPlot.pdf` | COPIED | 5263 B |
| `Figure_S1_Absolute_Grip.pdf` | COPIED | 243803 B |
| `Figure_S2_BMI_Paradox_Source.pdf` | COPIED | 472991 B |
| `Stage1_Figure_PCA.pdf` | COPIED | 8249 B |
| `Stage1_Figure_PLSDA_Permutation.pdf` | COPIED | 34715 B |
| `Stage1_Figure_PLSDA_Score.pdf` | COPIED | 8302 B |
| `Stage1_Figure_Univariate_Boxplots.pdf` | COPIED | 97922 B |
| `Stage1_Figure_VIP_Lollipop.pdf` | COPIED | 6226 B |
| `Fig_S3_LeaveOneOut_Raw_Nature.pdf` | COPIED | 14941 B |
| `Fig_S4_Scatter_Clean_Nature.pdf` | COPIED | 7343 B |
| `Fig_S5_Forest_Clean_Nature.pdf` | COPIED | 7700 B |
| `Fig_S6_Funnel_Clean_Nature.pdf` | COPIED | 6556 B |
| `Fig_S7_LeaveOneOut_Clean_Nature.pdf` | COPIED | 7947 B |
| `Figure_2.9_Regression_ForestPlot.pdf` | COPIED | 5284 B |
| `Fig_S3_LeaveOneOut_Raw_Nature.pdf` | dup-of:Fig_S3_LeaveOneOut_Raw_Nature.pdf | /Users/bing/MS/bmj |
| `Fig_S4_Scatter_Clean_Nature.pdf` | dup-of:Fig_S4_Scatter_Clean_Nature.pdf | /Users/bing/MS/bmj |
| `Fig_S5_Forest_Clean_Nature.pdf` | dup-of:Fig_S5_Forest_Clean_Nature.pdf | /Users/bing/MS/bmj |
| `Fig_S6_Funnel_Clean_Nature.pdf` | dup-of:Fig_S6_Funnel_Clean_Nature.pdf | /Users/bing/MS/bmj |
| `Fig_S7_LeaveOneOut_Clean_Nature.pdf` | dup-of:Fig_S7_LeaveOneOut_Clean_Nature.pdf | /Users/bing/MS/bmj |
| `Figure_2.6_Adult_RelativeGrip_TG.pdf` | dup-of:Figure_2.6_Adult_RelativeGrip_TG.pdf | /Users/bing/MS/bmj |
| `Figure_2.8_Adult_RelativeASM_TG.pdf` | dup-of:Figure_2.8_Adult_RelativeASM_TG.pdf | /Users/bing/MS/bmj |
| `Figure_2.9_Regression_ForestPlot.pdf` | dup-of:Figure_2.9_Regression_ForestPlot.pdf | /Users/bing/MS/bmj |
| `Figure_S1_Absolute_Grip.pdf` | dup-of:Figure_S1_Absolute_Grip.pdf | /Users/bing/MS/bmj |
| `Figure_S2_BMI_Paradox_Source.pdf` | dup-of:Figure_S2_BMI_Paradox_Source.pdf | /Users/bing/MS/bmj |
| `Stage1_Figure_PCA.pdf` | dup-of:Stage1_Figure_PCA.pdf | /Users/bing/MS/bmj |
| `Stage1_Figure_PLSDA_Permutation.pdf` | dup-of:Stage1_Figure_PLSDA_Permutation.pdf | /Users/bing/MS/bmj |
| `Stage1_Figure_PLSDA_Score.pdf` | dup-of:Stage1_Figure_PLSDA_Score.pdf | /Users/bing/MS/bmj |
| `Stage1_Figure_Univariate_Boxplots.pdf` | dup-of:Stage1_Figure_Univariate_Boxplots.pdf | /Users/bing/MS/bmj |
| `Stage1_Figure_VIP_Lollipop.pdf` | dup-of:Stage1_Figure_VIP_Lollipop.pdf | /Users/bing/MS/bmj |
| `Supplementary_Figures.pdf` | COPIED | 876521 B |
| `bmjdocument.pdf` | COPIED | 1111076 B |

## 未纳入本仓库的内容（刻意排除）

| 类 | 例子 | 原因 |
|---|---|---|
| 内部工作笔记 | 源 md 的勘误索引 / 待办 / 稿件口径讨论段 | 含未发表稿件的内部不一致记录与投稿策略，**公开发布会自我暴露** |
| 环境快照 | `1.RDate`（20/31 MB）、`.RDataTmp`（14.7 MB） | 体积大、含会话内存、可被脚本重跑 |
| 命令历史 | `.Rhistory`（MS 根 23 KB / M/ 24 KB） | 历史命令明文（凭证泄露面）；已实测无真凭证，但仍排除 |
| 备份链 | `*_backup_20260913_*`（tex / xlsx / csv / bib × 数十个） | 噪音 + 旧口径结论，易被误引 |
| 稿件 | `bmjdocument.tex/.pdf`、`Supplementary_Figures.tex`、`STROBE*.docx`、`coverletter.docx`、`references.bib` | 稿件单独管理（见 README §9），且投稿前不宜公开 |
| 大体积中间对象 | `Table_S4_Sensitivity_SurveyWeighted.rds`(18 MB)、`WTSAF2YR_sensitivity.rds`(19 MB) 等 | 可由脚本重跑，0 信息增量 |
| 编译产物 | `*.aux/.log/.out/.blg`、`.DS_Store` | build artifact |
