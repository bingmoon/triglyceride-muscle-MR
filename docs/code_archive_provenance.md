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
| `code/01_stage1_nhanes.R` | 2 | `L283-724` | 440 | 18996 | `a0d5e255a764533d43924eb4e67df30b` |
| `code/02_stage2_mr.R` | 3 | `L726-1279` | 552 | 19852 | `f639b74280b49017cc8482f286fa8a02` |
| `code/03_stage2_mr_ukb_grip.R` | 4 | `L1282-1557` | 274 | 12152 | `4b91075e6aa8b360fbe0be1949c7dcb6` |
| `code/04_supplementary_analyses.R` | 5 | `L1559-1655` | 95 | 3955 | `5397d7e07d29e8b3d2eed6361becd4af` |
| `code/05_stage1_nhanes_weighted_FIXED.R` | 6 | `L1665-1768` | 102 | 5143 | `acbfe5ea8f07d443dde657cd21b296ec` |
| `code/06_power_analysis.R` | 7 | `L1772-1842` | 69 | 3306 | `191d382e680bd67f083eec45e29e4cb7` |
| `code/07_mr_presso_bidirectional.R` | 8 | `L1846-1949` | 102 | 4618 | `28e25760139af53ac59f4ba3a9d28c37` |
| `code/08_abs_grip_bmi_weighted.R` | 9 | `L1967-2146` | 178 | 6110 | `8cf2f6253500a195a01d7511a20a82cc` |
| `code/09_wtsaf2yr_sensitivity.R` | 10 | `L2170-2303` | 132 | 5418 | `045a06985014fa38b61c94d45f183a96` |
| `code/10_wtsaf2yr_sensitivity_trueN.R` | standalone | `/Users/bing/MS/results/WTSAF2YR_sensitivity_trueN.R (verbatim copy)` | 136 | 5593 | `a591abe2ee04cbab3af5b6b04ab0bc6c` |

> 源 md **块 1**（前临床代谢组学阶段）已随「三阶段 → 两阶段」重构移出本仓库；**块号保持原编号**，故现存块号为 2–10（连续但非从 1 起）。
>
> 每个 `.R` 文件开头的 3 行 `# ---` 注释由归档工具添加（含来源块与正文 md5），**代码正文逐字未改**。

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
