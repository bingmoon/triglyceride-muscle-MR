# triglyceride-muscle-MR

**循环甘油三酯（TG）与骨骼肌质量/肌力的因果关系研究 —— 两阶段分析代码归档（NHANES 2011–2018 / 两样本 MR）。**

> 本仓库是论文方法学的**可复现代码归档**。所有结果表与图均由 `code/` 下的脚本生成，可按 §4 顺序重跑得到。

---

## 1. 项目简介

本研究分两阶段检验「TG 升高 → 骨骼肌质量与肌力下降」这一假设：阶段一用 NHANES 2011–2018 四周期（G/H/I/J）复杂抽样加权回归做人群层面验证；阶段二用两样本与多变量孟德尔随机化（MR/MVMR）做因果推断与反向方向检验。

---

## 2. 环境要求

| 项 | 值 |
|---|---|
| 语言 | R ≥ 4.2（本项目开发于 **R 4.5.2**，macOS） |
| 网络 | 阶段一（`nhanesA` 下载 NHANES）、阶段二（IEU OpenGWAS，**需凭证**）均需联网 |

### 2.1 依赖包（**按安装来源分列 —— 只装 CRAN 会失败**）

```r
# ---- CRAN ----
install.packages(c(
  "dplyr", "tidyr", "tibble", "stringr", "jsonlite",
  "ggplot2", "ggpubr", "ggrepel", "ggsci",
  "broom", "survey", "nhanesA", "openxlsx",
  "MVMR", "BiocManager"
))

# ---- GitHub-only（不写来源读者必装不上）----
remotes::install_github("MRCIEU/TwoSampleMR")   # 两样本 MR / harmonise / 异质性
remotes::install_github("MRCIEU/ieugwasr")      # OpenGWAS 客户端 + 混杂筛查
remotes::install_github("rondolab/MR-PRESSO")   # MR-PRESSO 水平多效性
```

> 阶段一同时使用了 NHANES 官方 `survey` 设计与 `broom::tidy()` 取模型系数。

---

## 3. 数据来源

| 阶段 | 数据层 | 来源与获取方式 | 是否随仓库提供 |
|---|---|---|---|
| 一 | NHANES 2011–2018（周期 G/H/I/J）：DEMO / BMX / TRIGLY / FASTQX / MGX / DXX | R 包 `nhanesA::nhanes()` **在线下载**（NHANES 官方公开数据） | ⛔ 不随附（脚本自动下载）；派生中间对象见 `data/derived/` |
| 二 | GWAS 汇总统计：TG = `ieu-b-111`、BMI = `ieu-b-40`、握力 = `ukb-b-10215` | **IEU OpenGWAS**（需 token，见 §5） | ⛔ 不随附（需 token 下载） |
| 二 | 混杂筛查命中列表（PhenoScanner / ieugwasr） | 由 `ieugwasr` 查询 OpenGWAS 生成 | ✅ `data/phenoscanner_sig_confounders_filtered.csv` |

**随附的派生数据**（`data/derived/`，小体积、可由脚本重跑，随附以便免下载核对）：

| 文件 | 内容 |
|---|---|
| `NHANES_Adult_Grip_Main.rds` | 阶段一握力队列（G+H 周期，Age ≥ 20）已清洗数据框 |
| `NHANES_Adult_DXA_All.rds` | 阶段一 DXA 队列（G/H/I/J 四周期，Age ≥ 20）已清洗数据框 |

---

## 4. 运行顺序

```text
code/01_stage1_nhanes.R                NHANES 主队列：拉数 → 握力/DXA 队列 → 回归 → 图 → 汇总
      │                                （⛔ 该脚本的回归为**未加权 lm**；加权正确路线见 05）
code/02_stage2_mr.R                    阶段二：两样本 MR + MVMR（旧结局 ieu-b-39）
code/03_stage2_mr_ukb_grip.R           阶段二（主用）：结局换 ukb-b-10215，从零生成图表
      │
      ├─ code/04_supplementary_analyses.R  补充分析合集：加权回归 / MR-PRESSO / 条件 F
      ├─ code/05_stage1_nhanes_weighted_FIXED.R  ⭐阶段一**加权修正版**（svydesign + svyglm，取代 01 的未加权）
      ├─ code/06_power_analysis.R           MR 功效分析（把「未发现」升级为「可排除」）
      ├─ code/07_mr_presso_bidirectional.R  MR-PRESSO 双向（NbDistribution=10000）
      ├─ code/08_abs_grip_bmi_weighted.R    绝对握力+BMI 加权回归（带诊断）
      └─ code/09_wtsaf2yr_sensitivity.R     空腹子样本权重（WTSAF2YR）敏感性分析
         code/10_wtsaf2yr_sensitivity_trueN.R 同上，N 用 nobs(fit) 的**权威版**
```

> 编号已统一为两阶段体系：`01`–`10` 连续；文件名中的 `stage1` = NHANES 人群层、`stage2` = MR 因果层。

### ⚠️ 4.1 运行前**必须**修改的硬编码路径

所有脚本都以绝对路径写死工作目录，**共 19 处、分布于 9 个脚本**（`01/02/03/05/06/07/08/09/10`；`04` 无绝对路径）。其中 `01_stage1_nhanes.R` 同时出现 `/Users/bing/MSM`（L11）与 `/Users/bing/MS`（L349）。请统一改成你本地的仓库根：

```bash
# macOS / Linux
sed -i '' 's#/Users/bing/MSM*#/你的/仓库/路径#g' code/*.R      # macOS
sed -i    's#/Users/bing/MSM*#/你的/仓库/路径#g'  code/*.R      # GNU sed
```

同时确保仓库根下存在脚本期望的输出目录：`results/`（脚本多数带 `dir.create("results", showWarnings = FALSE)`）。

> 🚧 **known limitation**：上述硬编码路径**尚未参数化**（未改为 `here::here()` / 相对路径）。在其修复前，本仓库的定位是**可审计的代码归档**（每行代码可追溯，见 §3 与 `docs/code_archive_provenance.*`），而非开箱即跑的复现包。

### 4.2 未固定 / 需注意的随机性与版本漂移

| 项 | 说明 |
|---|---|
| `set.seed(2026)` | `03`、`05`、`07` 已固定；**置换检验与 MR-PRESSO 的末位数字仍可能随 R/包版本抖动**（`07` 的 MR-PRESSO 用 `NbDistribution = 10000`，重复运行 outlier 集合可能略有差异） |
| NHANES 下载 | 脚本每次运行都会从 NHANES 官网**重新下载**当前存档文件；NHANES 会发布修订版，样本量可能在极小量级上漂移 |
| OpenGWAS | GWAS 汇总统计为**版本化数据集**（`ieu-b-*` / `ukb-b-*` 固定 ID），但服务端偶有元数据修订 |

### 4.3 ⛔ 不可混用的口径（**读结果前必读**）

| # | 口径 | 正确用法 |
|---|---|---|
| 1 | **加权 vs 未加权**（阶段一） | 正式结论一律用 **加权**（`svyglm`，`05`/`08`/`09`/`10`）。`01` 的 `lm` 版本为早期版本，其绝对握力+BMI 估计（β≈−0.108，P=0.6995）**已作废**，加权真值为 **β=0.3679, P=0.400** |
| 2 | **DXA 队列 2 周期 vs 4 周期** | 权威口径 = **四周期 G/H/I/J，n=4,644**（`01` 早期仅取 G/H → 掉到 2,508；`05` 已补齐权重） |
| 3 | **正向 MR 的结局 GWAS** | 主用 **`ukb-b-10215`**（握力）。`02` 使用的旧结局 **`ieu-b-39` 已废弃**，不要引用其数字 |
| 4 | **正向 MR 工具集** | 主分析 **112 SNP**（`Table_S1`）；清洗后 **23 SNP**（`Table_S7`）；MVMR 联合集 **438 SNP**（`Table_S3`）。三者**不是同一集合**，勿混算 |
| 5 | **MR-PRESSO 工具集** | `07` 实跑为**全工具集（280 SNP）**；与主分析 112 SNP 不同集，稿件中已声明为敏感性分析 |
| 6 | **功效分析的 SE 口径** | `06_power_analysis.R` 读的是 `Table_S1`/`Table_S7` 里**已报出的 IVW SE**——经 SNP 级重算验证，该 SE 是**随机效应**口径 （`SE_RE = SE_FE x sqrt(Q/Q_df)`，112 集在 1e-12 内吻合：SE_FE=0.0043576, Q=577.2641, df=111）。固定效应 SE 小 **2.28 倍**；MR-PRESSO 离群校正 SE=0.0070615 属**第三口径**。三口径的「可排除上限」**不可混列**（0.0094 / 0.0163 / 0.0204 SD）。 |

---

## 5. 凭证设置（IEU OpenGWAS）

阶段二需要 OpenGWAS token。脚本中**显式**读取环境变量 `OPENGWAS_JWT` 的是 `code/07_mr_presso_bidirectional.R`；`code/02`、`03` 通过 `ieugwasr` / `TwoSampleMR` 使用同一环境变量。

```r
# 步骤 1：注册并生成 token
#   访问 https://api.opengwas.io/ 注册 → 在 Account 页面生成 JWT
# 步骤 2：写入 ~/.Renviron（R 启动时自动读取）
usethis::edit_r_environ()
#   追加一行（⚠️ 不要用真实 JWT 前缀写示例，否则会被密钥扫描器/CI 拦下）：
#   OPENGWAS_JWT=<在此粘贴你的 JWT token>
# 步骤 3：保存后**重启 R 会话**才生效
# 步骤 4：校验
Sys.getenv("OPENGWAS_JWT")        # 非空
ieugwasr::user()                  # 能返回账号信息即有效
```

> **凭证卫生**：`.Renviron` / `.Rhistory` / `*.RDate` 已在 `.gitignore` 中。**请勿**把 token 粘贴进任何 `.R` 脚本或 notebook —— 本归档的脚本一律只读环境变量。

---

## 6. 预期产出

Analysis outputs (tables and figures) are available in the manuscript's supplementary materials and are not included in this repository.

### 6.1 中间对象

脚本会在工作目录生成 `*.rds`（模型对象、队列快照）。本仓库**只随附 2 个体积小的派生物**（见 §3），其余（10 MB 级队列快照）可由脚本重跑得到，未纳入。

---

## 7. 复现验证

**三步核对法**：① 按 §4 顺序跑脚本 → ② 用下表比对关键数值 → ③ 数值一致即复现成功（图的形态差异可接受，数值须一致）。

### 7.1 复现锚点表

| # | 阶段 | 指标 | 期望值 | 出处（补充材料） |
|---|---|---|---|---|
| 1 | 一 | 握力队列样本量 | **n = 4,180**（G 2,017 / H 2,163） | `Stage2_Results_Summary.txt`（⚠️ 该文件名为**旧三阶段命名**遗留：旧 Stage 2 = NHANES 层，现对应 `stage1`；补充材料文件名未改，以遵守「不改补充材料」约束） |
| 2 | 一 | 相对握力 ~ Log2_TG（未加权 lm） | β = −0.05064, P = 5.168e−33, R² = 0.3899 | `Stage2_Results_Summary.txt` |
| 3 | 一 | DXA 队列样本量 | **n = 4,644**（四周期） | `Stage2_Results_Summary.txt` |
| 4 | 一 | 相对 ASM ~ Log2_TG（未加权 lm） | β = −0.01177, P = 8.328e−106, R² = 0.6331 | `Stage2_Results_Summary.txt` |
| 5 | 一 | **加权**相对握力（MEC 权重） | β = −0.0598, P < 0.001, n = 4,184 | `Table_S12` |
| 6 | 一 | **加权**相对 ASM（MEC 权重） | β = −0.0119, P < 0.001, n = 2,508 | `Table_S12` |
| 7 | 一 | **加权**绝对握力 + BMI | β = **0.3679**, SE = 0.4300, P = **0.400**, n = 4,180 | `Table_S4` / `Table_S12` |
| 8 | 二 | 正向 UVMR（112 SNP）IVW | β = 0.00088, P = 0.929 | `Table_S1` |
| 9 | 二 | 正向清洗集（23 SNP）IVW | β = −0.00711, P = 0.780 | `Table_S7` |
| 10 | 二 | 反向 MR（167 SNP）IVW | β = −0.07550, P = 0.0519 | `Table_S5` |
| 11 | 二 | 正向异质性（清洗集 IVW） | Q = 51.19, Q_df = 22, **P = 4.04e−4**（存在显著异质性 → 主分析用随机效应 IVW） | `Table_S9` |
| 12 | 二 | 反向异质性 | Q = 1109.00, Q_df = 166, P = 3.76e−139 | 未随附（本仓库复算产物） |
| 13 | 二 | MR-PRESSO 正向（112 SNP） | 原始 β = 0.00088, P = 0.9296；离群校正后 β = 0.00244, P = 0.7307 | `MRPRESSO_forward_112SNP_main.csv` |
| 14 | 二 | MR-Egger 截距（清洗集） | intercept = 0.00087, P = 0.480（无水平多效性证据） | `Table_S10` |
| 15 | 二 | MR 功效（可排除的效应下限） | 112 SNP → 可排除 \|β\| > 0.0204 SD；23 SNP → > 0.0569 SD | `Table_Power_Analysis_UVMR.csv` |

### 7.2 **负向 / 阴性结果也列出**（审稿人最常复核的部分）

| 指标 | 期望值（阴性） | 出处（补充材料） |
|---|---|---|
| 正向 MR（112 SNP）IVW | P = 0.929 —— **未发现 TG → 握力的因果证据** | `Table_S1` |
| 正向 MR 清洗集（23 SNP） | P = 0.780 —— 仍为阴性 | `Table_S7` |
| 反向 MR（167 SNP） | P = 0.0519 —— **未达显著**（边界） | `Table_S5` |
| MVMR 清洗后（8 SNP） | P = 0.345 —— 阴性 | `Table_S8` |
| 绝对握力 + BMI（加权） | P = 0.400 —— 加 BMI 后不再显著 | `Table_S4` |

**已知的不可完全复现项**（诚实披露）：

1. **`07` MR-PRESSO** 的离群 SNP 集合对 `NbDistribution` 敏感，重复运行可能有 ±1 个 SNP 的抖动。
2. **硬编码绝对路径（known limitation）**：脚本以原作者机器路径写死工作目录（共 19 处、分布于 9 个脚本，详见 §4.1）。**本仓库不保证「开箱即跑（out-of-the-box reproduction）」**——运行前必须先按 §4.1 替换路径，否则脚本会在拉数/写盘阶段失败。
3. **补充材料文件名沿用旧三阶段命名**：如 `Stage2_Results_Summary.txt`（旧 Stage 2 = 现 `stage1` NHANES 层）。文件内容与正文两阶段体系一致，仅文件名未改（遵守「不改补充材料」约束）。
4. **Table 1 的两项描述统计无法由随附派生数据复现**（2026-09 复核）：握力队列 BMI SD 稿件报 **6.8**，`data/derived/NHANES_Adult_Grip_Main.rds` 复算为 **6.97**（未加权）/ **6.92**（加权）；DXA 队列 TG IQR 稿件报 **63–138**，`NHANES_Adult_DXA_All.rds` 复算为 **62–139**（未加权）/ **63–139**（加权）。两项在未加权与加权两种口径下**均不可复现**，且随附 rds 与作者本地镜像（`M/1.RDate`）逐字节一致（`identical() = TRUE`），故**非 NHANES 数据集版本漂移所致**。其余 Table 1 数字（N、年龄、女性 n%、握力 TG 98 (68–144)、DXA BMI 28.7 (6.9)）均可逐字复现。稿件按作者裁定**保留原值**，并在 Table 1 表注**具体披露**复算值与不可复现事实（“…were not reproducible under unweighted, MEC-weighted, or alternative quantile conventions…”）；补充材料 S1–S12 与图件均未改动。

---

## 8. 引用方式

```bibtex
@misc{tang2026triglyceride,
  title  = {Triglyceride and skeletal muscle mass/strength: a two-stage
            analysis (NHANES 2011--2018, two-sample MR)},
  author = {Tang, Bing},
  year   = {2026},
  note   = {Code archive},
  howpublished = {\url{https://github.com/bingmoon/triglyceride-muscle-MR}}
}
```

**方法学依赖（请一并引用其原始文献）**：

| 软件 | 引用 |
|---|---|
| TwoSampleMR | Hemani G, et al. *eLife* 2018;7:e34408 |
| MR-PRESSO | Verbanck M, et al. *Nat Genet* 2018;50:693–698 |
| MVMR | Sanderson E, et al. *Am J Epidemiol* 2019;188:1568–1577 |
| nhanesA / NHANES | CDC/NCHS National Health and Nutrition Examination Survey |

---

## 9. 联系方式

- **通讯作者**：Bing Tang
- **问题反馈**：请开 [GitHub Issues](https://github.com/bingmoon/triglyceride-muscle-MR/issues)
- **稿件状态**：分析对应的论文稿件单独管理，不在本仓库内。

---

## 附：仓库结构

```text
triglyceride-muscle-MR/
├── README.md
├── LICENSE                 (MIT)
├── .gitignore
├── code/                   (10 个 R 脚本；01→10 见 §4，两阶段连续编号)
├── data/
│   ├── phenoscanner_sig_confounders_filtered.csv
│   └── derived/            (2 个小体积中间对象)
└── docs/
    ├── code_archive_provenance.{md,json}   (每个代码文件的来源 md 行号 + 正文 md5)
    └── release_audit_card.png             (发布包审计卡)
```
