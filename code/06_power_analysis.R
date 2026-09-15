# ---
# Archived from: 4.MSM（ST003803） 1.md  (block 10, block 6 deprecated)
# Body md5: 9b6170272fde857b27079450d05fbed1   | 3-line header added by release tooling
# ---
#!/usr/bin/env Rscript
# =====================================================================
# MS_power_analysis_FINAL.R
# 项目 : 甘油三酯(TG) x 骨骼肌（NHANES + 两样本 MR），MS 手稿
# 用途 : MR 功效分析 —— 把"未发现效应"升级为"排除了大于 X 的效应"
# 版本 : R 4.5.2 (macOS)
# ---------------------------------------------------------------------
# 输入（全部为既有磁盘产物；本脚本不重跑任何 MR、不下载任何数据）:
#   M/Table_S1_UVMR_Raw_Results.csv          112-SNP 全集 UVMR (IVW) 汇总统计量
#   M/Table_S7_UVMR_Clean_Results.csv         23-SNP 混杂清洗集 UVMR (IVW) 汇总统计量
#   M/Table_S2_Harmonised_SNPs_UVMR.csv       112-SNP SNP 级 harmonised 数据（算 Cochran Q）
#   M/Table_S9_Heterogeneity_Clean.csv        23-SNP 集 Cochran Q
#   results/MRPRESSO_forward_112SNP_main.csv  MR-PRESSO 离群校正估计（112 集）
# 输出:
#   MS_power_analysis.csv / MS_power_analysis.rds / MS_power_curve.csv
#   （写两份：/Users/bing/MS/results/ 与 会话 results/）
# ---------------------------------------------------------------------
# 方法（功效分析，非新分析）:
#   给定已观测 IVW 估计的 SE，双侧 alpha=0.05、power=80% 时
#     最小可检出效应 MDES = (z_0.975 + z_0.80) * SE = 2.801585 * SE
#   两种"可排除上限"（语义不同，分别报告，禁止混用）:
#     (A) ExclAbove_CIbound = 95% CI 中最远离 0 的一端 = 观测可排除边界
#     (B) MDES_80power      = 80% power 最小可检出效应  = 设计灵敏度边界
#   异质性口径: 随机效应 SE_RE = SE_FE * sqrt(Q / Q_df)
# 参考: Burgess 2011 Int J Epidemiol 40(3):755-764 (doi 10.1093/ije/dyr036)
#       Brion 2013  Int J Epidemiol 42(5):1497-1501 (doi 10.1093/ije/dyt179, PMID 24159078)
# 运行: Rscript MS_power_analysis_FINAL.R
# =====================================================================

M   <- "/Users/bing/MS/M"
RES <- "/Users/bing/MS/results"
SES <- "/Users/bing/MemOmics-Agent/results/memomics-ab3a936b/results"
dir.create(RES, showWarnings = FALSE)
dir.create(SES, showWarnings = FALSE)

z_a <- qnorm(0.975)
z_b <- qnorm(0.80)
K   <- z_a + z_b
cat("=== [0] 常量 ===")
cat(sprintf("z_0.975 = %.6f ; z_0.80 = %.6f ; K = %.6f", z_a, z_b, K))

# ---------------------------------------------------------------- 1
cat("=== [1] 读取既有 UVMR 汇总统计量（不重跑 MR）===")
s1  <- read.csv(file.path(M, "Table_S1_UVMR_Raw_Results.csv"),  stringsAsFactors = FALSE)
s7  <- read.csv(file.path(M, "Table_S7_UVMR_Clean_Results.csv"), stringsAsFactors = FALSE)
ivw1 <- s1[s1$method == "Inverse variance weighted", ]
ivw7 <- s7[s7$method == "Inverse variance weighted", ]
stopifnot(nrow(ivw1) == 1, nrow(ivw7) == 1, ivw1$nsnp == 112, ivw7$nsnp == 23)
cat(sprintf("112-SNP IVW : b = %.10f  SE = %.10f  P = %.6f  nsnp = %d",
            ivw1$b, ivw1$se, ivw1$pval, ivw1$nsnp))
cat(sprintf(" 23-SNP IVW : b = %.10f  SE = %.10f  P = %.6f  nsnp = %d",
            ivw7$b, ivw7$se, ivw7$pval, ivw7$nsnp))

# ---------------------------------------------------------------- 2
# 管道自校验：由 SNP 级 harmonised 数据独立重算 112 集 IVW（b 必须逐位吻合），
# 并判定表内 SE 属固定效应还是随机效应
cat("=== [2] 管道自校验 + 异质性口径判定（112 集）===")
h   <- read.csv(file.path(M, "Table_S2_Harmonised_SNPs_UVMR.csv"), stringsAsFactors = FALSE)
bx  <- h$beta.exposure; by <- h$beta.outcome; sey <- h$se.outcome
ok  <- is.finite(bx) & is.finite(by) & is.finite(sey) & bx != 0
bx <- bx[ok]; by <- by[ok]; sey <- sey[ok]
n_h  <- length(bx)
bj   <- by / bx
sej  <- sey / abs(bx)
w    <- 1 / sej^2
b_fe <- sum(w * bj) / sum(w)
se_fe112 <- sqrt(1 / sum(w))
Q_112 <- sum(w * (bj - b_fe)^2)
df_112 <- n_h - 1
Qp_112 <- pchisq(Q_112, df_112, lower.tail = FALSE)
infl_112 <- sqrt(Q_112 / df_112)
se_re112 <- se_fe112 * infl_112
cat(sprintf("SNP 级 n = %d ; 重算 IVW b = %.10f", n_h, b_fe))
stopifnot(abs(b_fe - ivw1$b) < 1e-12)
cat("PASS: SNP 级重算 b 与 Table_S1 逐位吻合 -> 管道正确")
cat(sprintf("112 集 Cochran Q = %.4f ; Q_df = %d ; Q_p = %.6g ; sqrt(Q/df) = %.6f",
            Q_112, df_112, Qp_112, infl_112))
cat(sprintf("SE_fixed = %.10f ; SE_random = %.10f ; Table_S1 记录 SE = %.10f",
            se_fe112, se_re112, ivw1$se))
stopifnot(abs(se_re112 - ivw1$se) < 1e-12)
cat("PASS: Table_S1 记录的 SE 与随机效应 SE 逐位吻合 -> 稿件报出的 112 集 IVW 为随机效应口径")

# ---------------------------------------------------------------- 3
cat("=== [3] 23 集 Cochran Q（读 Table_S9）===")
s9 <- read.csv(file.path(M, "Table_S9_Heterogeneity_Clean.csv"), stringsAsFactors = FALSE)
q9 <- s9[s9$method == "Inverse variance weighted", ]
stopifnot(nrow(q9) == 1)
Q_23 <- q9$Q; df_23 <- q9$Q_df; Qp_23 <- q9$Q_pval
infl_23 <- sqrt(Q_23 / df_23)
se_fe23 <- ivw7$se / infl_23         # 由随机效应 SE 反解固定效应 SE（派生量）
cat(sprintf("23 集 Cochran Q = %.4f ; Q_df = %d ; Q_p = %.6g ; sqrt(Q/df) = %.6f",
            Q_23, df_23, Qp_23, infl_23))
cat(sprintf("23 集 SE_random(记录) = %.10f ; SE_fixed(派生) = %.10f", ivw7$se, se_fe23))

# ---------------------------------------------------------------- 4
cat("=== [4] MR-PRESSO 离群校正估计（112 集，读磁盘）===")
mp  <- read.csv(file.path(RES, "MRPRESSO_forward_112SNP_main.csv"), stringsAsFactors = FALSE)
mpc <- mp[mp[["MR.Analysis"]] == "Outlier-corrected", ]
stopifnot(nrow(mpc) == 1)
b_mp <- mpc[["Causal.Estimate"]]; se_mp <- mpc[["Sd"]]
cat(sprintf("MR-PRESSO 校正 b = %.10f  SE = %.10f  P = %.6f", b_mp, se_mp, mpc[["P.value"]]))

# ---------------------------------------------------------------- 5
cat("=== [5] MDES / 95% CI / 可排除上限 ===")
rowfun <- function(lab, nsnp, b, se, est) {
  lo <- b - z_a * se; hi <- b + z_a * se
  data.frame(Set = lab, nSNP = nsnp,
             Beta = signif(b, 6), SE = signif(se, 6),
             MDES_80power = signif(K * se, 4),
             CI_low = signif(lo, 4), CI_high = signif(hi, 4),
             ExclAbove_CIbound = signif(max(abs(lo), abs(hi)), 4),
             Estimator = est, stringsAsFactors = FALSE)
}
pw <- rbind(
  rowfun("UVMR Full (112 SNP) - random-effects (as reported)", 112, ivw1$b, ivw1$se,     "random-effects IVW"),
  rowfun("UVMR Full (112 SNP) - fixed-effect (sensitivity)",    112, ivw1$b, se_fe112,    "fixed-effect IVW"),
  rowfun("UVMR Full (112 SNP) - MR-PRESSO outlier-corrected",   112, b_mp,    se_mp,      "MR-PRESSO corrected"),
  rowfun("UVMR Clean (23 SNP) - random-effects (as reported)",   23, ivw7$b, ivw7$se,     "random-effects IVW"),
  rowfun("UVMR Clean (23 SNP) - fixed-effect (derived)",         23, ivw7$b, se_fe23,     "fixed-effect IVW")
)
print(pw, row.names = FALSE)

# ---------------------------------------------------------------- 6
cat("=== [6] 关键结论（供正文引用）===")
for (i in seq_len(nrow(pw))) {
  cat(sprintf("%-52s SE=%.6f | 可排除 |beta| > %.4f SD | 80%% power 可检出 %.4f SD",
              pw$Set[i], pw$SE[i], pw$ExclAbove_CIbound[i], pw$MDES_80power[i]))
}

# ---------------------------------------------------------------- 7
cat("=== [7] Power 曲线 ===")
power_at <- function(tb, se) pnorm(qnorm(0.025) - abs(tb)/se) + pnorm(qnorm(0.975) - abs(tb)/se)
grid <- seq(0, 0.15, by = 0.005)
pc <- data.frame(TrueEffect_SD = grid,
                 Power_Full_112SNP = round(sapply(grid, power_at, se = ivw1$se), 3),
                 Power_Clean_23SNP = round(sapply(grid, power_at, se = ivw7$se), 3))
print(pc[pc$TrueEffect_SD %in% c(0.02, 0.03, 0.04, 0.05, 0.06, 0.08, 0.10), ], row.names = FALSE)

# ---------------------------------------------------------------- 8
cat("=== [8] 落盘 ===")
for (d in c(RES, SES)) {
  write.csv(pw, file.path(d, "MS_power_analysis.csv"), row.names = FALSE)
  write.csv(pc, file.path(d, "MS_power_curve.csv"), row.names = FALSE)
  saveRDS(list(
      table = pw, power_curve = pc,
      heterogeneity = list(Q_112 = Q_112, df_112 = df_112, Q_p_112 = Qp_112,
                           Q_23  = Q_23,  df_23  = df_23,  Q_p_23  = Qp_23,
                           inflation_112 = infl_112, inflation_23 = infl_23),
      reported_SE_estimator = "random-effects IVW (verified for 112-SNP set to 1e-12)",
      source = list(S1 = file.path(M, "Table_S1_UVMR_Raw_Results.csv"),
                    S7 = file.path(M, "Table_S7_UVMR_Clean_Results.csv"),
                    S2 = file.path(M, "Table_S2_Harmonised_SNPs_UVMR.csv"),
                    S9 = file.path(M, "Table_S9_Heterogeneity_Clean.csv"),
                    MRPRESSO = file.path(RES, "MRPRESSO_forward_112SNP_main.csv")),
      method = "MDES = (qnorm(0.975)+qnorm(0.80))*SE ; ExclAbove = max(|CI_lo|,|CI_hi|) ; SE_RE = SE_FE*sqrt(Q/Q_df)",
      constant_K = K),
    file.path(d, "MS_power_analysis.rds"))
  cat(sprintf("written -> %s/MS_power_analysis.{csv,rds} + MS_power_curve.csv", d))
}
cat("DONE")
