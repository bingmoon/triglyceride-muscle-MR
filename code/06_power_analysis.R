# ---
# Archived from: 4.MSM（ST003803） 1.md  (L1772-1842, block 7)
# Body md5: 191d382e680bd67f083eec45e29e4cb7   | 3-line header added by release tooling
# ---
# =====================================================================
# MSM 项目 - 补充分析 A：MR 功效分析（把"未发现"升级为"可排除"）
# 【新增块 — 应对审稿人对阴性 MR 的必问质疑："是否只是 power 不足？"】
# 原理：
#   给定已观测的 IVW 标准误 SE，双侧 α=0.05、power=80% 时
#   最小可检出效应 |b|_min = (z_{0.975} + z_{0.80}) × SE = 2.8016 × SE
#   "可排除上限" = 观测估计 95% CI 中最远离 0 的一端
# 结果（2026-09-13 实测，结局 ukb-b-10215 握力，N=461,026）：
#   112-SNP 全集 : 可排除 > 0.0204 SD 的 TG 因果效应（80% power 可检出 0.0278）
#   23-SNP 纯净集: 可排除 > 0.0569 SD 的 TG 因果效应（80% power 可检出 0.0712）
# =====================================================================
setwd("/Users/bing/MS")

out_dir <- "/Users/bing/MS/M"
dir.create(out_dir, showWarnings = FALSE)

cat("\n========== 1. 读取已有 UVMR 结果 ==========\n")
s1 <- read.csv(file.path(out_dir, "Table_S1_UVMR_Raw_Results.csv"))
s7 <- read.csv(file.path(out_dir, "Table_S7_UVMR_Clean_Results.csv"))
ivw_r <- s1[s1$method == "Inverse variance weighted", ]
ivw_c <- s7[s7$method == "Inverse variance weighted", ]

cat(sprintf("全集  IVW: b=%.4f, SE=%.4f, P=%.4f (nsnp=%d)\n",
            ivw_r$b, ivw_r$se, ivw_r$pval, ivw_r$nsnp))
cat(sprintf("纯净集 IVW: b=%.4f, SE=%.4f, P=%.4f (nsnp=%d)\n",
            ivw_c$b, ivw_c$se, ivw_c$pval, ivw_c$nsnp))

cat("\n========== 2. 计算最小可检出效应与可排除上限 ==========\n")
z_a <- qnorm(0.975)   # 1.959964
z_b <- qnorm(0.80)    # 0.8416212
calc_power <- function(b, se, nsnp, lab) {
  min_detect <- (z_a + z_b) * se                 # 80% power 下的最小可检出 |b|
  ci_lo <- b - z_a * se; ci_hi <- b + z_a * se   # 95% CI
  data.frame(Set = lab, nSNP = nsnp,
             Beta = round(b, 4), SE = round(se, 4),
             MinDetectable_80pct = round(min_detect, 4),
             CI_low = round(ci_lo, 4), CI_high = round(ci_hi, 4),
             ExcludableAbove = round(max(abs(ci_lo), abs(ci_hi)), 4),
             row.names = NULL)
}
pw <- rbind(
  calc_power(ivw_r$b, ivw_r$se, ivw_r$nsnp, "UVMR Full (112 SNP)"),
  calc_power(ivw_c$b, ivw_c$se, ivw_c$nsnp, "UVMR Clean (23 SNP)")
)
print(pw)

cat("\n========== 3. Power 曲线（不同真实效应量下） ==========\n")
power_at <- function(true_b, se) {
  pnorm(qnorm(0.025) - abs(true_b) / se) + pnorm(qnorm(0.975) - abs(true_b) / se)
}
grid <- seq(0, 0.15, by = 0.005)
pw_curve <- data.frame(
  TrueEffect_SD = grid,
  Power_Full_112SNP = round(sapply(grid, power_at, se = ivw_r$se), 3),
  Power_Clean_23SNP = round(sapply(grid, power_at, se = ivw_c$se), 3)
)
print(pw_curve[pw_curve$TrueEffect_SD %in% c(0.02, 0.03, 0.04, 0.05, 0.06, 0.08, 0.10), ])

write.csv(pw, file.path(out_dir, "Table_Power_Analysis_UVMR.csv"), row.names = FALSE)
write.csv(pw_curve, file.path(out_dir, "Table_Power_Curve.csv"), row.names = FALSE)

cat("\n========== 4. 论文可用表述（供 Results/Limitations 直接引用） ==========\n")
cat(sprintf("With 112 instruments (N = 461,026), our analysis had 80%% power to detect
a causal effect of TG on grip strength of %.4f SD or larger (per 1-SD TG).
The observed estimate (beta = %.4f, 95%% CI %.4f to %.4f) excludes causal
effects larger than %.4f SD.\n",
            pw$MinDetectable_80pct[1], pw$Beta[1], pw$CI_low[1], pw$CI_high[1], pw$ExcludableAbove[1]))

cat("\n✅ 功效分析完成，Table_Power_Analysis_UVMR.csv + Table_Power_Curve.csv 已落盘。\n")
