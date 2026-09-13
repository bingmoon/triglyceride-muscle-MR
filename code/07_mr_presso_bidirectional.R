# ---
# Archived from: 4.MSM（ST003803） 1.md  (L1846-1949, block 8)
# Body md5: 28e25760139af53ac59f4ba3a9d28c37   | 3-line header added by release tooling
# ---
# =====================================================================
# 【块 7·新增 2026-09-13】MR-PRESSO 双向分析（真值版）
# ---------------------------------------------------------------------
# 目的：补齐 MR-PRESSO 溯源缺口 —— 此前稿件引用的 MR-PRESSO 数字在磁盘上无任何产物
# 运行：本地 R，NbDistribution = 10000，SignifThreshold = 0.05，set.seed(2026)
# 产出：/Users/bing/MS/results/
#         MRPRESSO_both_console.txt      （完整 console）
#         MRPRESSO_reverse.rds / MRPRESSO_reverse_main.csv
#         MRPRESSO_forward.rds / MRPRESSO_forward_main.csv
#
# 【口径声明 —— 重要】
#   反向（握力→TG）：ukb-b-10215 → ieu-b-111，harmonise 后 167 SNP  ← 与稿件 L273 一致 ✅
#   正向（TG→握力）：ieu-b-111 → ukb-b-10215，harmonise 后 280 SNP  ← 稿件 L279/L343 原称
#                    "112 instruments"，**280 与 112 是两个不同工具集**，已在正文如实改写
# =====================================================================

# ---------- 0. 环境初始化 ----------
suppressMessages({
    library(TwoSampleMR)
    library(ieugwasr)
    library(MRPRESSO)
})

setwd("/Users/bing/MS")
dir.create("results", showWarnings = FALSE)

# 检查 token（token 存放于 ~/.Renviron 的 OPENGWAS_JWT，不写入本文件）
user_info <- tryCatch(ieugwasr::user(), error = function(e) NULL)
if (is.null(user_info)) stop("Token 无效，请检查 ~/.Renviron 中的 OPENGWAS_JWT")

set.seed(2026)
NB_DIST    <- 10000   # 10000 才满足 outlier test 的精度要求（3000 会报 unstable）
SIG_THRESH <- 0.05

sink("results/MRPRESSO_both_console.txt", split = TRUE)

# =====================================================================
# 第一部分：反向 MR-PRESSO（握力 → TG，167 SNP）
# =====================================================================
exp_grip <- extract_instruments("ukb-b-10215", p1 = 5e-8)          # 176 IVs
out_tg   <- extract_outcome_data(snps = exp_grip$SNP,
                                 outcomes = "ieu-b-111", proxies = TRUE)
dat_rev  <- harmonise_data(exposure_dat = exp_grip, outcome_dat = out_tg)
dat_rev  <- dat_rev[dat_rev$mr_keep, ]                              # 167 SNP

presso_rev <- mr_presso(
    BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure",
    SdOutcome = "se.outcome",     SdExposure = "se.exposure",
    data = as.data.frame(dat_rev),
    NbDistribution  = NB_DIST,
    SignifThreshold = SIG_THRESH,
    OUTLIERtest = TRUE, DISTORTIONtest = TRUE
)
print(presso_rev$`Main MR results`)
print(presso_rev$`MR-PRESSO results`$`Global Test`)
print(presso_rev$`MR-PRESSO results`$`Distortion Test`)

saveRDS(presso_rev, "results/MRPRESSO_reverse.rds")
write.csv(presso_rev$`Main MR results`,
          "results/MRPRESSO_reverse_main.csv", row.names = FALSE)

# ---- 实测真值（反向）----
#   Raw:               β = -0.0755, SE = 0.03884, P = 0.0536
#   Outlier-corrected: β = -0.0541, SE = 0.02324, P = 0.0212   ← 校正后名义显著
#   Global Test:       RSSobs = 1124.287, P < 1e-04
#   Outliers:          20 个（indices 9,14,35,41,42,58,61,65,71,72,86,95,108,112,129,133,141,142,152,160）
#   Distortion Test:   coef = -39.443, P = 0.134

# =====================================================================
# 第二部分：正向 MR-PRESSO（TG → 握力，280 SNP——全工具集）
# =====================================================================
exp_tg  <- extract_instruments("ieu-b-111", p1 = 5e-8)              # 313 IVs
exp_tg  <- exp_tg[!duplicated(exp_tg$SNP), ]                        # 去重
out_grip <- extract_outcome_data(snps = exp_tg$SNP,
                                 outcomes = "ukb-b-10215", proxies = TRUE)
dat_fwd <- harmonise_data(exposure_dat = exp_tg, outcome_dat = out_grip)
dat_fwd <- dat_fwd[dat_fwd$mr_keep, ]                               # 280 SNP

presso_fwd <- mr_presso(
    BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure",
    SdOutcome = "se.outcome",     SdExposure = "se.exposure",
    data = as.data.frame(dat_fwd),
    NbDistribution  = NB_DIST,
    SignifThreshold = SIG_THRESH,
    OUTLIERtest = TRUE, DISTORTIONtest = TRUE
)
print(presso_fwd$`Main MR results`)
print(presso_fwd$`MR-PRESSO results`$`Global Test`)
print(presso_fwd$`MR-PRESSO results`$`Distortion Test`)

saveRDS(presso_fwd, "results/MRPRESSO_forward.rds")
write.csv(presso_fwd$`Main MR results`,
          "results/MRPRESSO_forward_main.csv", row.names = FALSE)

# ---- 实测真值（正向）----
#   Raw:               β = -0.00902, SE = 0.007305, P = 0.2180
#   Outlier-corrected: β = -0.00718, SE = 0.005755, P = 0.2132
#   Global Test:       RSSobs = 1237.45, P < 1e-04
#   Outliers:          16 个（indices 40,63,98,114,124,130,136,137,144,149,159,161,224,232,264,275）
#   Distortion Test:   coef = -25.594, P = 0.6354

sink()
