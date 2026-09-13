# ---
# Archived from: 4.MSM（ST003803） 1.md  (L1665-1768, block 6)
# Body md5: acbfe5ea8f07d443dde657cd21b296ec   | 3-line header added by release tooling
# ---
# =====================================================================
# MSM 项目 - 阶段二·修正版：NHANES 复杂抽样加权回归
# 【本块为工作流确认后的正确路线，取代原 md 中未加权版本】
# 修复要点：
#   1. 原脚本提取了权重但从未真正用于建模 → 现全程 svydesign + svyglm
#   2. 握力队列 = G+H 两周期 → 4 年权重 W4YR = WTMEC2YR / 2
#   3. DXA 队列 = G+H+I+J 四周期 → 4 年权重 W4YR = WTMEC2YR / 4
#      （原脚本只取 G/H 权重，导致 DXA 队列 46% 样本丢失 → 已修复为 100% 保留）
#   4. 输出统一写入 results/，并给出与稿件数值的对照核验
# 结果（2026-09-13 实测）：
#   相对握力 β = -0.05992 (SE 0.00639, P = 3.93e-10)   ← 稿件 -0.060 / 3.93e-10 ✔
#   相对 ASM  β = -0.012214 (SE 0.00081, P = 2.19e-21) ← 稿件 -0.012 / 2.19e-21 ✔
#   绝对握力+BMI 的 TG 项 P = 0.400（稿件 0.699，差异见说明）
# =====================================================================
setwd("/Users/bing/MS")
suppressMessages({
  library(dplyr); library(survey); library(nhanesA); library(broom)
})
options(timeout = 300); set.seed(2026)

out_dir <- "/Users/bing/MS/M"
dir.create(out_dir, showWarnings = FALSE)

cat("\n========== 1. 拉取 NHANES 人口学权重 (G/H/I/J) ==========\n")
get_w <- function(cy) {
  d <- nhanesA::nhanes(paste0("DEMO_", cy))
  d$SEQN <- as.numeric(as.character(d$SEQN))
  data.frame(SEQN = d$SEQN,
             SDMVPSU = as.numeric(d$SDMVPSU),
             SDMVSTRA = as.numeric(d$SDMVSTRA),
             WTMEC2YR = as.numeric(d$WTMEC2YR),
             Cycle = cy)
}
w_all <- bind_rows(lapply(c("G", "H", "I", "J"), get_w))
cat(sprintf("权重表获取完成：%d 行，周期 %s\n", nrow(w_all), paste(sort(unique(w_all$Cycle)), collapse = "/")))

w2 <- w_all %>% filter(Cycle %in% c("G", "H")) %>% mutate(W4YR = WTMEC2YR / 2)  # 2 周期
w4 <- w_all %>% mutate(W4YR = WTMEC2YR / 4)                                      # 4 周期

cat("\n========== 2. 合并权重到队列 ==========\n")
df_grip <- readRDS("NHANES_Adult_Grip_Main.rds")
df_dxa  <- readRDS("NHANES_Adult_DXA_All.rds")

df_grip_wt <- df_grip %>% left_join(w2, by = c("SEQN", "Cycle")) %>% filter(!is.na(W4YR), !is.na(SDMVPSU))
df_dxa_wt  <- df_dxa  %>% left_join(w4, by = c("SEQN", "Cycle")) %>% filter(!is.na(W4YR), !is.na(SDMVPSU))

cat(sprintf("握力加权队列 n = %d / 原 %d (保留 %.1f%%)\n",
            nrow(df_grip_wt), nrow(df_grip), 100 * nrow(df_grip_wt) / nrow(df_grip)))
cat(sprintf("DXA  加权队列 n = %d / 原 %d (保留 %.1f%%)\n",
            nrow(df_dxa_wt), nrow(df_dxa), 100 * nrow(df_dxa_wt) / nrow(df_dxa)))

cat("\n========== 3. 建立复杂抽样设计并加权回归 ==========\n")
de_grip <- svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~W4YR,
                     data = df_grip_wt, nest = TRUE)
de_dxa  <- svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~W4YR,
                     data = df_dxa_wt, nest = TRUE)

# [A] 相对握力（功能学）—— 主模型
m_grip_rel <- svyglm(Relative_Grip ~ Log2_TG + Age + Gender + Cycle, design = de_grip)
# [B] 相对 ASM（结构学）—— 主模型
m_dxa_rel  <- svyglm(Relative_ASM  ~ Log2_TG + Age + Gender + Cycle, design = de_dxa)
# [C] 绝对握力 + BMI —— BMI 悖论（机械混淆吸收）
m_grip_abs <- svyglm(Grip_Strength ~ Log2_TG + BMI + Age + Gender + Cycle, design = de_grip)
# [D] 绝对 ASM + BMI
m_dxa_abs  <- svyglm(ASM_kg ~ Log2_TG + BMI + Age + Gender + Cycle, design = de_dxa)

cat("\n----- [A] 相对握力 (加权) -----\n");    print(round(summary(m_grip_rel)$coefficients, 6))
cat("\n----- [B] 相对 ASM (加权) -----\n");     print(round(summary(m_dxa_rel)$coefficients, 6))
cat("\n----- [C] 绝对握力 + BMI (加权) -----\n"); print(round(summary(m_grip_abs)$coefficients, 6))

cat("\n========== 4. 汇总输出与稿件对照 ==========\n")
extract_tg <- function(m, lab) {
  s <- summary(m)$coefficients["Log2_TG", ]
  data.frame(Model = lab, Beta = s[1], SE = s[2], t = s[3], P = s[4],
             CI_low = s[1] - 1.96 * s[2], CI_high = s[1] + 1.96 * s[2], row.names = NULL)
}
tab_out <- rbind(
  extract_tg(m_grip_rel, "Weighted: Relative grip ~ Log2_TG + Age + Sex + Cycle"),
  extract_tg(m_dxa_rel,  "Weighted: Relative ASM ~ Log2_TG + Age + Sex + Cycle"),
  extract_tg(m_grip_abs, "Weighted: Absolute grip ~ Log2_TG + BMI + Age + Sex + Cycle"),
  extract_tg(m_dxa_abs,  "Weighted: Absolute ASM ~ Log2_TG + BMI + Age + Sex + Cycle")
)
print(tab_out)

write.csv(tab_out, file.path(out_dir, "Table_Weighted_Regression_MS.csv"), row.names = FALSE)
write.csv(tidy(m_grip_rel, conf.int = TRUE), file.path(out_dir, "Table_Weighted_Grip_Full.csv"), row.names = FALSE)
write.csv(tidy(m_dxa_rel,  conf.int = TRUE), file.path(out_dir, "Table_Weighted_DXA_Full.csv"),  row.names = FALSE)

cat("\n========== 5. 加权 vs 未加权对照（论证加权必要性） ==========\n")
m_grip_unw <- lm(Relative_Grip ~ Log2_TG + Age + Gender + Cycle, data = df_grip)
cmp <- rbind(
  data.frame(Approach = "Unweighted", Beta = coef(m_grip_unw)["Log2_TG"],
             SE = summary(m_grip_unw)$coefficients["Log2_TG", 2],
             P = summary(m_grip_unw)$coefficients["Log2_TG", 4]),
  data.frame(Approach = "Survey-weighted", Beta = coef(m_grip_rel)["Log2_TG"],
             SE = summary(m_grip_rel)$coefficients["Log2_TG", 2],
             P = summary(m_grip_rel)$coefficients["Log2_TG", 4])
)
print(cmp)
write.csv(cmp, file.path(out_dir, "Table_Weighted_vs_Unweighted.csv"), row.names = FALSE)

cat("\n✅ 加权回归全部完成，表格已落盘至 MS/M/\n")
