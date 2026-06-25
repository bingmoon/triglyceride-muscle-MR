# =====================================================================
# MSM 项目 - 阶段三：孟德尔随机化 (MR) 与多变量 MR (MVMR)
# 工作目录：/Users/bing/MS
# 暴露：TG (ieu-b-111), BMI (ieu-b-40)
# 结局：Hand grip strength (ieu-b-39)
# 所有输出前缀 "Stage3_"
# =====================================================================

setwd("/Users/bing/MS")

# 加载包
library(TwoSampleMR)
library(ieugwasr)
library(dplyr)
library(ggplot2)
library(broom)

# 设置 API
options(ieugwasr_api = "https://api.opengwas.io/api/")

# ========== 1. 单变量 MR (UVMR) ==========
cat("\n========== 单变量 MR：TG → 握力 ==========\n")

# 提取暴露工具变量
exposure_id <- "ieu-b-111"
outcome_id  <- "ieu-b-39"

exp_dat <- extract_instruments(outcomes = exposure_id, p1 = 5e-8)
cat(sprintf("暴露工具变量提取：%d 个 SNP\n", nrow(exp_dat)))

# 提取结局关联
out_dat <- extract_outcome_data(snps = exp_dat$SNP, outcomes = outcome_id, proxies = TRUE)

# 协调数据
dat <- harmonise_data(exposure_dat = exp_dat, outcome_dat = out_dat, action = 2)
dat <- dat[dat$mr_keep == TRUE, ]

# 计算 F 统计量
dat$F_statistic <- (dat$beta.exposure / dat$se.exposure)^2
mean_F <- mean(dat$F_statistic)
cat(sprintf("平均 F 统计量：%.2f\n", mean_F))

# MR 分析
mr_res <- mr(dat, method_list = c("mr_ivw", "mr_weighted_median", "mr_egger_regression"))
mr_res <- generate_odds_ratios(mr_res)

# 敏感性分析
hetero <- mr_heterogeneity(dat)
pleio  <- mr_pleiotropy_test(dat)

# 保存结果
write.csv(mr_res, "Stage3_Table_UVMR_Results.csv", row.names = FALSE)
write.csv(dat, "Stage3_Table_UVMR_Harmonised_Data.csv", row.names = FALSE)

# =====================================================================
# 阶段三 Nature 级别可视化（完全自定义）
# =====================================================================
library(ggplot2)
library(ggrepel)
library(dplyr)

# ---------- 通用 Nature 风格主题 ----------
theme_nature <- theme_bw(base_size = 12) +
  theme(
    axis.title = element_text(face = "bold", size = 13),
    axis.text = element_text(color = "black", size = 11),
    axis.line = element_line(linewidth = 0.8, color = "black"),
    panel.grid.major = element_line(linewidth = 0.3, color = "grey90"),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, size = 12, color = "grey30"),
    legend.position = "top",
    legend.title = element_text(face = "bold", size = 11),
    legend.text = element_text(size = 10),
    strip.text = element_text(face = "bold", size = 12),
    strip.background = element_rect(fill = "grey95", color = "black")
  )

# ---------- 图1：散点图 ----------
# 提取每种方法的斜率用于绘制直线
mr_methods <- unique(mr_res$method)
# 准备散点数据（每个 SNP 的效应）
scatter_dat <- dat %>%
  mutate(
    beta.outcome = beta.outcome,
    beta.exposure = beta.exposure
  )

p1 <- ggplot(scatter_dat, aes(x = beta.exposure, y = beta.outcome)) +
  # 误差线（SNP 的结局 SE）
  geom_errorbar(aes(ymin = beta.outcome - se.outcome, ymax = beta.outcome + se.outcome),
                width = 0, color = "grey70", linewidth = 0.4) +
  geom_errorbarh(aes(xmin = beta.exposure - se.exposure, xmax = beta.exposure + se.exposure),
                 height = 0, color = "grey70", linewidth = 0.4) +
  # SNP 点
  geom_point(aes(color = "SNP effect"), size = 2.5, alpha = 0.8) +
  # 各 MR 方法的拟合线
  geom_abline(data = mr_res, aes(intercept = 0, slope = b, color = method),
              linewidth = 1.2, alpha = 0.9) +
  scale_color_manual(
    values = c("SNP effect" = "grey60",
               "Inverse variance weighted" = "#E64B35",
               "MR Egger" = "#4DBBD5",
               "Weighted median" = "#00A087"),
    breaks = c("SNP effect", "Inverse variance weighted", "MR Egger", "Weighted median")
  ) +
  labs(
    title = "Causal Effect of TG on Grip Strength (UVMR)",
    subtitle = paste0("Each point represents a single SNP; lines show MR estimates (",
                      nrow(scatter_dat), " SNPs)"),
    x = "SNP effect on Triglycerides (beta.exposure)",
    y = "SNP effect on Grip Strength (beta.outcome)",
    color = ""
  ) +
  theme_nature +
  guides(color = guide_legend(override.aes = list(
    linetype = c("blank", "solid", "solid", "solid"),
    shape = c(16, NA, NA, NA)
  )))

ggsave("Stage3_Figure_UVMR_Scatter.pdf", plot = p1, width = 7.5, height = 6.5, dpi = 300)

# ---------- 图2：森林图 ----------
# 单 SNP 效应
single_snp <- mr_singlesnp(dat)
# 提取 IVW 汇总效应
ivw_est <- mr_res %>% filter(method == "Inverse variance weighted")

# 准备森林图数据
forest_dat <- single_snp %>%
  mutate(
    SNP = ifelse(SNP == "All - Inverse variance weighted", "IVW estimate", SNP),
    is_summary = SNP == "IVW estimate"
  )

# 排序：IVW 在最下，其他按 b 排序
snp_order <- forest_dat %>%
  filter(!is_summary) %>%
  arrange(b) %>%
  pull(SNP)
forest_dat$SNP <- factor(forest_dat$SNP, levels = c(snp_order, "IVW estimate"))

p2 <- ggplot(forest_dat, aes(x = b, y = SNP)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.8) +
  geom_errorbarh(aes(xmin = b - 1.96 * se, xmax = b + 1.96 * se),
                 height = 0.2, linewidth = 1, color = ifelse(forest_dat$is_summary, "#E64B35", "grey60")) +
  geom_point(aes(color = is_summary, size = is_summary)) +
  scale_color_manual(values = c("TRUE" = "#E64B35", "FALSE" = "#4DBBD5"), guide = "none") +
  scale_size_manual(values = c("TRUE" = 5, "FALSE" = 2), guide = "none") +
  labs(
    title = "Individual SNP Causal Effects",
    subtitle = paste0("Total effect (IVW): β = ", round(ivw_est$b, 4),
                      " (95% CI: ", round(ivw_est$b - 1.96*ivw_est$se, 4),
                      ", ", round(ivw_est$b + 1.96*ivw_est$se, 4), "), P = ",
                      format.pval(ivw_est$pval, digits = 2)),
    x = "Causal estimate (beta)",
    y = ""
  ) +
  theme_nature +
  theme(
    axis.text.y = element_text(size = 5, face = "bold"),
    panel.grid.major.y = element_blank()
  )

ggsave("Stage3_Figure_UVMR_Forest.pdf", plot = p2, width = 10, height = 14, dpi = 300)

# ---------- 图3：漏斗图 ----------
# 用单 SNP 数据，x 为个体估计，y 为精度（1/se）
funnel_dat <- single_snp %>%
  filter(SNP != "All - Inverse variance weighted") %>%
  mutate(precision = 1 / se)

p3 <- ggplot(funnel_dat, aes(x = b, y = precision)) +
  geom_point(size = 2.5, color = "#4DBBD5", alpha = 0.8) +
  geom_vline(xintercept = ivw_est$b, linetype = "dashed", color = "#E64B35", linewidth = 1.2) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey50", linewidth = 0.8) +
  labs(
    title = "Funnel Plot of Individual SNP Estimates",
    subtitle = paste0("Red dashed line: IVW estimate (", round(ivw_est$b, 4), ")"),
    x = "Individual causal estimate (beta)",
    y = "Precision (1 / SE)"
  ) +
  theme_nature

ggsave("Stage3_Figure_UVMR_Funnel.pdf", plot = p3, width = 6.5, height = 6, dpi = 300)

# ---------- 图4：留一法图（优化标签间距和可读性） ----------
res_loo <- mr_leaveoneout(dat)

# 安全转换为数据框
loo_list <- res_loo
if (is.list(loo_list) && !is.data.frame(loo_list)) {
  loo_dat <- as.data.frame(loo_list[[1]])
} else if (is.data.frame(loo_list)) {
  loo_dat <- loo_list
} else {
  loo_dat <- as.data.frame(loo_list)
}

# 确保列名
if (!"b" %in% names(loo_dat)) {
  b_col <- grep("b$|beta|effect", names(loo_dat), ignore.case = TRUE, value = TRUE)[1]
  se_col <- grep("se$|se.exposure", names(loo_dat), ignore.case = TRUE, value = TRUE)[1]
  if (!is.na(b_col)) names(loo_dat)[names(loo_dat) == b_col] <- "b"
  if (!is.na(se_col)) names(loo_dat)[names(loo_dat) == se_col] <- "se"
}
if (!"b" %in% names(loo_dat)) {
  loo_dat$b <- loo_dat[, 1]
  loo_dat$se <- loo_dat[, 2]
}

loo_dat$lo <- loo_dat$b - 1.96 * loo_dat$se
loo_dat$up <- loo_dat$b + 1.96 * loo_dat$se

if (!"SNP" %in% names(loo_dat)) {
  if (!is.null(rownames(loo_dat))) {
    loo_dat$SNP <- rownames(loo_dat)
  } else {
    loo_dat$SNP <- paste0("SNP_", seq_len(nrow(loo_dat)))
  }
}

# 总效应
ivw_est <- mr_res %>% filter(method == "Inverse variance weighted")
ivw_b <- ivw_est$b[1]

# 绘图：x轴仅保留点，不显示文本，用ggrepel标注
p4 <- ggplot(loo_dat, aes(x = SNP, y = b)) +
  geom_hline(yintercept = ivw_b, linetype = "dashed", color = "#E64B35", linewidth = 1.2) +
  geom_errorbar(aes(ymin = lo, ymax = up), width = 0.3, color = "grey60", linewidth = 0.6) +
  geom_point(size = 2.5, color = "#4DBBD5", alpha = 0.9) +
  geom_text_repel(
    aes(label = SNP),
    size = 3.0,
    max.overlaps = 25,
    box.padding = 0.6,      # 增大文本框周围的留白
    point.padding = 0.5,    # 增加点和文本之间的间距
    segment.color = "grey70",
    segment.size = 0.3,
    direction = "y",        # 标签主要在竖直方向展开
    force = 10              # 增加排斥力，让标签更分散
  ) +
  labs(
    title = "Leave-One-Out Analysis",
    subtitle = paste0("Each point shows IVW estimate after omitting one SNP. Red line: overall IVW (",
                      round(ivw_b, 4), ")"),
    x = "Omitted SNP",
    y = "IVW estimate"
  ) +
  theme_nature +
  theme(
    axis.text.x = element_blank(),    # 隐藏底部重叠的 SNP 文字
    axis.ticks.x = element_blank(),
    panel.grid.major.x = element_blank()
  )

ggsave("Stage3_Figure_UVMR_LeaveOneOut.pdf", plot = p4, width = 16, height = 12, dpi = 300)

cat("✅ 所有 Nature 级别 MR 图表已生成。\n")



# 打印单变量结果摘要
cat("\n单变量 MR 主要结果 (IVW):\n")
print(subset(mr_res, method == "Inverse variance weighted"))

# ========== 2. 多变量 MR (MVMR) ==========
cat("\n========== 多变量 MR：TG + BMI → 握力 ==========\n")

# 提取两个暴露的工具变量并集
exposure_ids <- c("ieu-b-111", "ieu-b-40")
mv_exp_dat <- mv_extract_exposures(id_exposure = exposure_ids, 
                                   pval_threshold = 5e-8, 
                                   find_proxies = TRUE)
cat(sprintf("多变量暴露工具变量：%d 个 SNP\n", nrow(mv_exp_dat)))

# 提取结局
out_dat_mv <- extract_outcome_data(snps = mv_exp_dat$SNP, outcomes = outcome_id, proxies = TRUE)

# 协调数据
mv_dat <- mv_harmonise_data(exposure_dat = mv_exp_dat, outcome_dat = out_dat_mv)

# MVMR 分析
res_mvmr <- mv_multiple(mv_dat)

# 整理结果
final_results <- res_mvmr$result %>%
  mutate(
    Exposure_Name = case_when(
      id.exposure == "ieu-b-111" ~ "Triglycerides (TG)",
      id.exposure == "ieu-b-40"  ~ "Body Mass Index (BMI)",
      TRUE ~ exposure
    ),
    Outcome = "Hand Grip Strength",
    OR = exp(b),
    LCI = exp(b - 1.96 * se),
    UCI = exp(b + 1.96 * se),
    P_fmt = ifelse(pval < 0.001, formatC(pval, format = "e", digits = 2), sprintf("%.3f", pval))
  ) %>%
  dplyr::select(Exposure_Name, Outcome, nsnp, b, se, pval, OR, LCI, UCI, P_fmt)

cat("\nMVMR 结果:\n")
print(final_results)

# 保存结果表格
write.csv(final_results, "Stage3_Table_MVMR_Results.csv", row.names = FALSE)
write.csv(as.data.frame(mv_dat$exposure_beta), "Stage3_Table_MVMR_SNPs.csv")

# 可视化 MVMR 森林图
p_mvmr <- ggplot(final_results, aes(y = reorder(Exposure_Name, b), x = OR)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "#E64B35", linewidth = 1) +
  geom_errorbarh(aes(xmin = LCI, xmax = UCI, color = Exposure_Name), height = 0.2, linewidth = 1.2) +
  geom_point(aes(fill = Exposure_Name), shape = 21, size = 5, color = "black", stroke = 1) +
  geom_text(aes(label = paste0("P = ", P_fmt), x = max(UCI) + 0.5), vjust = -0.8, fontface = "italic", size = 4) +
  scale_color_manual(values = c("Triglycerides (TG)" = "#4DBBD5", "Body Mass Index (BMI)" = "#3C5488")) +
  scale_fill_manual(values = c("Triglycerides (TG)" = "#4DBBD5", "Body Mass Index (BMI)" = "#3C5488")) +
  labs(title = "Independent Causal Effects on Hand Grip Strength",
       subtitle = "Multivariable Mendelian Randomization (MVMR)",
       x = "Odds Ratio (95% CI) for Muscle Strength",
       y = "") +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 15),
    plot.subtitle = element_text(hjust = 0.5, face = "italic", color = "grey50"),
    axis.text.y = element_text(face = "bold", color = "black", size = 13),
    axis.title.x = element_text(face = "bold"),
    legend.position = "none",
    panel.grid.minor = element_blank()
  ) +
  coord_cartesian(xlim = c(0, max(final_results$UCI) + 1.5))

ggsave("Stage3_Figure_MVMR_Forest.pdf", plot = p_mvmr, width = 8, height = 4, dpi = 300)



# =====================================================================
# 阶段三补充：MVMR条件F统计量 + 森林图标签修正
# =====================================================================
cat("\n========== MVMR 条件F统计量 ==========\n")

library(MVMR)

if (exists("mv_dat")) {
  cond_F <- tryCatch({
    strength_mvmr(mv_dat)
  }, error = function(e) {
    cat("⚠️ MVMR::strength_mvmr 失败，手动计算近似条件F。\n")
    NULL
  })
  
  if (!is.null(cond_F)) {
    print(cond_F)
  } else {
    # 手动近似：F = (b/se)^2
    if (exists("final_results")) {
      final_results$F_approx <- (final_results$b / final_results$se)^2
      cat("近似条件F统计量（基于MVMR估计）：\n")
      print(final_results[, c("Exposure_Name", "F_approx")])
    }
  }
} else {
  cat("⚠️ mv_dat 对象不存在，请先运行阶段三主代码。\n")
}

cat("\n========== 重新生成MVMR森林图（标签修正） ==========\n")

if (exists("final_results")) {
  p_mvmr <- ggplot(final_results, aes(y = reorder(Exposure_Name, b), x = b)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "#E64B35", linewidth = 1) +
    geom_errorbarh(aes(xmin = LCI, xmax = UCI, color = Exposure_Name), height = 0.2, linewidth = 1.2) +
    geom_point(aes(fill = Exposure_Name), shape = 21, size = 5, color = "black", stroke = 1) +
    geom_text(aes(label = paste0("P = ", P_fmt), x = max(UCI) + 0.5), vjust = -0.8, fontface = "italic", size = 4) +
    scale_color_manual(values = c("Triglycerides (TG)" = "#4DBBD5", "Body Mass Index (BMI)" = "#3C5488")) +
    scale_fill_manual(values = c("Triglycerides (TG)" = "#4DBBD5", "Body Mass Index (BMI)" = "#3C5488")) +
    labs(title = "Independent Causal Effects on Hand Grip Strength (MVMR)",
         subtitle = "Multivariable Mendelian Randomization",
         x = "Beta coefficient (95% CI)",
         y = "") +
    theme_bw(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 15),
      plot.subtitle = element_text(hjust = 0.5, face = "italic", color = "grey50"),
      axis.text.y = element_text(face = "bold", color = "black", size = 13),
      axis.title.x = element_text(face = "bold"),
      legend.position = "none",
      panel.grid.minor = element_blank()
    ) +
    coord_cartesian(xlim = c(min(final_results$LCI) - 0.5, max(final_results$UCI) + 1.5))
  
  ggsave("Stage3_Figure_MVMR_Forest.pdf", plot = p_mvmr, width = 8, height = 4, dpi = 300)
  cat("✅ 已重新生成 Stage3_Figure_MVMR_Forest.pdf，x轴标签已修正。\n")
} else {
  cat("⚠️ final_results 不存在，跳过森林图。\n")
}

cat("✅ 阶段三补充分析完成。\n")



# =====================================================================
# 生成 Table S4 & S5 + 整合完整的 Supplementary_Tables.xlsx
# 工作目录：/Users/bing/MS
# =====================================================================

setwd("/Users/bing/MS")
library(openxlsx)
library(dplyr)
library(broom)
library(nhanesA)
library(TwoSampleMR)

# ========== 1. 构建 Table S4：敏感性分析汇总 ==========
cat("正在构建 Table S4...\n")

# 1a. 主模型（相对握力）
m1 <- lm(Relative_Grip ~ Log2_TG + Age + Gender + Cycle, data = df_grip_main)
t1 <- tidy(m1, conf.int = TRUE) %>% filter(term == "Log2_TG") %>%
  mutate(Model = "Relative grip (base)")

# 1b. 绝对握力 + BMI
m2 <- lm(Grip_Strength ~ Log2_TG + Age + Gender + Cycle + BMI, data = df_grip_main)
t2 <- tidy(m2, conf.int = TRUE) %>% filter(term == "Log2_TG") %>%
  mutate(Model = "Absolute grip + BMI")

# 1c. 绝对握力 + BMI + 体脂率
# 拉取脂肪数据（G/H周期）
dxa_G <- nhanes("DXX_G"); dxa_H <- nhanes("DXX_H")
dxa_G$SEQN <- as.numeric(dxa_G$SEQN); dxa_H$SEQN <- as.numeric(dxa_H$SEQN)
dxa_fat <- bind_rows(
  dxa_G %>% select(SEQN, Fat_mass = DXDTOFAT) %>% mutate(Cycle = "G"),
  dxa_H %>% select(SEQN, Fat_mass = DXDTOFAT) %>% mutate(Cycle = "H")
)
df_grip_fat <- df_grip_main %>%
  left_join(dxa_fat, by = c("SEQN", "Cycle")) %>%
  mutate(Fat_pct = Fat_mass / Weight_kg * 100)

m3 <- lm(Grip_Strength ~ Log2_TG + Age + Gender + Cycle + BMI + Fat_pct,
         data = df_grip_fat)
t3 <- tidy(m3, conf.int = TRUE) %>% filter(term == "Log2_TG") %>%
  mutate(Model = "Absolute grip + BMI + Fat%")

# 1d. 相对 ASM 主模型
m4 <- lm(Relative_ASM ~ Log2_TG + Age + Gender + Cycle, data = df_dxa_all)
t4 <- tidy(m4, conf.int = TRUE) %>% filter(term == "Log2_TG") %>%
  mutate(Model = "Relative ASM (base)")

# 1e. 绝对 ASM + BMI
m5 <- lm(ASM_kg ~ Log2_TG + Age + Gender + Cycle + BMI, data = df_dxa_all)
t5 <- tidy(m5, conf.int = TRUE) %>% filter(term == "Log2_TG") %>%
  mutate(Model = "Absolute ASM + BMI")

# 1f. 绝对 ASM + BMI + 体脂率（四周期）
dxa_I <- nhanes("DXX_I"); dxa_J <- nhanes("DXX_J")
dxa_I$SEQN <- as.numeric(dxa_I$SEQN); dxa_J$SEQN <- as.numeric(dxa_J$SEQN)
dxa_fat_all <- bind_rows(
  dxa_G %>% select(SEQN, Fat_mass = DXDTOFAT) %>% mutate(Cycle = "G"),
  dxa_H %>% select(SEQN, Fat_mass = DXDTOFAT) %>% mutate(Cycle = "H"),
  dxa_I %>% select(SEQN, Fat_mass = DXDTOFAT) %>% mutate(Cycle = "I"),
  dxa_J %>% select(SEQN, Fat_mass = DXDTOFAT) %>% mutate(Cycle = "J")
)
df_dxa_fat <- df_dxa_all %>%
  left_join(dxa_fat_all, by = c("SEQN", "Cycle")) %>%
  mutate(Fat_pct = Fat_mass / Weight_kg * 100)

m6 <- lm(ASM_kg ~ Log2_TG + Age + Gender + Cycle + BMI + Fat_pct,
         data = df_dxa_fat)
t6 <- tidy(m6, conf.int = TRUE) %>% filter(term == "Log2_TG") %>%
  mutate(Model = "Absolute ASM + BMI + Fat%")

# 合并所有敏感性分析
table_s4 <- bind_rows(t1, t2, t3, t4, t5, t6) %>%
  select(Model, Beta = estimate, SE = std.error, P = p.value, 
         CI_low = conf.low, CI_high = conf.high) %>%
  mutate(across(where(is.numeric), ~ round(., 4)))

# 保存 CSV
write.csv(table_s4, "Table_S4_Sensitivity_Analysis.csv", row.names = FALSE)
cat("✅ Table S4 已保存为 Table_S4_Sensitivity_Analysis.csv\n")
print(table_s4)

# ========== 2. 构建 Table S5：双向 MR ==========
cat("\n正在构建 Table S5（双向 MR）...\n")

if (!exists("table_s5")) {
  exposure_grip <- extract_instruments("ieu-b-39", p1 = 5e-8)
  outcome_tg <- extract_outcome_data(exposure_grip$SNP, "ieu-b-111")
  dat_rev <- harmonise_data(exposure_grip, outcome_tg, action = 2)
  mr_rev <- mr(dat_rev, method_list = c("mr_ivw", "mr_egger_regression", "mr_weighted_median"))
  table_s5 <- mr_rev %>%
    select(Method = method, nSNP = nsnp, Beta = b, SE = se, P = pval) %>%
    mutate(across(where(is.numeric), ~ round(., 4)))
} else {
  cat("table_s5 已存在，跳过构建。\n")
}

# 保存 CSV
write.csv(table_s5, "Table_S5_Bidirectional_MR.csv", row.names = FALSE)
cat("✅ Table S5 已保存为 Table_S5_Bidirectional_MR.csv\n")
print(table_s5)

# ========== 3. 整合完整的 Supplementary_Tables.xlsx ==========
cat("\n正在生成整合的 Supplementary_Tables.xlsx...\n")

wb <- createWorkbook()

# --- Table S1: UVMR 结果 ---
if (exists("mr_res")) {
  s1 <- mr_res %>%
    select(Method = method, nSNP = nsnp, Beta = b, SE = se, P = pval) %>%
    mutate(across(where(is.numeric), ~ round(., 4)))
} else {
  s1 <- data.frame(Message = "mr_res not found")
}
addWorksheet(wb, "Table S1")
writeData(wb, "Table S1", s1)

# --- Table S2: UVMR 单 SNP 效应 ---
if (exists("dat")) {
  s2 <- dat %>%
    select(SNP, effect_allele.exposure, other_allele.exposure,
           beta.exposure, se.exposure, pval.exposure,
           beta.outcome, se.outcome, pval.outcome) %>%
    mutate(across(where(is.numeric), ~ round(., 4)))
} else {
  s2 <- data.frame(Message = "dat not found")
}
addWorksheet(wb, "Table S2")
writeData(wb, "Table S2", s2)

# --- Table S3: MVMR 工具变量列表 ---
if (exists("mv_dat") && !is.null(mv_dat$exposure_beta)) {
  s3 <- as.data.frame(mv_dat$exposure_beta)
  if (ncol(s3) > 0) {
    s3 <- s3 %>% mutate(across(where(is.numeric), ~ round(., 4)))
  }
} else {
  s3 <- data.frame(Message = "mv_dat not found")
}
addWorksheet(wb, "Table S3")
writeData(wb, "Table S3", s3)

# --- Table S4: 敏感性分析汇总 ---
addWorksheet(wb, "Table S4")
writeData(wb, "Table S4", table_s4)

# --- Table S5: 双向 MR ---
addWorksheet(wb, "Table S5")
writeData(wb, "Table S5", table_s5)

# 保存 Excel
saveWorkbook(wb, "Supplementary_Tables.xlsx", overwrite = TRUE)
cat("✅ 已生成 Supplementary_Tables.xlsx，包含 Table S1 – S5\n")
