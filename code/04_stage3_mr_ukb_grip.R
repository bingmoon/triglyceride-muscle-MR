# ---
# Archived from: 4.MSM（ST003803） 1.md  (L1282-1557, block 4)
# Body md5: 4b91075e6aa8b360fbe0be1949c7dcb6   | 3-line header added by release tooling
# ---


# =====================================================================
# MSM 项目 - 阶段三：孟德尔随机化 (MR) 绝对纯净全流程脚本
# 暴露：TG (ieu-b-111), BMI (ieu-b-40)
# 真实结局：Hand grip strength (ukb-b-10215, UK Biobank)
# 特性：不依赖任何本地历史数据，从零生成所有图表与表格
# =====================================================================

# ---------- 0. 环境与目录初始化 ----------
setwd("/Users/bing/MS")
out_dir <- "/Users/bing/MS/M"
dir.create(out_dir, showWarnings = FALSE)

suppressMessages({
  library(TwoSampleMR)
  library(ieugwasr)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
})
options(ieugwasr_api = "https://api.opengwas.io/api/")
set.seed(2026) # 保证多变量排版的绝对可重复性

cat("\n========== 1. 提取暴露与真实结局数据 ==========\n")
# 1.1 提取 TG 和 BMI 暴露工具变量
cat("提取暴露数据...\n")
exp_dat <- extract_instruments("ieu-b-111", p1 = 5e-8)
mv_exp_dat <- mv_extract_exposures(c("ieu-b-111", "ieu-b-40"), pval_threshold = 5e-8, find_proxies = TRUE)

# 1.2 提取正确结局 (握力: ukb-b-10215)
cat("提取结局数据 (ukb-b-10215)...\n")
outcome_grip <- extract_outcome_data(snps = unique(mv_exp_dat$SNP), outcomes = "ukb-b-10215", proxies = TRUE)

# 1.3 数据对齐 (Harmonise)
dat <- harmonise_data(exposure_dat = exp_dat, outcome_dat = outcome_grip)
mv_dat <- mv_harmonise_data(mv_exp_dat, outcome_grip)

# 计算原始 F 统计量
dat$F_statistic <- (dat$beta.exposure / dat$se.exposure)^2


cat("\n========== 2. 生成 Table S6：PhenoScanner 混杂在线排查 ==========\n")
cat("正在通过 API 实时查询 SNP 关联表型 (此步可能耗时 1-2 分钟，请耐心等待)...\n")

snps_to_check <- unique(dat$SNP)
results_list <- list()
snps_batches <- split(snps_to_check, ceiling(seq_along(snps_to_check) / 50))

# 加入容错机制的在线抓取
for(i in seq_along(snps_batches)) {
  res <- tryCatch(
    phewas(variants = snps_batches[[i]], pval = 1e-5), 
    error = function(e) { cat(sprintf("批次 %d 查询超时，跳过\n", i)); return(NULL) }
  )
  if(!is.null(res)) results_list[[i]] <- res
  Sys.sleep(1) # 防止 API 频率限制
}
all_phewas_results <- bind_rows(results_list)

# 匹配论文预设的混杂关键字
confounder_keywords <- c("physical activity", "smoking", "alcohol", "diabetes", "glucose", "insulin", "hypertension", "HDL cholesterol", "LDL cholesterol", "triglycerides", "body mass index", "BMI", "waist circumference")

# 筛选显著混杂 SNPs (P < 5e-6)
sig_confounders_unique <- all_phewas_results %>% 
  filter(p < 5e-6) %>%
  filter(grepl(paste(confounder_keywords, collapse = "|"), trait, ignore.case = TRUE)) %>% 
  group_by(rsid) %>% slice_min(p, n = 1) %>% ungroup()

confounder_snps <- unique(sig_confounders_unique$rsid)
cat(sprintf("查杀完成！发现 %d 个潜在混杂 SNP。\n", length(confounder_snps)))

# 👉 输出 Table S6
write.csv(sig_confounders_unique, file.path(out_dir, "Table_S6_Excluded_Confounder_SNPs.csv"), row.names = FALSE)


cat("\n========== 3. MR 因果推断 (原始与纯净版) ==========\n")
# 3.1 原始模型 (含混杂)
mr_res_raw <- mr(dat, method_list = c("mr_ivw", "mr_weighted_median", "mr_egger_regression"))
res_mvmr_raw <- mv_multiple(mv_dat)

# 3.2 纯净模型 (剔除混杂)
# 采用安全的 Base R 子集过滤，防止格式崩塌
dat_clean <- subset(dat, !(SNP %in% confounder_snps))
mr_res_clean <- mr(dat_clean, method_list = c("mr_ivw", "mr_weighted_median", "mr_egger_regression"))

mv_exp_clean <- subset(mv_exp_dat, !(SNP %in% confounder_snps))
out_clean <- subset(outcome_grip, !(SNP %in% confounder_snps))
mv_dat_clean <- mv_harmonise_data(exposure_dat = mv_exp_clean, outcome_dat = out_clean)
res_mvmr_clean <- mv_multiple(mv_dat_clean)

# 👉 输出 Table S7 (Clean UVMR) 与 Table S8 (Clean MVMR)
write.csv(mr_res_clean, file.path(out_dir, "Table_S7_UVMR_Clean_Results.csv"), row.names = FALSE)
write.csv(res_mvmr_clean$result, file.path(out_dir, "Table_S8_MVMR_Clean_Results.csv"), row.names = FALSE)


cat("\n========== 4. 敏感性分析 (异质性与多效性) ==========\n")
# 基于 Clean 数据的敏感性检验
hetero_clean <- mr_heterogeneity(dat_clean)
pleio_clean <- mr_pleiotropy_test(dat_clean)

# 👉 输出 Table S9 (Heterogeneity) 与 Table S10 (Pleiotropy)
write.csv(hetero_clean, file.path(out_dir, "Table_S9_Heterogeneity_Clean.csv"), row.names = FALSE)
write.csv(pleio_clean, file.path(out_dir, "Table_S10_Pleiotropy_Clean.csv"), row.names = FALSE)


cat("\n========== 5. 生成 Table S11：Top 10 核心 SNPs (生物学机制验证) ==========\n")
# 提取 F 统计量最高的 10 个纯净 SNPs 供文献核对靶基因
top_10_snps <- dat_clean %>% 
  arrange(desc(F_statistic)) %>% 
  head(10) %>% 
  dplyr::select(SNP, chr.exposure, pos.exposure, effect_allele.exposure, beta.exposure, pval.exposure, F_statistic)

# 👉 输出 Table S11
write.csv(top_10_snps, file.path(out_dir, "Table_S11_Top10_SNPs_By_Fstat.csv"), row.names = FALSE)


cat("\n========== 6. 顶级医学期刊标准高阶可视化 (基于 Clean 数据) ==========\n")
theme_nature <- theme_classic(base_size = 13) + 
  theme(
    axis.title = element_text(face = "bold", size = 14, color = "black"),
    axis.text = element_text(color = "black", size = 12),
    axis.line = element_line(linewidth = 1, color = "black"),
    plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, size = 12, color = "grey40"),
    legend.position = "top", legend.title = element_blank(),
    panel.grid = element_blank()
  )

sci_colors <- c("Inverse variance weighted" = "#E64B35", "MR Egger" = "#4DBBD5", "Weighted median" = "#00A087", "SNP effect" = "grey50")



# =====================================================================
# MSM 项目 - 阶段 6 及以后：顶级期刊标准可视化与附属图表补全
# 备注：运行前请确保当前环境中已存在 dat, dat_clean, mr_res_clean, mr_res_raw, mv_exp_dat 等基础数据
# =====================================================================
library(ggplot2)
library(ggrepel)
library(dplyr)
library(TwoSampleMR)

out_dir <- "/Users/bing/MS/M"

# 提取 Clean 数据的 IVW 结果，供后续图表基准线使用
ivw_est <- mr_res_clean %>% filter(method == "Inverse variance weighted")

cat("\n========== 6.1 散点图 (Fig S4 - 极简防弹版) ==========\n")
# 强制清理内存中残留的旧图表对象，防止图例冲突
if(exists("p_scatter")) rm(p_scatter)

p_scatter <- ggplot(dat_clean, aes(x = beta.exposure, y = beta.outcome)) +
  geom_errorbar(aes(ymin = beta.outcome - se.outcome, ymax = beta.outcome + se.outcome), color = "grey80", width = 0) +
  geom_errorbar(aes(xmin = beta.exposure - se.exposure, xmax = beta.exposure + se.exposure), color = "grey80", width = 0) +
  geom_point(color = "grey50", size = 3, alpha = 0.8) +
  geom_abline(data = mr_res_clean, aes(slope = b, color = method), intercept = 0, linewidth = 1.2) +
  scale_color_manual(values = c("Inverse variance weighted" = "#E64B35", 
                                "MR Egger" = "#4DBBD5", 
                                "Weighted median" = "#00A087")) +
  labs(title = "Causal Effect of TG on Grip Strength", 
       subtitle = sprintf("Cleaned Instruments (N = %d SNPs)", nrow(dat_clean)), 
       x = "SNP effect on Triglycerides", 
       y = "SNP effect on Grip Strength") +
  theme_classic(base_size = 13) +
  theme(legend.title = element_blank(), legend.position = "top")

ggsave(file.path(out_dir, "Fig_S4_Scatter_Clean_Nature.pdf"), plot = p_scatter, width = 7.5, height = 6.5)
cat("✅ Fig S4 散点图已成功生成并保存！\n")


cat("\n========== 6.2 漏斗图 (Fig S6) ==========\n")
funnel_dat <- mr_singlesnp(dat_clean) %>% filter(SNP != "All - Inverse variance weighted") %>% mutate(precision = 1 / se)

p_funnel <- ggplot(funnel_dat, aes(x = b, y = precision)) +
  geom_point(size = 3, color = "#4DBBD5", alpha = 0.8, shape = 16) +
  geom_vline(xintercept = ivw_est$b, linetype = "dashed", color = "#E64B35", linewidth = 1.2) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "black", linewidth = 0.8) +
  labs(title = "Funnel Plot of Individual SNP Estimates", 
       subtitle = sprintf("Dashed line represents IVW estimate (%.4f)", ivw_est$b), 
       x = "Individual causal estimate (beta)", 
       y = "Precision (1 / SE)") +
  theme_classic(base_size = 13) + 
  theme(plot.title = element_text(face="bold", hjust=0.5), plot.subtitle = element_text(hjust=0.5))

ggsave(file.path(out_dir, "Fig_S6_Funnel_Clean_Nature.pdf"), plot = p_funnel, width = 7, height = 6.5)
cat("✅ Fig S6 漏斗图已成功生成并保存！\n")


cat("\n========== 6.3 留一法森林图 (Fig S7 - 终极修复版) ==========\n")
# 去掉 [[1]] 防止提取报错
loo_res <- mr_leaveoneout(dat_clean)
loo_dat <- as.data.frame(loo_res) %>%
  mutate(b_val = b, 
         se_val = se, 
         lo = b_val - 1.96 * se_val, 
         up = b_val + 1.96 * se_val, 
         SNP_label = SNP)

p_loo <- ggplot(loo_dat, aes(x = SNP_label, y = b_val)) +
  geom_hline(yintercept = ivw_est$b, linetype = "dashed", color = "#E64B35", linewidth = 1.2) +
  geom_errorbar(aes(ymin = lo, ymax = up), width = 0.2, color = "grey50", linewidth = 0.8) +
  geom_point(size = 3, color = "#00A087", alpha = 0.9, shape = 21, fill = "white", stroke = 1.2) +
  geom_text_repel(aes(label = SNP_label), size = 3.5, max.overlaps = 30, box.padding = 0.8, segment.color = "grey80", direction = "y") +
  labs(title = "Leave-One-Out Sensitivity Analysis", 
       subtitle = "Assessing robustness to single pleiotropic variants", 
       x = "Omitted SNP", 
       y = "IVW estimate after omission") +
  theme_classic(base_size = 13) + 
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        axis.line.x = element_blank(), 
        plot.title = element_text(face="bold", hjust=0.5), 
        plot.subtitle = element_text(hjust=0.5))

ggsave(file.path(out_dir, "Fig_S7_LeaveOneOut_Clean_Nature.pdf"), plot = p_loo, width = 14, height = 8)
cat("✅ Fig S7 留一法图已成功生成并保存！\n")


cat("\n========== 7. 补全附属文件：基础 MR 表格与原始留一法图 (S1, S2, S3) ==========\n")
write.csv(mr_res_raw, file.path(out_dir, "Table_S1_UVMR_Raw_Results.csv"), row.names = FALSE)
write.csv(dat, file.path(out_dir, "Table_S2_Harmonised_SNPs_UVMR.csv"), row.names = FALSE)
write.csv(mv_exp_dat, file.path(out_dir, "Table_S3_MVMR_Instrumental_Variables.csv"), row.names = FALSE)

loo_dat_raw <- as.data.frame(mr_leaveoneout(dat)) %>%
  mutate(b_val = b, se_val = se, lo = b_val - 1.96 * se_val, up = b_val + 1.96 * se_val, SNP_label = SNP)
ivw_est_raw <- mr_res_raw %>% filter(method == "Inverse variance weighted")

p_loo_raw <- ggplot(loo_dat_raw, aes(x = SNP_label, y = b_val)) +
  geom_hline(yintercept = ivw_est_raw$b, linetype = "dashed", color = "#E64B35", linewidth = 1.2) +
  geom_errorbar(aes(ymin = lo, ymax = up), width = 0.2, color = "grey50", linewidth = 0.8) +
  geom_point(size = 3, color = "#00A087", alpha = 0.9, shape = 21, fill = "white", stroke = 1.2) +
  labs(title = "Leave-One-Out Analysis (Raw)", 
       subtitle = "Univariable MR of TG on Grip Strength", 
       x = "Omitted SNP", 
       y = "IVW estimate after omission") +
  theme_classic(base_size = 13) + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank())

ggsave(file.path(out_dir, "Fig_S3_LeaveOneOut_Raw_Nature.pdf"), plot = p_loo_raw, width = 14, height = 8)
cat("✅ 基础 MR 附件表格与 Fig S3 原始留一法图已全部落盘！\n")


cat("\n========== 8. 补全附属文件：双向 MR (Table S5) ==========\n")
cat("正在运行双向 MR: 握力(暴露) -> TG(结局)...\n")
exp_bidir <- tryCatch({ extract_instruments("ukb-b-10215", p1 = 5e-8) }, error = function(e) { NULL })

if(!is.null(exp_bidir) && nrow(exp_bidir) > 0) {
  out_bidir <- extract_outcome_data(snps = exp_bidir$SNP, outcomes = "ieu-b-111", proxies = TRUE)
  if(!is.null(out_bidir) && nrow(out_bidir) > 0) {
    dat_bidir <- harmonise_data(exposure_dat = exp_bidir, outcome_dat = out_bidir)
    res_bidir <- mr(dat_bidir, method_list = c("mr_ivw", "mr_weighted_median", "mr_egger_regression"))
    write.csv(res_bidir, file.path(out_dir, "Table_S5_Bidirectional_MR.csv"), row.names = FALSE)
    cat("✅ 双向 MR 运行成功，已保存为 Table S5。\n")
  } else {
    cat("⚠️ TG 数据集中未匹配到足量 SNPs。\n")
  }
} else {
  cat("⚠️ 无法提取握力的显著工具变量。\n")
}

cat("\n🎉 数据工程全剧终！所有图表与表格数据完美闭环！\n")

# 补充 Fig S5：Clean UVMR 森林图
cat("\n========== 补充 Fig S5：Clean UVMR 森林图 ==========\n")
res_single_clean <- mr_singlesnp(dat_clean)
p_forest_clean <- mr_forest_plot(res_single_clean)
pdf(file.path(out_dir, "Fig_S5_Forest_Clean_Nature.pdf"), width = 10, height = 14)
print(p_forest_clean[[1]] + 
      labs(title = "Individual SNP Effects on Grip Strength (Clean Instruments)") + 
      theme_classic(base_size = 12) + 
      theme(plot.title = element_text(face = "bold", hjust = 0.5)))
dev.off()
cat("✅ Fig S5 森林图已生成。\n")

