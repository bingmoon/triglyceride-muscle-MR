# =====================================================================
# MSM 项目 - 阶段一：代谢组学（ST003803）
# 工作目录：/Users/bing/MS
# 说明：所有图表命名前缀 "Stage1_"，表示阶段一产出。
# =====================================================================

# 设置工作目录
setwd("/Users/bing/MS")

# 加载必要的包
library(jsonlite)
library(dplyr)
library(stringr)
library(impute)      # 需要先安装 BiocManager::install("impute")
library(ropls)
library(ggplot2)
library(ggpubr)
library(ggsci)
library(tidyr)

# ========== 1. 读取并清洗表达矩阵 ==========
cat("\n📥 正在读取代谢组数据 MSdata_ST003803_1.txt...\n")

file_path <- "MSdata_ST003803_1.txt"
if(!file.exists(file_path)) stop("❌ 找不到文件，请确认 MSdata_ST003803_1.txt 在 MS 文件夹下。")

raw_df <- read.delim(file_path, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
rownames(raw_df) <- make.unique(as.character(raw_df[[1]]))
expr_data <- raw_df[, -1]

# 移除可能的 RefMet_name 列
if("RefMet_name" %in% colnames(expr_data)) {
  expr_data <- expr_data[, colnames(expr_data) != "RefMet_name"]
}

# 获取分组信息（从 Metabolomics Workbench API 在线拉取）
cat("🌐 正在从 Metabolomics Workbench 获取分组信息...\n")
url_3803 <- "https://www.metabolomicsworkbench.org/rest/study/study_id/ST003803/factors"
raw_meta <- bind_rows(fromJSON(url_3803))

meta_3803 <- raw_meta %>%
  mutate(
    mb_sample_id = trimws(ifelse("sample_id" %in% colnames(raw_meta), sample_id, mb_sample_id)),
    local_sample_id = trimws(local_sample_id),
    RawGroup = str_match(factors, "Treatment:\\s*([^|\\s]+)")[,2],
    RawGroup = ifelse(is.na(RawGroup), factors, RawGroup),
    Group = case_when(
      grepl("control", RawGroup, ignore.case = TRUE) ~ "Control",
      grepl("hindlimb\\+taurine", RawGroup, ignore.case = TRUE) ~ "DMA_Taurine",
      grepl("hindlimb", RawGroup, ignore.case = TRUE) ~ "DMA",
      TRUE ~ "Unknown"
    )
  ) %>%
  dplyr::select(mb_sample_id, local_sample_id, Group)

cat("分组拉取成功，样本分布：\n")
print(table(meta_3803$Group))

# 匹配样本
sample_cols <- colnames(expr_data)
if(any(sample_cols %in% meta_3803$local_sample_id)) {
  match_col <- "local_sample_id"
} else if(any(sample_cols %in% meta_3803$mb_sample_id)) {
  match_col <- "mb_sample_id"
} else {
  stop("❌ 表达矩阵列名与云端样本 ID 无法匹配！")
}

meta_matched <- meta_3803[meta_3803[[match_col]] %in% sample_cols, ]
expr_matched <- expr_data[, meta_matched[[match_col]]]
cat(sprintf("样本匹配成功，共 %d 个样本。\n", ncol(expr_matched)))

# 转为数值矩阵，处理0值
expr_mat <- as.matrix(expr_matched)
class(expr_mat) <- "numeric"
expr_mat[expr_mat == 0] <- NA

# 质控：剔除缺失率 > 30% 的特征
missing_rate <- rowMeans(is.na(expr_mat))
expr_mat_filtered <- expr_mat[missing_rate <= 0.3, ]
cat(sprintf("缺失率过滤后保留 %d 个代谢物（原始 %d 个）。\n", 
            nrow(expr_mat_filtered), nrow(expr_mat)))

# TIC 归一化
tic <- colSums(expr_mat_filtered, na.rm = TRUE)
expr_mat_tic <- sweep(expr_mat_filtered, 2, tic, "/") * median(tic, na.rm = TRUE)

# KNN 插补
set.seed(2026)
k_val <- min(10, ncol(expr_mat_tic) - 1)
cat("正在进行 KNN 缺失值插补...\n")
imputed_res <- impute.knn(expr_mat_tic, k = k_val)$data

# Log2 转换
min_val <- min(imputed_res[imputed_res > 0], na.rm = TRUE)
imputed_log2 <- log2(imputed_res + (min_val / 2))

# 组装最终矩阵
data_clean <- as.data.frame(t(imputed_log2))
data_clean$SampleID <- rownames(data_clean)
data_clean <- data_clean %>%
  left_join(meta_matched %>% dplyr::select(all_of(match_col), Group), 
            by = c("SampleID" = match_col)) %>%
  dplyr::select(SampleID, Group, everything())

data_clean$Group <- factor(data_clean$Group, levels = c("Control", "DMA", "DMA_Taurine"))

saveRDS(data_clean, "Stage1_data_clean.rds")
cat("清洗完毕，矩阵已保存为 Stage1_data_clean.rds。\n")

# ========== 2. PCA 无监督分析 ==========
cat("\n🔬 正在进行 PCA 分析...\n")
expr_for_pca <- data_clean[, 3:ncol(data_clean)]
pca_res <- prcomp(expr_for_pca, center = TRUE, scale. = TRUE)

pca_scores <- as.data.frame(pca_res$x)
pca_scores$Group <- data_clean$Group

variance <- round(100 * pca_res$sdev^2 / sum(pca_res$sdev^2), 1)
pc1_lab <- paste0("PC1 (", variance[1], "%)")
pc2_lab <- paste0("PC2 (", variance[2], "%)")

sci_colors <- c("Control" = "#E64B35", "DMA" = "#4DBBD5", "DMA_Taurine" = "#00A087")

p_pca <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = Group, fill = Group)) +
  geom_point(size = 4, alpha = 0.85, shape = 21, stroke = 1.2, color = "white") +
  stat_ellipse(geom = "polygon", alpha = 0.1, linetype = "dashed", linewidth = 0.8) +
  scale_fill_manual(values = sci_colors) +
  scale_color_manual(values = sci_colors) +
  labs(x = pc1_lab, y = pc2_lab) +
  theme_classic(base_size = 15) +
  theme(
    axis.title = element_text(face = "bold", color = "black"),
    axis.text = element_text(color = "black", size = 12),
    axis.line = element_line(linewidth = 1, color = "black"),
    legend.position = c(0.85, 0.85),
    legend.title = element_blank(),
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.text = element_text(size = 12, face = "bold")
  )

ggsave("Stage1_Figure_PCA.pdf", plot = p_pca, width = 6.5, height = 5)
cat("PCA 图已保存为 Stage1_Figure_PCA.pdf。\n")

# ========== 3. PLS-DA 与置换检验 ==========
cat("\n🚀 正在进行 PLS-DA 分析（含 1000 次置换检验）...\n")
X <- as.matrix(data_clean[, 3:ncol(data_clean)])
Y <- data_clean$Group

plsda_model <- opls(x = X, y = Y, 
                    algoC = "nipals", predI = 2, orthoI = 0,
                    permI = 1000, fig.pdfC = "none", info.txtC = "none")

# 得分图
plsda_scores <- as.data.frame(plsda_model@scoreMN)
colnames(plsda_scores) <- c("t1", "t2")
plsda_scores$Group <- Y

p_plsda_score <- ggplot(plsda_scores, aes(x = t1, y = t2, color = Group, fill = Group)) +
  geom_point(size = 4, alpha = 0.85, shape = 21, stroke = 1.2, color = "white") +
  stat_ellipse(geom = "polygon", alpha = 0.1, linetype = "dashed", linewidth = 0.8) +
  scale_fill_manual(values = sci_colors) +
  scale_color_manual(values = sci_colors) +
  labs(title = "PLS-DA Score Plot", x = "t1", y = "t2") +
  theme_classic(base_size = 15) +
  theme(
    axis.title = element_text(face = "bold", color = "black"),
    axis.text = element_text(color = "black", size = 12),
    axis.line = element_line(linewidth = 1, color = "black"),
    legend.position = c(0.85, 0.85),
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

ggsave("Stage1_Figure_PLSDA_Score.pdf", plot = p_plsda_score, width = 6.5, height = 5)

# 置换检验图
pdf("Stage1_Figure_PLSDA_Permutation.pdf", width = 6, height = 5)
plot(plsda_model, typeVc = "permutation")
dev.off()

# VIP 值提取与棒棒糖图
vip_values <- getVipVn(plsda_model)
vip_df <- data.frame(Metabolite = names(vip_values), VIP = vip_values) %>%
  arrange(desc(VIP))

vip_core <- vip_df %>% filter(VIP > 1)
top_vip <- head(vip_core, 20)
top_vip$Metabolite <- factor(top_vip$Metabolite, levels = rev(top_vip$Metabolite))

p_vip <- ggplot(top_vip, aes(x = Metabolite, y = VIP)) +
  geom_segment(aes(xend = Metabolite, y = 0, yend = VIP), color = "grey50", linewidth = 1) +
  geom_point(size = 4, color = "#E64B35") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "blue") +
  coord_flip() +
  labs(title = "Top Metabolites by VIP Score (VIP > 1)", x = "", y = "VIP Score") +
  theme_classic(base_size = 14) +
  theme(
    axis.text.y = element_text(color = "black", face = "bold"),
    axis.title.x = element_text(face = "bold", color = "black"),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

ggsave("Stage1_Figure_VIP_Lollipop.pdf", plot = p_vip, width = 6, height = 7)
write.csv(vip_df, "Stage1_Table_VIP_Scores.csv", row.names = FALSE)
cat("PLS-DA 及 VIP 图表已保存。\n")

# ========== 4. 单变量验证（Top 5 代谢物） ==========
cat("\n🔍 正在绘制 Top 5 核心代谢物箱线图...\n")
top_targets <- c("Succinate", "Glycerol", "Isoleucine", "Carnitine", "Glutamate")

# 确保列名存在
missing <- setdiff(top_targets, colnames(data_clean))
if(length(missing) > 0) {
  for(i in seq_along(top_targets)) {
    match_col <- grep(top_targets[i], colnames(data_clean), ignore.case = TRUE, value = TRUE)
    if(length(match_col) > 0) top_targets[i] <- match_col[1]
  }
}

data_sub <- data_clean %>% dplyr::select(SampleID, Group, all_of(top_targets))
data_long <- data_sub %>%
  pivot_longer(cols = all_of(top_targets), names_to = "Metabolite", values_to = "Expression")
data_long$Metabolite <- factor(data_long$Metabolite, levels = top_targets)

my_comps <- list(c("Control","DMA"), c("DMA","DMA_Taurine"), c("Control","DMA_Taurine"))

p_box <- ggviolin(data_long, x = "Group", y = "Expression", fill = "Group",
                  palette = sci_colors, add = "boxplot", 
                  add.params = list(fill = "white", width = 0.1), trim = FALSE) +
  geom_jitter(width = 0.1, size = 1.5, alpha = 0.6, color = "black") +
  facet_wrap(~ Metabolite, scales = "free_y", nrow = 1) +
  stat_compare_means(comparisons = my_comps, method = "wilcox.test", 
                     label = "p.signif", tip.length = 0.02) +
  labs(y = "Relative Abundance (Log2 TIC)", x = "") +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", color = "black"),
    axis.title.y = element_text(face = "bold", color = "black"),
    strip.text = element_text(face = "bold", size = 12, color = "black"),
    strip.background = element_rect(fill = "grey90", color = "black"),
    legend.position = "none",
    panel.grid.major.x = element_blank()
  )

ggsave("Stage1_Figure_Univariate_Boxplots.pdf", plot = p_box, width = 12, height = 5)
cat("单变量箱线图已保存。\n")

cat("\n🎉 阶段一全部完成！所有图表和表格均以 'Stage1_' 开头命名。\n")
