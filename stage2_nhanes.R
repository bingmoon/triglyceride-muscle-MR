# =====================================================================
# 项目名称：MSM (Muscle, Strength & Metabolism)
# 分析阶段：二、NHANES 多周期临床流行病学分析 (成年人全年龄段 Age >= 20)
# 核心目标：验证脂毒性 (TG) 与相对握力 / 相对 DXA 的独立负相关，揭示 BMI 悖论
# =====================================================================

setwd("/Users/bing/MSM")

# 确保必要的包都已安装并加载
suppressMessages({
  if(!require(nhanesA)) install.packages("nhanesA", ask = FALSE)
  if(!require(dplyr)) install.packages("dplyr", ask = FALSE)
  if(!require(ggplot2)) install.packages("ggplot2", ask = FALSE)
  if(!require(ggpubr)) install.packages("ggpubr", ask = FALSE)
  if(!require(broom)) install.packages("broom", ask = FALSE)
  if(!require(tibble)) install.packages("tibble", ask = FALSE)
})

library(nhanesA)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(broom)
library(tibble)

# ========== 1. 拉取所有周期的原始数据 ==========
cat("\n📥 [1/6] 正在拉取 2011-2018 四个周期的原始数据 (请保持网络畅通)...\n")

# G周期 (2011-2012)
demo_G  <- tryCatch(nhanes("DEMO_G"), error=function(e) NULL); bmx_G  <- tryCatch(nhanes("BMX_G"), error=function(e) NULL); trig_G  <- tryCatch(nhanes("TRIGLY_G"), error=function(e) NULL)
fast_G  <- tryCatch(nhanes("FASTQX_G"), error=function(e) NULL); grip_G <- tryCatch(nhanes("MGX_G"), error=function(e) NULL); dxa_G   <- tryCatch(nhanes("DXX_G"), error=function(e) NULL)

# H周期 (2013-2014)
demo_H  <- tryCatch(nhanes("DEMO_H"), error=function(e) NULL); bmx_H  <- tryCatch(nhanes("BMX_H"), error=function(e) NULL); trig_H  <- tryCatch(nhanes("TRIGLY_H"), error=function(e) NULL)
fast_H  <- tryCatch(nhanes("FASTQX_H"), error=function(e) NULL); grip_H <- tryCatch(nhanes("MGX_H"), error=function(e) NULL); dxa_H   <- tryCatch(nhanes("DXX_H"), error=function(e) NULL)

# I周期 (2015-2016, 无握力)
demo_I  <- tryCatch(nhanes("DEMO_I"), error=function(e) NULL); bmx_I  <- tryCatch(nhanes("BMX_I"), error=function(e) NULL); trig_I  <- tryCatch(nhanes("TRIGLY_I"), error=function(e) NULL)
fast_I  <- tryCatch(nhanes("FASTQX_I"), error=function(e) NULL); dxa_I  <- tryCatch(nhanes("DXX_I"), error=function(e) NULL)

# J周期 (2017-2018, 无握力)
demo_J  <- tryCatch(nhanes("DEMO_J"), error=function(e) NULL); bmx_J  <- tryCatch(nhanes("BMX_J"), error=function(e) NULL); trig_J  <- tryCatch(nhanes("TRIGLY_J"), error=function(e) NULL)
fast_J  <- tryCatch(nhanes("FASTQX_J"), error=function(e) NULL); dxa_J  <- tryCatch(nhanes("DXX_J"), error=function(e) NULL)


# ========== 2. 构建握力主队列 (G+H 周期, Age >= 20) ==========
cat("\n⚙️ [2/6] 正在构建握力功能队列...\n")

build_grip_cycle <- function(demo, bmx, trig, fast, grip, cycle_name) {
  if(is.null(demo) | is.null(grip)) return(data.frame()) 
  
  # 强制 SEQN 为数值型
  demo <- demo %>% mutate(SEQN = as.numeric(as.character(SEQN)))
  bmx <- bmx %>% mutate(SEQN = as.numeric(as.character(SEQN)))
  trig <- trig %>% mutate(SEQN = as.numeric(as.character(SEQN)))
  fast <- fast %>% mutate(SEQN = as.numeric(as.character(SEQN)))
  grip <- grip %>% mutate(SEQN = as.numeric(as.character(SEQN)))
  
  df <- demo %>%
    dplyr::select(SEQN, RIDAGEYR, RIAGENDR) %>%
    rename(Age = RIDAGEYR, Gender = RIAGENDR) %>%
    # 极度鲁棒的性别转换
    mutate(
      Gender = as.character(Gender),
      Gender = case_when(
        Gender %in% c("1", "Male", "MALE", "male") ~ "Male",
        Gender %in% c("2", "Female", "FEMALE", "female") ~ "Female",
        TRUE ~ NA_character_
      ),
      Gender = factor(Gender, levels = c("Male", "Female"))
    ) %>%
    left_join(bmx %>% dplyr::select(SEQN, BMXBMI, BMXWT), by = "SEQN") %>%
    rename(BMI = BMXBMI, Weight_kg = BMXWT) %>%
    left_join(trig %>% dplyr::select(SEQN, LBXTR), by = "SEQN") %>%
    rename(Triglycerides = LBXTR) %>%
    left_join(fast %>% dplyr::select(SEQN, PHAFSTHR), by = "SEQN") %>%
    left_join(grip %>% dplyr::select(SEQN, MGDCGSZ), by = "SEQN") %>%
    rename(Grip_Strength = MGDCGSZ) %>%
    # 核心修改：Age >= 20 拓展到全体成年人
    filter(Age >= 20, PHAFSTHR >= 8, Triglycerides < 500) %>%
    mutate(
      Log2_TG = log2(Triglycerides),
      Relative_Grip = Grip_Strength / Weight_kg,
      Cycle = cycle_name
    ) %>%
    filter(!is.na(Relative_Grip), !is.na(Triglycerides), !is.na(BMI), !is.na(Gender))
  
  return(df)
}

df_grip_main <- bind_rows(
  build_grip_cycle(demo_G, bmx_G, trig_G, fast_G, grip_G, "G"),
  build_grip_cycle(demo_H, bmx_H, trig_H, fast_H, grip_H, "H")
)
cat(sprintf("✅ 握力队列构建成功，包含 %d 名成年人 (Age 20-80+)\n", nrow(df_grip_main)))


# ========== 3. 构建 DXA 结构队列 (G+H+I+J 周期, Age >= 20) ==========
cat("\n⚙️ [3/6] 正在构建 DXA 结构队列...\n")

build_dxa_cycle <- function(demo, bmx, trig, fast, dxa, cycle_name) {
  if(is.null(demo) | is.null(dxa)) return(data.frame()) 
  
  # 强制 SEQN 为数值型
  demo <- demo %>% mutate(SEQN = as.numeric(as.character(SEQN)))
  bmx <- bmx %>% mutate(SEQN = as.numeric(as.character(SEQN)))
  trig <- trig %>% mutate(SEQN = as.numeric(as.character(SEQN)))
  fast <- fast %>% mutate(SEQN = as.numeric(as.character(SEQN)))
  dxa <- dxa %>% mutate(SEQN = as.numeric(as.character(SEQN)))
  
  df <- demo %>%
    dplyr::select(SEQN, RIDAGEYR, RIAGENDR) %>%
    rename(Age = RIDAGEYR, Gender = RIAGENDR) %>%
    mutate(
      Gender = as.character(Gender),
      Gender = case_when(
        Gender %in% c("1", "Male", "MALE", "male") ~ "Male",
        Gender %in% c("2", "Female", "FEMALE", "female") ~ "Female",
        TRUE ~ NA_character_
      ),
      Gender = factor(Gender, levels = c("Male", "Female"))
    ) %>%
    left_join(bmx %>% dplyr::select(SEQN, BMXBMI, BMXWT), by = "SEQN") %>%
    rename(BMI = BMXBMI, Weight_kg = BMXWT) %>%
    left_join(trig %>% dplyr::select(SEQN, LBXTR), by = "SEQN") %>%
    rename(Triglycerides = LBXTR) %>%
    left_join(fast %>% dplyr::select(SEQN, PHAFSTHR), by = "SEQN") %>%
    left_join(dxa %>% dplyr::select(SEQN, any_of(c("DXDLALE", "DXDRALE", "DXDLLLE", "DXDRLLE"))), by = "SEQN") %>%
    # 核心修改：Age >= 20 (由于 NHANES DXA 限制，实际获取到的是 20-59 岁人群)
    filter(Age >= 20, PHAFSTHR >= 8, Triglycerides < 500) %>%
    mutate(
      Log2_TG = log2(Triglycerides),
      ASM_kg = ifelse("DXDLALE" %in% names(.) & !is.na(DXDLALE) & !is.na(DXDRALE) & !is.na(DXDLLLE) & !is.na(DXDRLLE),
                       (DXDLALE + DXDRALE + DXDLLLE + DXDRLLE) / 1000, NA_real_),
      Relative_ASM = ASM_kg / Weight_kg,
      Cycle = cycle_name
    ) %>%
    filter(!is.na(Relative_ASM), !is.na(Gender))
  
  return(df)
}

df_dxa_all <- bind_rows(
  build_dxa_cycle(demo_G, bmx_G, trig_G, fast_G, dxa_G, "G"),
  build_dxa_cycle(demo_H, bmx_H, trig_H, fast_H, dxa_H, "H"),
  build_dxa_cycle(demo_I, bmx_I, trig_I, fast_I, dxa_I, "I"),
  build_dxa_cycle(demo_J, bmx_J, trig_J, fast_J, dxa_J, "J")
)
cat(sprintf("✅ DXA 队列构建成功，包含 %d 名中青年 (Age 20-59)\n", nrow(df_dxa_all)))


# ========== 4. 智能自适应回归模型验证 ==========
cat("\n📊 [4/6] 正在执行多因素线性回归分析...\n")

# [A] 握力回归 (相对握力 ~ TG + 年龄 + 性别 + 周期)
has_cycles_grip <- length(unique(df_grip_main$Cycle)) > 1
form_grip <- "Relative_Grip ~ Log2_TG + Age + Gender"
if(has_cycles_grip) form_grip <- paste(form_grip, "+ Cycle")

model_grip <- lm(as.formula(form_grip), data = df_grip_main)
cat(">>> 握力功能学模型 (Age >= 20):\n")
print(summary(model_grip)$coefficients)

# [B] DXA 回归 (相对 ASM ~ TG + 年龄 + 性别 + 周期)
has_cycles_dxa <- length(unique(df_dxa_all$Cycle)) > 1
form_dxa <- "Relative_ASM ~ Log2_TG + Age + Gender"
if(has_cycles_dxa) form_dxa <- paste(form_dxa, "+ Cycle")

model_dxa <- lm(as.formula(form_dxa), data = df_dxa_all)
cat("\n>>> DXA 结构学模型 (Age 20-59):\n")
print(summary(model_dxa)$coefficients)


# ========== 5. 高质量图表生成 ==========
cat("\n🎨 [5/6] 正在生成顶级期刊格式的可视化图表...\n")

theme_nature <- theme_classic(base_size=12) +
  theme(axis.title=element_text(face="bold",size=13),
        axis.text=element_text(color="black",size=11),
        axis.line=element_line(linewidth=0.8,color="black"),
        panel.grid=element_blank(),
        plot.title=element_text(face="bold",size=14,hjust=0.5),
        plot.subtitle=element_text(hjust=0.5,size=11,color="grey30"),
        legend.position="top",legend.title=element_blank(),
        legend.text=element_text(size=11),
        strip.text=element_text(face="bold",size=12,color="white"),
        strip.background=element_rect(fill="#333333",color="black"))

nature_colors <- c("G"="#E64B35","H"="#4DBBD5","I"="#00A087","J"="#3C5488")

# 图 1：握力散点图 (功能学)
p1 <- ggplot(df_grip_main, aes(x=Log2_TG, y=Relative_Grip, color=Cycle)) +
  geom_point(alpha=0.5, size=2) + 
  geom_smooth(method="lm", aes(group=1), color="black", linewidth=1.2, fill="grey70") +
  stat_cor(method="spearman", label.x.npc="left", label.y.npc="bottom", size=4.5, fontface="bold") +
  scale_color_manual(values=nature_colors) +
  labs(title="Fasting TG and Relative Grip Strength in Adults", 
       subtitle=sprintf("NHANES 2011–2014 (Age ≥ 20, N = %d)", nrow(df_grip_main)),
       x=expression(bold("Log"[2]*" Triglycerides")), y="Relative Grip Strength (kg/kg)") + theme_nature
ggsave("Figure_2.6_Adult_RelativeGrip_TG.pdf", p1, width=7.5, height=6, dpi=300)

# 图 2：DXA 分面散点图 (结构学)
df_dxa_plot <- df_dxa_all %>%
  mutate(Cohort = ifelse(Cycle %in% c("G","H"), "Discovery (2011–2014)", "Replication (2015–2018)"))
p2 <- ggplot(df_dxa_plot, aes(x=Log2_TG, y=Relative_ASM, color=Cycle)) +
  geom_point(alpha=0.5, size=2) + 
  geom_smooth(method="lm", aes(group=1), color="black", linewidth=1.2, fill="grey70") +
  facet_wrap(~Cohort, scales="free_x", nrow=1) +
  stat_cor(method="spearman", label.x.npc="left", label.y.npc="bottom", size=4, fontface="bold") +
  scale_color_manual(values=nature_colors) +
  labs(title="Fasting TG and DXA-Derived Relative ASM", 
       subtitle=sprintf("NHANES 2011–2018 (Age 20-59, N = %d)", nrow(df_dxa_all)),
       x=expression(bold("Log"[2]*" Triglycerides")), y="Relative ASM (kg/kg)") + theme_nature
ggsave("Figure_2.8_Adult_RelativeASM_TG.pdf", p2, width=12, height=5.5, dpi=300)

# 图 3：森林图 (多因素校正独立性验证)
coef_df <- broom::tidy(model_grip, conf.int = TRUE) %>% 
  filter(!term %in% c("(Intercept)", "CycleH")) %>%
  mutate(term_label = case_when(
    term=="Log2_TG" ~ "Triglycerides (Log2)",
    term=="Age" ~ "Age (per year)",
    term=="GenderFemale" ~ "Female vs Male",
    TRUE ~ term))

p3 <- ggplot(coef_df, aes(x=estimate, y=reorder(term_label, estimate))) +
  geom_vline(xintercept=0, linetype="dashed", color="grey50", linewidth=0.8) +
  geom_errorbarh(aes(xmin=conf.low, xmax=conf.high), height=0.2, linewidth=1.2, color="#2C3E50") +
  geom_point(size=4.5, shape=21, fill="#3498DB", color="black", stroke=0.8) +
  geom_text(aes(label=sprintf("β = %.3f   P = %.3f", estimate, p.value)),
            x = max(coef_df$conf.high) + max(abs(coef_df$conf.high))*0.2, hjust=0, size=4, fontface="italic") +
  labs(title="Independent Predictors of Relative Muscle Quality",
       subtitle="Multivariable Linear Regression in Adults (Age ≥ 20)",x="β (95% CI)",y="") +
  theme_nature + theme(axis.text.y=element_text(face="bold",size=12)) +
  coord_cartesian(xlim = c(min(coef_df$conf.low)-0.05, max(coef_df$conf.high)+0.15))
ggsave("Figure_2.9_Regression_ForestPlot.pdf", p3, width=8.5, height=3.5, dpi=300)

# 图 4：解密 BMI 悖论
p4a <- ggplot(df_grip_main, aes(x=BMI, y=Log2_TG)) +
  geom_point(alpha=0.4, size=1.8, color="#3C5488") + geom_smooth(method="lm",color="black",linewidth=1,fill="grey70") +
  stat_cor(method="spearman",label.x.npc="left",label.y.npc="top",size=4,fontface="bold") +
  labs(title="A: Higher BMI is Associated with Elevated TG",x="BMI (kg/m²)",y=expression(bold("Log"[2]*" Triglycerides"))) + theme_nature
p4b <- ggplot(df_grip_main, aes(x=BMI, y=Grip_Strength)) +
  geom_point(alpha=0.4, size=1.8, color="#E64B35") + geom_smooth(method="lm",color="black",linewidth=1,fill="grey70") +
  stat_cor(method="spearman",label.x.npc="left",label.y.npc="top",size=4,fontface="bold") +
  labs(title="B: Higher BMI Confers Absolute Strength Advantage",x="BMI (kg/m²)",y="Absolute Grip Strength (kg)") + theme_nature
p4 <- ggarrange(p4a, p4b, ncol=2)
ggsave("Figure_S2_BMI_Paradox_Source.pdf", p4, width=11, height=5, dpi=300)


# ========== 6. 数据及汇总保存 ==========
cat("\n💾 [6/6] 正在保存最终队列数据和回归结果...\n")

write.csv(broom::tidy(model_grip, conf.int=TRUE), "Table_2.9_Adult_Grip_Regression.csv", row.names=FALSE)
write.csv(broom::tidy(model_dxa, conf.int=TRUE), "Table_2.8_Adult_DXA_Regression.csv", row.names=FALSE)
saveRDS(df_grip_main, "NHANES_Adult_Grip_Main.rds")
saveRDS(df_dxa_all, "NHANES_Adult_DXA_All.rds")

cat("\n🎉 [完美收官] 阶段二全部代码执行完毕！所有图表已重新生成。\n")

library(dplyr)
library(broom)

# 汇总容器
summary_lines <- c()

add_line <- function(...) {
  summary_lines <<- c(summary_lines, paste(...))
}

add_line("==============================================")
add_line("阶段二：NHANES 多周期分析结果汇总")
add_line("==============================================")
add_line("")

# 1. 握力队列
if (exists("df_grip_main")) {
  add_line(">>> 握力功能队列 (Age >= 20, G+H)")
  add_line("样本量: ", nrow(df_grip_main))
  add_line("性别分布:")
  add_line(paste(capture.output(table(df_grip_main$Gender)), collapse="\n"))
  add_line("周期分布:")
  add_line(paste(capture.output(table(df_grip_main$Cycle)), collapse="\n"))
  add_line("年龄范围: ", paste(range(df_grip_main$Age), collapse=" - "))
  add_line("TG (mg/dL) 范围: ", paste(range(df_grip_main$Triglycerides), collapse=" - "))
  add_line("")
}

# 2. 握力回归
if (exists("model_grip")) {
  add_line(">>> 握力回归模型 (Relative_Grip ~ Log2_TG + Age + Gender + Cycle)")
  sum_grip <- summary(model_grip)
  coef_grip <- tidy(model_grip, conf.int = TRUE)
  add_line(paste(capture.output(print(coef_grip, row.names = FALSE)), collapse="\n"))
  add_line(sprintf("R-squared: %.4f   Adj R-squared: %.4f", sum_grip$r.squared, sum_grip$adj.r.squared))
  tg <- coef_grip %>% filter(term == "Log2_TG")
  if (nrow(tg) > 0) {
    add_line(sprintf("核心效应 Log2_TG：β = %.5f (95%% CI: %.5f, %.5f), P = %.3e", 
                     tg$estimate, tg$conf.low, tg$conf.high, tg$p.value))
  }
  add_line("")
}

# 3. DXA 队列
if (exists("df_dxa_all")) {
  add_line(">>> DXA 结构队列 (Age 20-59, 四周期)")
  add_line("样本量: ", nrow(df_dxa_all))
  add_line("性别分布:")
  add_line(paste(capture.output(table(df_dxa_all$Gender)), collapse="\n"))
  add_line("周期分布:")
  add_line(paste(capture.output(table(df_dxa_all$Cycle)), collapse="\n"))
  add_line("年龄范围: ", paste(range(df_dxa_all$Age), collapse=" - "))
  add_line("TG (mg/dL) 范围: ", paste(range(df_dxa_all$Triglycerides), collapse=" - "))
  add_line("")
}

# 4. DXA 回归
if (exists("model_dxa")) {
  add_line(">>> DXA 回归模型 (Relative_ASM ~ Log2_TG + Age + Gender + Cycle)")
  sum_dxa <- summary(model_dxa)
  coef_dxa <- tidy(model_dxa, conf.int = TRUE)
  add_line(paste(capture.output(print(coef_dxa, row.names = FALSE)), collapse="\n"))
  add_line(sprintf("R-squared: %.4f   Adj R-squared: %.4f", sum_dxa$r.squared, sum_dxa$adj.r.squared))
  tg <- coef_dxa %>% filter(term == "Log2_TG")
  if (nrow(tg) > 0) {
    add_line(sprintf("核心效应 Log2_TG：β = %.5f (95%% CI: %.5f, %.5f), P = %.3e", 
                     tg$estimate, tg$conf.low, tg$conf.high, tg$p.value))
  }
  add_line("")
}

add_line("==============================================")
add_line("结果汇总完毕")
add_line("==============================================")

# 统一输出到屏幕
cat(paste(summary_lines, collapse = "\n"), "\n")

# 保存为 txt 文件
writeLines(summary_lines, "Stage2_Results_Summary.txt")
message("✅ 结果已保存至：Stage2_Results_Summary.txt")

# =====================================================================
# 补充生成：Figure S1 - 绝对握力与 TG 的欺骗性正相关
# =====================================================================
setwd("/Users/bing/MS")
library(ggplot2)
library(ggpubr)

theme_nature <- theme_classic(base_size = 12) +
  theme(
    axis.title = element_text(face = "bold", size = 13),
    axis.text = element_text(color = "black", size = 11),
    axis.line = element_line(linewidth = 0.8, color = "black"),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, size = 11, color = "grey30"),
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 11)
  )

nature_colors <- c("G" = "#E64B35", "H" = "#4DBBD5")

p_s1 <- ggplot(df_grip_main, aes(x = Log2_TG, y = Grip_Strength, color = Cycle)) +
  geom_point(alpha = 0.5, size = 2) +
  geom_smooth(method = "lm", aes(group = 1), color = "black", linewidth = 1.2, fill = "grey70") +
  stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top", size = 4.5, fontface = "bold") +
  scale_color_manual(values = nature_colors) +
  labs(
    title = "Deceptive Correlation: Absolute Grip Strength vs. TG",
    subtitle = paste("NHANES 2011–2014, Age ≥ 20, N =", nrow(df_grip_main)),
    x = expression(bold("Log"[2]*" Triglycerides")),
    y = "Absolute Combined Grip Strength (kg)"
  ) +
  theme_nature

ggsave("Figure_S1_Absolute_Grip.pdf", plot = p_s1, width = 7.5, height = 6, dpi = 300)
cat("✅ Figure_S1_Absolute_Grip.pdf 已保存。\n")

# =====================================================================
# 阶段二补充：体脂率校正敏感性分析
# =====================================================================
cat("\n========== 体脂率校正敏感性分析 ==========\n")

# ===== 新增：20-59 岁亚组基础模型 =====
df_grip_20_59 <- df_grip_main %>% filter(Age <= 59)
m_grip_20_59 <- lm(Relative_Grip ~ Log2_TG + Age + Gender + Cycle, data = df_grip_20_59)
cat("20-59岁握力基础模型 (Log2_TG):\n")
print(summary(m_grip_20_59)$coefficients["Log2_TG", ])
cat(sprintf("20-59岁亚组样本量: %d\n", nrow(df_grip_20_59)))

library(nhanesA)

# 读取 DXA 全身脂肪量 (DXDTOFAT) 并计算体脂率
add_fat_data <- function(df, cycles) {
  df_fat <- data.frame()
  for (cy in cycles) {
    dxa <- nhanes(paste0("DXX_", cy))
    dxa$SEQN <- as.numeric(dxa$SEQN)
    if ("DXDTOFAT" %in% names(dxa)) {
      dxa_sub <- dxa %>% select(SEQN, Fat_mass = DXDTOFAT)
    } else if ("DXXTRFAT" %in% names(dxa)) {
      dxa_sub <- dxa %>% select(SEQN, Fat_mass = DXXTRFAT)
    } else {
      dxa_sub <- dxa %>% select(SEQN) %>% mutate(Fat_mass = NA)
    }
    dxa_sub$Cycle <- cy
    df_fat <- bind_rows(df_fat, dxa_sub)
  }
  df_out <- df %>% left_join(df_fat, by = c("SEQN", "Cycle"))
  df_out <- df_out %>% mutate(Fat_pct = Fat_mass / Weight_kg * 100)
  return(df_out)
}

# 握力队列体脂率校正模型
if (exists("df_grip_main")) {
  df_grip_fat <- add_fat_data(df_grip_main, c("G", "H"))
  model_grip_fat <- lm(Grip_Strength ~ Log2_TG + Age + Gender + Cycle + BMI + Fat_pct,
                       data = df_grip_fat)
  cat("\n握力模型（绝对握力，+ BMI + 体脂率）：\n")
  tg <- tidy(model_grip_fat) %>% filter(term == "Log2_TG")
  print(tg)
} else {
  cat("⚠️ df_grip_main 不存在，请先运行阶段二主代码。\n")
}

# DXA队列体脂率校正模型
if (exists("df_dxa_all")) {
  df_dxa_fat <- add_fat_data(df_dxa_all, c("G","H","I","J"))
  model_dxa_fat <- lm(ASM_kg ~ Log2_TG + Age + Gender + Cycle + BMI + Fat_pct,
                      data = df_dxa_fat)
  cat("\nDXA模型（绝对ASM，+ BMI + 体脂率）：\n")
  tg_dxa <- tidy(model_dxa_fat) %>% filter(term == "Log2_TG")
  print(tg_dxa)
} else {
  cat("⚠️ df_dxa_all 不存在，请先运行阶段二主代码。\n")
}

cat("✅ 体脂率敏感性分析完成。\n")
