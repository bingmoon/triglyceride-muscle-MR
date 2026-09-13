# ---
# Archived from: 4.MSM（ST003803） 1.md  (L1967-2146, block 9)
# Body md5: 8cf2f6253500a195a01d7511a20a82cc   | 3-line header added by release tooling
# ---
# =====================================================================
# 【块 8·新增 2026-09-13】绝对握力 + BMI 加权回归（带诊断版）
# ---------------------------------------------------------------------
# 目的：为稿件数字 β=0.368, P=0.400 建立可溯源产物（此前该数字仅来自口述，磁盘无输出）
# 运行：本地 R，NHANES G/H 两周期合并，svyglm + Taylor series linearization
# 产出：/Users/bing/MS/results/
#         weighted_abs_grip_BMI.rds
#         weighted_abs_grip_BMI.csv
#
# 【口径声明】
#   样本：NHANES 2011-2012 (G) + 2013-2014 (H) 合并，n = 4180
#         （与 MS 主数据集 NHANES_Adult_Grip_Main.rds 的 4180×11 一致 ✅）
#   权重：WTMEC2YR / 2（两周期合并，2 年权重折半，标准做法）
#   设计：ids = SDMVPSU, strata = SDMVSTRA, nest = TRUE
#   模型：Grip(kg, MGDCGSZ 绝对握力) ~ Log2_TG + Age + Gender + Cycle + BMI
# =====================================================================

suppressMessages({
    library(nhanesA)
    library(dplyr)
    library(survey)
    library(broom)
})

setwd("/Users/bing/MS")
dir.create("results", showWarnings = FALSE)

# ========== 1. 提取数据 ==========
cat("\n[1/4] 提取 NHANES G/H 数据...\n")

fetch_cycle <- function(cy) {
    demo <- nhanes(paste0("DEMO_", cy))
    bmx  <- nhanes(paste0("BMX_", cy))
    trig <- nhanes(paste0("TRIGLY_", cy))
    fast <- nhanes(paste0("FASTQX_", cy))
    grip <- nhanes(paste0("MGX_", cy))

    for (d in list(demo, bmx, trig, fast, grip)) {
        d$SEQN <- as.numeric(as.character(d$SEQN))
    }

    df <- demo %>%
        transmute(
            SEQN,
            Age = as.numeric(RIDAGEYR),
            Gender_raw = as.character(RIAGENDR),
            SDMVPSU = as.numeric(SDMVPSU),
            SDMVSTRA = as.numeric(SDMVSTRA),
            WTMEC2YR = as.numeric(WTMEC2YR)
        ) %>%
        left_join(bmx  %>% transmute(SEQN, BMI = as.numeric(BMXBMI),
                                     Weight_kg = as.numeric(BMXWT)), by = "SEQN") %>%
        left_join(trig %>% transmute(SEQN, TG = as.numeric(LBXTR)), by = "SEQN") %>%
        left_join(fast %>% transmute(SEQN, FastHr = as.numeric(PHAFSTHR)), by = "SEQN") %>%
        left_join(grip %>% transmute(SEQN, Grip = as.numeric(MGDCGSZ)), by = "SEQN") %>%
        mutate(Cycle = cy)

    return(df)
}

df_G <- fetch_cycle("G")
df_H <- fetch_cycle("H")
df_all <- bind_rows(df_G, df_H)

# ========== 诊断 1：原始 Gender_raw 的取值 ==========
cat("\n========== 诊断 1: Gender_raw 取值 ==========\n")
print(table(df_all$Gender_raw, useNA = "ifany"))
# 实测：Female 10072 / Male 9859（无 NA）

# 更鲁棒的性别转换
df_all <- df_all %>%
    mutate(
        Gender = case_when(
            Gender_raw %in% c("1", "Male", "MALE", "male") ~ "Male",
            Gender_raw %in% c("2", "Female", "FEMALE", "female") ~ "Female",
            TRUE ~ NA_character_
        )
    )

# ========== 2. 过滤 ==========
cat("\n[2/4] 过滤样本...\n")

df_all <- df_all %>%
    filter(
        Age >= 20,
        FastHr >= 8,
        TG < 500,
        !is.na(Grip), !is.na(TG), !is.na(BMI),
        !is.na(Gender), !is.na(WTMEC2YR),
        !is.na(SDMVPSU), !is.na(SDMVSTRA)
    ) %>%
    mutate(
        Log2_TG = log2(TG),
        Gender = factor(Gender, levels = c("Male", "Female")),
        Cycle = factor(Cycle, levels = c("G", "H")),
        WTMEC2YR_adj = WTMEC2YR / 2
    )

cat(sprintf("合并后样本量：%d\n", nrow(df_all)))
# 实测：4180

# ========== 诊断 2-4：检查 Gender / Cycle 水平与 NA ==========
cat("\n========== 诊断 2: Gender 分布 ==========\n")
print(table(df_all$Gender, useNA = "ifany"))
# 实测：Male 2063 / Female 2117

cat("\n========== 诊断 3: Cycle 分布 ==========\n")
print(table(df_all$Cycle, useNA = "ifany"))
# 实测：G 2017 / H 2163

cat("\n========== 诊断 4: 各变量 NA 计数 ==========\n")
print(colSums(is.na(df_all[, c("Age", "Gender", "BMI", "Log2_TG", "Grip",
                               "SDMVPSU", "SDMVSTRA", "WTMEC2YR_adj")])))
# 实测：全部为 0

# ========== 3. 若 Gender/Cycle 只有 1 个水平则提前退出 ==========
if (length(unique(na.omit(df_all$Gender))) < 2) {
    stop("Gender 只有一个水平，无法跑回归。请检查上方诊断输出。")
}
if (length(unique(na.omit(df_all$Cycle))) < 2) {
    stop("Cycle 只有一个水平，无法跑回归。请检查上方诊断输出。")
}

# ========== 4. 构建 survey design ==========
cat("\n[3/4] 构建 survey design...\n")

grip_design <- svydesign(
    ids = ~SDMVPSU,
    strata = ~SDMVSTRA,
    weights = ~WTMEC2YR_adj,
    data = df_all,
    nest = TRUE
)

# ========== 5. 加权回归 ==========
cat("\n[4/4] 加权回归: Grip ~ Log2_TG + Age + Gender + Cycle + BMI\n")

model_abs_bmi <- svyglm(
    Grip ~ Log2_TG + Age + Gender + Cycle + BMI,
    design = grip_design
)

cat("\n========== 结果 ==========\n")
print(tidy(model_abs_bmi, conf.int = TRUE))

tg_coef <- tidy(model_abs_bmi) %>% filter(term == "Log2_TG")
cat(sprintf("\n>>> TG 效应: β = %.4f, SE = %.4f, P = %.4f\n",
            tg_coef$estimate, tg_coef$std.error, tg_coef$p.value))

# ========== 6. 落盘 ==========
saveRDS(model_abs_bmi, "results/weighted_abs_grip_BMI.rds")
write.csv(tidy(model_abs_bmi, conf.int = TRUE),
          "results/weighted_abs_grip_BMI.csv", row.names = FALSE)

cat("\n✅ 落盘: results/weighted_abs_grip_BMI.rds / .csv\n")

# ---- 实测真值（磁盘产物逐位核对，2026-09-13）----
# >>> TG 效应: β = 0.3679, SE = 0.4300, P = 0.3998
#     95% CI = [-0.5144, 1.2502]
#     模型全系数：
#       (Intercept)   94.2739  SE 2.8603   P = 2.41e-23
#       Log2_TG        0.3679  SE 0.4300   P = 4.00e-01   ← 目标量（稿件 β=0.368, P=0.400 ✅）
#       Age           -0.3801  SE 0.0164   P = 2.24e-19
#       GenderFemale -33.0977  SE 0.6020   P = 2.95e-29
#       CycleH        -0.1995  SE 0.7289   P = 7.86e-01
#       BMI            0.3707  SE 0.0466   P = 1.52e-08
#   n = 4180；权重 WTMEC2YR/2；Taylor series linearization
#
# 【方法学要点 —— 供 Limitations / 讨论引用】
#   1) BMI 是 TG 与握力之间的潜在中介/共同效应，纳入线性 BMI 属"过度调整"，
#      可能引入 collider/overadjustment bias。未加权、不含 BMI 的同一模型得
#      β=+0.743, P=0.0075（阳性）；含 BMI 后 β=0.368, P=0.400（阴性）
#      → "non-significant after BMI adjustment" 的结论完全由 BMI 调整驱动，必须如实披露。
#   2) TG 测于空腹子样本，NHANES 分析规范更宜用空腹子样本权重（WTSAF2YR）；
#      本模型用 MEC 体检权重（WTMEC2YR/2），属可辩护但需声明的选择。
#   3) 未调整种族/民族、吸烟、饮酒、体力活动 → 残余混杂。
#   4) 效应量口径：Log2_TG 的 β = TG 每翻一倍时握力变化（kg）；CI [-0.51, +1.25] 较宽，
#      对中等效应把握度不足（SE=0.43）。
