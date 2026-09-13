# =====================================================================
# WTSAF2YR 空腹权重敏感性分析（真实 N 版）
# 修正：run_model() 加入 nobs(fit) —— 用模型实际纳入量取代 design 样本量
# =====================================================================

suppressMessages({
  library(nhanesA)
  library(dplyr)
  library(survey)
  library(broom)
  library(tibble)
})

setwd("/Users/bing/MS")
dir.create("results", showWarnings = FALSE)

fetch_cycle <- function(cy) {
  demo <- nhanes(paste0("DEMO_", cy))
  bmx  <- nhanes(paste0("BMX_", cy))
  trig <- nhanes(paste0("TRIGLY_", cy))
  fast <- nhanes(paste0("FASTQX_", cy))
  grip <- nhanes(paste0("MGX_", cy))
  dxa  <- nhanes(paste0("DXX_", cy))

  for (d in list(demo, bmx, trig, fast, grip, dxa)) {
    if (!is.null(d)) d$SEQN <- as.numeric(as.character(d$SEQN))
  }

  # WTSAF2YR 从 TRIGLY 模块取
  wt_saf <- trig %>% transmute(SEQN, WTSAF2YR = as.numeric(WTSAF2YR))

  df <- demo %>%
    transmute(
      SEQN,
      Age = as.numeric(RIDAGEYR),
      Gender_raw = as.character(RIAGENDR),
      SDMVPSU = as.numeric(SDMVPSU),
      SDMVSTRA = as.numeric(SDMVSTRA),
      WTMEC2YR = as.numeric(WTMEC2YR)
    ) %>%
    left_join(wt_saf, by = "SEQN") %>%
    left_join(bmx  %>% transmute(SEQN, BMI = as.numeric(BMXBMI),
                                 Weight_kg = as.numeric(BMXWT)), by = "SEQN") %>%
    left_join(trig %>% transmute(SEQN, TG = as.numeric(LBXTR)), by = "SEQN") %>%
    left_join(fast %>% transmute(SEQN, FastHr = as.numeric(PHAFSTHR)), by = "SEQN") %>%
    left_join(grip %>% transmute(SEQN, Grip = as.numeric(MGDCGSZ)), by = "SEQN") %>%
    left_join(dxa  %>% transmute(
      SEQN,
      ASM_kg = (as.numeric(DXDLALE) + as.numeric(DXDRALE) +
                    as.numeric(DXDLLLE) + as.numeric(DXDRLLE)) / 1000
    ), by = "SEQN") %>%
    mutate(Cycle = cy)

  return(df)
}

df_all <- bind_rows(fetch_cycle("G"), fetch_cycle("H")) %>%
  mutate(
    Gender = case_when(
      Gender_raw %in% c("1", "Male") ~ "Male",
      Gender_raw %in% c("2", "Female") ~ "Female",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(Age >= 20, FastHr >= 8, TG < 500, !is.na(Gender)) %>%
  mutate(
    Log2_TG = log2(TG),
    Relative_Grip = Grip / Weight_kg,
    Relative_ASM = ASM_kg / Weight_kg,
    Gender = factor(Gender, levels = c("Male", "Female")),
    Cycle = factor(Cycle, levels = c("G", "H")),
    WTMEC2YR_adj = WTMEC2YR / 2,
    WTSAF2YR_adj = WTSAF2YR / 2
  )

cat(sprintf("合并后样本量：%d\n", nrow(df_all)))
cat(sprintf("WTMEC2YR 非缺失：%d\n", sum(!is.na(df_all$WTMEC2YR_adj))))
cat(sprintf("WTSAF2YR 非缺失：%d\n", sum(!is.na(df_all$WTSAF2YR_adj))))

# ---------- 通用跑模型（修正版：nobs(fit) 取真实纳入 N）----------
run_model <- function(df, weight_var, formula, label) {
  df_sub <- df[!is.na(df[[weight_var]]) & df[[weight_var]] > 0, ]
  design <- svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA,
                      weights = as.formula(paste0("~", weight_var)),
                      data = df_sub, nest = TRUE)
  fit <- svyglm(as.formula(formula), design = design)
  tg <- tidy(fit) %>% filter(term == "Log2_TG")
  n_design <- nrow(df_sub)      # design 样本量（旧口径，仅用于对比）
  n_model  <- nobs(fit)         # 模型真实纳入量（新口径）
  cat(sprintf("[%s] N_design=%d | N_model=%d (drop %d) | TG: β=%.4f, SE=%.4f, P=%.4f\n",
              label, n_design, n_model, n_design - n_model,
              tg$estimate, tg$std.error, tg$p.value))
  return(list(fit = fit, tg = tg, n = n_model, n_design = n_design))
}

cat("\n===== 模型 1: 相对握力 =====\n")
m1a <- run_model(df_all, "WTMEC2YR_adj", "Relative_Grip ~ Log2_TG + Age + Gender + Cycle", "MEC")
m1b <- run_model(df_all, "WTSAF2YR_adj", "Relative_Grip ~ Log2_TG + Age + Gender + Cycle", "Fasting")

cat("\n===== 模型 2: 相对 ASM =====\n")
m2a <- run_model(df_all, "WTMEC2YR_adj", "Relative_ASM ~ Log2_TG + Age + Gender + Cycle", "MEC")
m2b <- run_model(df_all, "WTSAF2YR_adj", "Relative_ASM ~ Log2_TG + Age + Gender + Cycle", "Fasting")

cat("\n===== 模型 3: 绝对握力 + BMI =====\n")
m3a <- run_model(df_all, "WTMEC2YR_adj", "Grip ~ Log2_TG + Age + Gender + Cycle + BMI", "MEC")
m3b <- run_model(df_all, "WTSAF2YR_adj", "Grip ~ Log2_TG + Age + Gender + Cycle + BMI", "Fasting")

# ---------- 汇总 ----------
summary_df <- tibble(
  Model = rep(c("Relative Grip", "Relative ASM", "Abs Grip + BMI"), each = 2),
  Weight = rep(c("MEC", "Fasting"), 3),
  N = c(m1a$n, m1b$n, m2a$n, m2b$n, m3a$n, m3b$n),
  Beta = c(m1a$tg$estimate, m1b$tg$estimate,
           m2a$tg$estimate, m2b$tg$estimate,
           m3a$tg$estimate, m3b$tg$estimate),
  SE = c(m1a$tg$std.error, m1b$tg$std.error,
         m2a$tg$std.error, m2b$tg$std.error,
         m3a$tg$std.error, m3b$tg$std.error),
  P = c(m1a$tg$p.value, m1b$tg$p.value,
        m2a$tg$p.value, m2b$tg$p.value,
        m3a$tg$p.value, m3b$tg$p.value)
)

cat("\n\n========== 汇总对比（真实 N）==========\n")
print(summary_df, n = 10)

write.csv(summary_df, "results/WTSAF2YR_sensitivity.csv", row.names = FALSE)
saveRDS(list(m1a = m1a$fit, m1b = m1b$fit,
             m2a = m2a$fit, m2b = m2b$fit,
             m3a = m3a$fit, m3b = m3b$fit,
             N = setNames(summary_df$N, paste(summary_df$Model, summary_df$Weight, sep = "|")),
             N_design = c(m1a$n_design, m1b$n_design, m2a$n_design,
                          m2b$n_design, m3a$n_design, m3b$n_design)),
        "results/WTSAF2YR_sensitivity.rds")

cat("\n✅ 落盘：results/WTSAF2YR_sensitivity.csv / .rds\n")
