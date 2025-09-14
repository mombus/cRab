# run_length_based_methods_clean.R
# Скрипт для запуска LBSPR и простых индикаторов состояния запаса

# ------------------------------ РАБОЧАЯ ДИРЕКТОРИЯ ----------------------------
setwd("C:/LBM/")  # <-- УКАЖИТЕ СВОЮ ПАПКУ С ДАННЫМИ

# ------------------------------ ОПЦИИ R ---------------------------------------
options(stringsAsFactors = FALSE)

# ------------------------------ УСТАНОВКА/ЗАГРУЗКА ПАКЕТОВ --------------------
need_pkgs <- c("tidyverse","lubridate")
to_install <- setdiff(need_pkgs, rownames(installed.packages()))
if (length(to_install) > 0) install.packages(to_install, dependencies = TRUE)
suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
})

# Установка LBSPR если нужно
if (!requireNamespace("LBSPR", quietly = TRUE)) {
  install.packages("LBSPR")
}

# ------------------------------ ПАРАМЕТРЫ ПОЛЬЗОВАТЕЛЯ ------------------------
sex_use   <- "M"   # пол для анализа
binw      <- 2     # ширина бина по длине (мм)
Lmin_mm   <- 40    # мин длина для анализа (мм)
Lmax_mm   <- 220   # макс длина для анализа (мм)

# Параметры зрелости (логистическая зрелость по длине, мм)
L50_mat_mm <- 73
L95_mat_mm <- 126

# Параметры селективности промысла (логистическая, мм)
SL50_mm <- 96
SL95_mm <- 112

# CV по длинам (используется LBSPR)
CVLinf <- 0.10
CVlen  <- 0.10

# ------------------------------ ЗАГРУЗКА ФАЙЛОВ --------------------------------
SURVEY  <- read.csv("SURVEYDATA.csv",  sep = ";")
FISHERY <- read.csv("FISHERYDATA.csv", sep = ";")
CATCH   <- read.csv("CATCH.csv",       sep = ";")
SINDEX  <- read.csv("SURVEY_INDEX.csv", sep = ";")
GROWTH  <- read.csv("ELEFAN_params_constrained.csv", sep = ",")
MK_PRI  <- read.csv("MK_prior_summary.csv", sep = ",")
M_PRI   <- read.csv("M_prior_summary.csv",  sep = ",")
LWCOEF  <- read.csv("LW_coeffs.csv", sep = ",")

# ------------------------------ ПРИВЕДЕНИЕ ТИПОВ ------------------------------
survey <- SURVEY %>%
  mutate(
    SEX      = if_else(SEX %in% c("M","male","m"), "M",
                if_else(SEX %in% c("F","female","f"), "F", NA_character_)),
    CARAPACE = as.numeric(CARAPACE),
    date     = ymd(paste(YEAR, MONTH, 15, sep = "-"))
  ) %>%
  filter(SEX == sex_use, is.finite(CARAPACE), CARAPACE > 0)

fishery <- FISHERY %>%
  mutate(
    SEX      = if_else(SEX %in% c("M","male","m"), "M",
                if_else(SEX %in% c("F","female","f"), "F", NA_character_)),
    CARAPACE = as.numeric(CARAPACE),
    date     = ymd(paste(YEAR, MONTH, 15, sep = "-"))
  ) %>%
  filter(SEX == sex_use, is.finite(CARAPACE), CARAPACE > 0)

CATCH <- CATCH %>%
  mutate(
    YEAR = as.integer(YEAR),
    CATCH_T = suppressWarnings(as.numeric(gsub(",", ".", CATCH_T)))
  ) %>%
  filter(is.finite(YEAR), is.finite(CATCH_T))

SINDEX <- SINDEX %>%
  mutate(YEAR = as.integer(YEAR),
         INDEX_T = as.numeric(INDEX_T)) %>%
  filter(is.finite(YEAR), is.finite(INDEX_T))

# ------------------------------ ПАРАМЕТРЫ РОСТА/ПРИОРЫ ------------------------
Linf_mm <- as.numeric(GROWTH$Linf[1])
K_vbgf  <- as.numeric(GROWTH$K[1])
MK_med  <- if ("MK_med" %in% names(MK_PRI)) MK_PRI$MK_med[1] else NA_real_
M_med   <- if ("M_med"  %in% names(M_PRI))  M_PRI$M_med[1]  else NA_real_

# Длина–масса
LW_use <- LWCOEF %>% filter(GROUP %in% c("SEX_M","ALL")) %>% slice(1)
a_w <- as.numeric(LW_use$a[1])
b_w <- as.numeric(LW_use$b[1])

# ------------------------------ ХЕЛПЕРЫ ---------------------------------------
mk_bins <- function(Lmin, Lmax, bw) {
  brks <- seq(floor(Lmin/bw)*bw, ceiling(Lmax/bw)*bw, by = bw)
  mids <- brks[-1] - bw/2
  list(breaks = brks, mids = mids)
}

bin_lengths <- function(x, breaks) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(rep(0L, length(breaks) - 1))
  
  # Фильтруем значения, которые не попадают в диапазон breaks
  x <- x[x >= min(breaks) & x <= max(breaks)]
  if (length(x) == 0) return(rep(0L, length(breaks) - 1))
  
  hist(x, breaks = breaks, plot = FALSE)$counts
}

# ------------------------------ BINS И КОМПОЗИЦИИ -----------------------------
bins <- mk_bins(Lmin_mm, Lmax_mm, binw)

# ------------------------------ LBSPR -----------------------------------------
cat("\n=====================================\n")
cat("       АНАЛИЗ МЕТОДОМ LBSPR\n")
cat("=====================================\n")

if (requireNamespace("LBSPR", quietly = TRUE) && is.finite(MK_med)) {
  library(LBSPR)
  
  # Фильтруем данные в пределах Lmin_mm и Lmax_mm
  fishery_filt <- fishery %>% 
    filter(CARAPACE >= Lmin_mm, CARAPACE <= Lmax_mm)
  
  years_fish <- sort(unique(fishery_filt$YEAR))
  
  # Создаем матрицу длинных композиций
  LF_mat <- sapply(years_fish, function(yr) {
    vals <- fishery_filt$CARAPACE[fishery_filt$YEAR == yr]
    bin_lengths(vals, breaks = bins$breaks)
  })
  
  # Проверяем, что есть данные
  if (ncol(LF_mat) > 0 && sum(LF_mat) > 0) {
    colnames(LF_mat) <- as.character(years_fish)
    write.csv(data.frame(Lmid = bins$mids, LF_mat), "LBSPR_length_comps.csv", row.names = FALSE)
    
    LBpar <- new("LB_pars")
    LBpar@Linf     <- Linf_mm
    LBpar@MK       <- MK_med
    LBpar@L50      <- L50_mat_mm
    LBpar@L95      <- L95_mat_mm
    LBpar@SL50     <- SL50_mm
    LBpar@SL95     <- SL95_mm
    LBpar@BinWidth <- binw
    LBpar@CVLinf   <- CVLinf
    LBpar@Species  <- "Chionoecetes opilio (M)"
    
    LBlen <- new("LB_lengths")
    LBlen@LMids  <- bins$mids
    LBlen@LData  <- LF_mat
    LBlen@Years  <- years_fish
    LBlen@NYears <- length(years_fish)
    
    Control <- list(modtype = "GTG")
    
    tryCatch({
      fit_lb <- LBSPRfit(LBpar, LBlen, Control = Control)
      
      SPR_by_year <- data.frame(
        YEAR = years_fish,
        SPR  = fit_lb@SPR,
        FM   = fit_lb@FM
      )
      write.csv(SPR_by_year, "LBSPR_results_by_year.csv", row.names = FALSE)
      
      # Расчет дополнительных метрик
      mean_SPR <- mean(fit_lb@SPR, na.rm=TRUE)
      mean_FM  <- mean(fit_lb@FM, na.rm=TRUE)
      last_SPR <- tail(fit_lb@SPR, 1)
      last_FM  <- tail(fit_lb@FM, 1)
      
      # Классификация состояния по SPR
      status_SPR <- case_when(
        mean_SPR < 0.2 ~ "Критический перелов",
        mean_SPR < 0.3 ~ "Сильный перелов",
        mean_SPR < 0.4 ~ "Умеренный перелов",
        TRUE ~ "Устойчивая эксплуатация"
      )
      
      cat("\n? LBSPR завершен успешно\n")
      cat("---------------------------------\n")
      cat("Средний SPR:", round(mean_SPR, 3), "\n")
      cat("Средний F/M:", round(mean_FM, 3), "\n")
      cat("Последний год SPR:", round(last_SPR, 3), "\n")
      cat("Последний год F/M:", round(last_FM, 3), "\n")
      cat("Статус по SPR:", status_SPR, "\n")
      
    }, error = function(e) {
      cat("? Ошибка в LBSPR:", conditionMessage(e), "\n")
    })
  } else {
    cat("? LBSPR: нет данных для анализа\n")
  }
} else {
  cat("? LBSPR: пакет не найден или нет MK_med.\n")
}

# ------------------------------ ПРОСТАЯ ОЦЕНКА НА ОСНОВЕ ДЛИН -----------------
cat("\n=====================================\n")
cat("    ИНДИКАТОРЫ СОСТОЯНИЯ ЗАПАСА\n")
cat("=====================================\n")

# Анализ для всего временного ряда
years_all <- sort(unique(fishery$YEAR))
indicators_by_year <- list()

for (yr in years_all) {
  fish_yr <- fishery %>% 
    filter(YEAR == yr, CARAPACE >= Lmin_mm, CARAPACE <= Lmax_mm)
  
  if (nrow(fish_yr) >= 30) {  # минимум 30 измерений
    indicators_by_year[[as.character(yr)]] <- data.frame(
      YEAR = yr,
      N = nrow(fish_yr),
      Lmean = mean(fish_yr$CARAPACE),
      Lmed = median(fish_yr$CARAPACE),
      L95 = quantile(fish_yr$CARAPACE, 0.95),
      Lmax = max(fish_yr$CARAPACE),
      Pmega = mean(fish_yr$CARAPACE >= 0.9 * Linf_mm) * 100  # % мегаспаунеров
    )
  }
}

if (length(indicators_by_year) > 0) {
  indicators_df <- bind_rows(indicators_by_year) %>%
    mutate(
      Lmean_Linf = Lmean / Linf_mm,
      L95_Linf = L95 / Linf_mm,
      Lmax_Linf = Lmax / Linf_mm,
      Lc_Lmat = SL50_mm / L50_mat_mm  # селективность относительно зрелости
    )
  
  write.csv(indicators_df, "Length_indicators_by_year.csv", row.names = FALSE)
  
  # Тренды последних лет
  n_recent <- min(5, nrow(indicators_df))
  recent_trend <- tail(indicators_df, n_recent)
  
  # Линейная регрессия для трендов
  if (nrow(recent_trend) >= 3) {
    trend_Lmean <- coef(lm(Lmean ~ YEAR, data = recent_trend))[2]
    trend_sign <- if_else(trend_Lmean > 0, "растет", "снижается")
    
    cat("\nТренды за последние", n_recent, "лет:\n")
    cat("---------------------------------\n")
    cat("Средняя длина", trend_sign, "на", round(trend_Lmean, 2), "мм/год\n")
    cat("Текущая Lmean/Linf:", round(tail(indicators_df$Lmean_Linf, 1), 3), "\n")
    cat("Доля мегаспаунеров:", round(tail(indicators_df$Pmega, 1), 1), "%\n")
  }
}

# Агрегированная оценка за последние 3 года
cat("\n=====================================\n")
cat("  ОЦЕНКА ЗА ПОСЛЕДНИЕ 3 ГОДА\n")
cat("=====================================\n")

n_last_years <- 3
yrs_recent <- tail(sort(unique(fishery$YEAR)), n_last_years)
fish_recent <- fishery %>% 
  filter(YEAR %in% yrs_recent, CARAPACE >= Lmin_mm, CARAPACE <= Lmax_mm)

if (nrow(fish_recent) > 0) {
  # Рассчитываем индикаторы
  Lmean <- mean(fish_recent$CARAPACE)
  Lmed <- median(fish_recent$CARAPACE)
  L25 <- quantile(fish_recent$CARAPACE, 0.25)
  L75 <- quantile(fish_recent$CARAPACE, 0.75)
  L95 <- quantile(fish_recent$CARAPACE, 0.95)
  Lmax_obs <- max(fish_recent$CARAPACE)
  
  # Оценка относительно Linf и пороговых значений
  Lmean_Linf <- Lmean / Linf_mm
  L95_Linf <- L95 / Linf_mm
  Lmax_Linf <- Lmax_obs / Linf_mm
  
  # Froese индикаторы (Froese, 2004)
  Pmat <- mean(fish_recent$CARAPACE >= L50_mat_mm) * 100  # % зрелых
  Popt <- mean(fish_recent$CARAPACE >= 0.9 * Linf_mm) * 100  # % оптимального размера
  Pmega <- mean(fish_recent$CARAPACE >= Linf_mm) * 100  # % мегаспаунеров
  
  # Cope & Punt (2009) индикаторы
  LBI <- Lmean / Linf_mm  # Length-based indicator
  
  # Простая оценка состояния
  status <- case_when(
    Lmean_Linf < 0.5 ~ "Критический перелов",
    Lmean_Linf < 0.6 ~ "Сильный перелов",
    Lmean_Linf < 0.7 ~ "Умеренный перелов",
    TRUE ~ "Устойчивая эксплуатация"
  )
  
  # Оценка по Froese (2004)
  froese_status <- case_when(
    Pmat < 50 ~ "Слишком много молоди в улове",
    Pmat > 90 & Popt < 10 ~ "Перелов крупных особей",
    Popt > 30 & Pmega > 10 ~ "Хорошее состояние",
    TRUE ~ "Требует внимания"
  )
  
  # Сохраняем результаты
  simple_assessment <- data.frame(
    Years = paste(yrs_recent, collapse="-"),
    N_samples = nrow(fish_recent),
    Lmean = round(Lmean, 1),
    Lmedian = round(Lmed, 1),
    L25 = round(L25, 1),
    L75 = round(L75, 1),
    L95 = round(L95, 1),
    Lmax_obs = round(Lmax_obs, 1),
    Linf = Linf_mm,
    Lmean_Linf = round(Lmean_Linf, 3),
    L95_Linf = round(L95_Linf, 3),
    Lmax_Linf = round(Lmax_Linf, 3),
    Pmat = round(Pmat, 1),
    Popt = round(Popt, 1),
    Pmega = round(Pmega, 1),
    Status = status,
    Froese_status = froese_status
  )
  
  write.csv(simple_assessment, "Simple_assessment.csv", row.names = FALSE)
  
  cat("\nИндикаторы состояния:\n")
  cat("---------------------------------\n")
  cat("Lmean/Linf:", round(Lmean_Linf, 3), "\n")
  cat("L95/Linf:", round(L95_Linf, 3), "\n")
  cat("% зрелых (>L50):", round(Pmat, 1), "%\n")
  cat("% оптимальных (>0.9*Linf):", round(Popt, 1), "%\n")
  cat("% мегаспаунеров (>Linf):", round(Pmega, 1), "%\n")
  cat("\nЗаключение:\n")
  cat("Статус по Lmean/Linf:", status, "\n")
  cat("Статус по Froese:", froese_status, "\n")
  
  # Сохраняем длинные композиции
  LF_recent <- bin_lengths(fish_recent$CARAPACE, breaks = bins$breaks)
  write.csv(data.frame(Lmid = bins$mids, Count = LF_recent), "Length_comps_recent.csv", row.names = FALSE)
  
  # Создаем сводную таблицу для отчета
  summary_table <- data.frame(
    Indicator = c("SPR (LBSPR)", "F/M (LBSPR)", "Lmean/Linf", "% зрелых", "% оптимальных"),
    Value = c(
      if (exists("last_SPR")) round(last_SPR, 3) else NA,
      if (exists("last_FM")) round(last_FM, 2) else NA,
      round(Lmean_Linf, 3),
      round(Pmat, 1),
      round(Popt, 1)
    ),
    Target = c("?0.4", "?1.0", "?0.7", "?90%", "?30%"),
    Status = c(
      if (exists("last_SPR")) if_else(last_SPR >= 0.4, "OK", "Перелов") else NA,
      if (exists("last_FM")) if_else(last_FM <= 1.0, "OK", "Перелов") else NA,
      if_else(Lmean_Linf >= 0.7, "OK", "Перелов"),
      if_else(Pmat >= 90, "OK", "Низкая"),
      if_else(Popt >= 30, "OK", "Низкая")
    )
  )
  
  write.csv(summary_table, "Assessment_summary_table.csv", row.names = FALSE)
  
} else {
  cat("? Нет данных для простой оценки\n")
}

# ------------------------------ ИТОГОВОЕ СООБЩЕНИЕ -----------------------------
cat("\n=====================================\n")
cat("        АНАЛИЗ ЗАВЕРШЕН\n")
cat("=====================================\n")
cat("Созданные файлы:\n")
cat("---------------------------------\n")

if (file.exists("LBSPR_results_by_year.csv")) {
  cat("? LBSPR_results_by_year.csv - результаты LBSPR по годам\n")
  cat("? LBSPR_length_comps.csv - длинные композиции для LBSPR\n")
}
if (file.exists("Length_indicators_by_year.csv")) {
  cat("? Length_indicators_by_year.csv - индикаторы по годам\n")
}
if (file.exists("Simple_assessment.csv")) {
  cat("? Simple_assessment.csv - итоговая оценка\n")
  cat("? Length_comps_recent.csv - длинные композиции за последние годы\n")
}
if (file.exists("Assessment_summary_table.csv")) {
  cat("? Assessment_summary_table.csv - сводная таблица результатов\n")
}

cat("\n=====================================\n")
cat("Рекомендации по интерпретации:\n")
cat("---------------------------------\n")
cat("SPR ? 0.4: устойчивый промысел\n")
cat("F/M ? 1.0: умеренная смертность\n")
cat("Lmean/Linf ? 0.7: здоровая структура\n")
cat("% зрелых ? 90%: хорошая селективность\n")
cat("% оптимальных ? 30%: сбалансированный промысел\n")
cat("=====================================\n")
