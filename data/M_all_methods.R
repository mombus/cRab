# =============================================================================
# ОЦЕНКА ПАРАМЕТРОВ ПРАЙЕРВ ДЛЯ КОЭФФИЦИЕНТА СМЕРТНОСТИ (M) И ОТНОШЕНИЯ M/K
# ДЛЯ КРАБА-СТРИГУНА (Chionoecetes opilio) В БАРЕНЦЕВОМ МОРЕ
# =============================================================================
#
# ЦЕЛЬ: Получить обоснованные априорные распределения для M и M/K 
#       для использования в length-based методах (LBSPR, LBB, LIME)
#
# МЕТОДЫ:
# 1. Pauly (1980) - температурно-зависимая формула
# 2. Then et al. (2015) - через максимальный возраст
# 3. Lorenzen (2005) - размерно-специфичная смертность
# 4. Z-minus-F - через длинную кривую и соотношение улов/биомасса
#
# ВХОДНЫЕ ДАННЫЕ (должны быть в рабочей директории):
# - SURVEYDATA.csv        - данные съемки (длины особей)
# - FISHERYDATA.csv       - данные промысла (длины особей)  
# - CATCH.csv             - годовые уловы (тонны)
# - SURVEY_INDEX.csv      - индексы биомассы по годам
# - ELEFAN_params_constrained.csv - параметры роста (Linf, K, t_anchor)
# - LW_coeffs.csv         - коэффициенты длина-масса
#
# =============================================================================

# -------------------------- 1. ЗАГРУЗКА БИБЛИОТЕК ----------------------------
suppressPackageStartupMessages({
  library(tidyverse)   # для обработки данных и визуализации
  library(lubridate)   # для работы с датами
})

# ------------------------ 2. НАСТРОЙКИ ПОЛЬЗОВАТЕЛЯ --------------------------

# Установка рабочей директории (замените на свою)
setwd("C:/LBM/")

## 2.1 Настройки для метода Pauly
# Температурная сетка для дна Баренцева моря (°C)
T_grid <- seq(1.0, 3.0, by = 0.5)

## 2.2 Настройки для метода Then
# Квантили для оценки максимального возраста
qvec_tmax <- c(0.99, 0.995, 0.999)

## 2.3 Настройки для длинной кривой
binw <- 2                    # ширина бина (мм)
min_counts_bin <- 5          # мин. количество в бине

## 2.4 Настройки для метода Lorenzen
L_sel_min <- 110             # мин. длина для полной улавливаемости (мм)

## 2.5 Настройки для метода Z-minus-F
q_median <- 0.3              # медиана для приора q (коэффициент уловистости)
q_cv <- 0.5                  # коэффициент вариации для q
nmc_q <- 5000                # число Монте-Карло итераций

## 2.6 Коэффициенты для Lorenzen
# Стандартные значения для рыб (для крабов могут отличаться!)
c_Lor <- 3.0                 # коэффициент c в M = c * W^(-d)
d_Lor <- 0.288               # показатель степени d

## 2.7 Выбор группы для длина-масса
LW_group_use <- "SEX_M"      # "SEX_M", "SEX_F" или "ALL"

# ------------------------ 3. ЗАГРУЗКА ДАННЫХ --------------------------------

cat("Загрузка данных...\n")

# 3.1 Основные данные по длинам
SURVEY  <- read.csv("SURVEYDATA.csv", sep = ";", stringsAsFactors = FALSE)
FISHERY <- read.csv("FISHERYDATA.csv", sep = ";", stringsAsFactors = FALSE)

# 3.2 Данные по уловам и индексам
CATCH <- read.csv("CATCH.csv", sep = ";", stringsAsFactors = FALSE)
SIND  <- read.csv("SURVEY_INDEX.csv", sep = ";", stringsAsFactors = FALSE)

# 3.3 Параметры роста и длина-масса
PAR_G <- read.csv("ELEFAN_params_constrained.csv", sep = ",", stringsAsFactors = FALSE)
LW    <- read.csv("LW_coeffs.csv", sep = ",", stringsAsFactors = FALSE)

# -------------------- 4. ПРЕДВАРИТЕЛЬНАЯ ОБРАБОТКА --------------------------

cat("Предварительная обработка данных...\n")

# 4.1 Очистка и стандартизация данных по длинам
clean_length_data <- function(df) {
  df %>%
    mutate(
      SEX = case_when(
        toupper(SEX) %in% c("M", "MALE") ~ "M",
        toupper(SEX) %in% c("F", "FEMALE") ~ "F",
        TRUE ~ NA_character_
      ),
      date = ymd(paste(YEAR, MONTH, 15, sep = "-")),
      CARAPACE = as.numeric(CARAPACE)
    ) %>%
    filter(!is.na(SEX), !is.na(CARAPACE), CARAPACE > 0)
}

survey <- clean_length_data(SURVEY)
fishery <- clean_length_data(FISHERY)

# 4.2 Очистка данных по уловам и индексам
CATCH <- CATCH %>%
  mutate(
    YEAR = as.integer(YEAR),
    CATCH_T = as.numeric(gsub(",", ".", CATCH_T))
  ) %>%
  filter(is.finite(YEAR), is.finite(CATCH_T), CATCH_T > 0)

SIND <- SIND %>%
  mutate(
    YEAR = as.integer(YEAR),
    INDEX_T = as.numeric(gsub(",", ".", INDEX_T))
  ) %>%
  filter(is.finite(YEAR), is.finite(INDEX_T), INDEX_T > 0)

# 4.3 Извлечение параметров роста
Linf <- as.numeric(PAR_G$Linf[1])
K    <- as.numeric(PAR_G$K[1])
t0   <- as.numeric(PAR_G$t_anchor[1])  # используем t_anchor как приближение t0

cat(sprintf("Параметры роста: Linf = %.1f мм, K = %.3f год??, t0 ? %.2f\n", 
            Linf, K, t0))

# 4.4 Выбор данных по самцам (основной объект промысла)
survey_m <- survey %>% filter(SEX == "M")
fish_m   <- fishery %>% filter(SEX == "M")

# ------------------------ 5. ВСПОМОГАТЕЛЬНЫЕ ФУНКЦИИ ------------------------

# 5.1 Обратное преобразование фон Берталанфи: длина > возраст
t_from_L <- function(L, Linf, K, t0 = 0) {
  Lc <- pmin(L, Linf * 0.999)  # защита от деления на ноль
  t0 - (1 / K) * log(1 - Lc / Linf)
}

# 5.2 Формула Pauly (1980) для естественной смертности
# Примечание: Linf должен быть в см! Конвертируем мм > см.
M_pauly <- function(Linf_mm, K, T) {
  Linf_cm <- Linf_mm / 10  # конвертация в см
  10^(-0.0066 - 0.279 * log10(Linf_cm) + 0.6543 * log10(K) + 0.4634 * log10(T))
}

# 5.3 Формула Then et al. (2015) через максимальный возраст
M_then_tmax <- function(tmax) {
  4.899 * tmax^(-0.916)
}

# 5.4 Размерно-специфичная смертность по Lorenzen
M_lorenzen_from_L <- function(L, a, b, c = 3.0, d = 0.288) {
  W <- a * (L^b)          # масса от длины
  c * (W^(-d))            # смертность от массы
}

# 5.5 Оценка Z через length-converted catch curve
Z_from_length_converted <- function(lengths, binw, Linf, K, t0 = 0, min_counts = 5) {
  lengths <- lengths[is.finite(lengths)]
  if (length(lengths) < 50) return(NA_real_)
  
  # Создание бинов
  Lmin <- floor(min(lengths) / binw) * binw
  Lmax <- ceiling(max(lengths) / binw) * binw
  if (Lmax <= Lmin) return(NA_real_)
  
  breaks <- seq(Lmin, Lmax, by = binw)
  if (length(breaks) < 5) return(NA_real_)
  
  # Гистограмма
  h <- hist(lengths, breaks = breaks, plot = FALSE)
  n <- h$counts
  if (sum(n) < 100) return(NA_real_)
  
  # Расчет возраста для середин бинов
  mids <- breaks[-1] - binw/2
  t_mid <- t_from_L(mids, Linf, K, t0)
  
  # Время пребывания в бине
  L_low  <- mids - binw/2
  L_high <- mids + binw/2
  t_low  <- t_from_L(L_low, Linf, K, t0)
  t_high <- t_from_L(L_high, Linf, K, t0)
  dt <- pmax(t_high - t_low, 1e-6)
  
  # Нисходящая ветвь (правее моды)
  i_mode <- which.max(n)
  idx <- seq(i_mode + 1, length(n))
  idx <- idx[n[idx] >= min_counts]
  if (length(idx) < 5) return(NA_real_)
  
  # Линейная регрессия в полулогарифмическом масштабе
  y <- log(n[idx] / dt[idx])
  x <- t_mid[idx]
  fit <- try(lm(y ~ x), silent = TRUE)
  if (inherits(fit, "try-error")) return(NA_real_)
  
  Z <- -as.numeric(coef(fit)[2])
  if (!is.finite(Z) || Z <= 0) return(NA_real_)
  
  return(Z)
}

# ------------------------ 6. РАСЧЕТ M ПО РАЗНЫМ МЕТОДАМ ----------------------

cat("Расчет естественной смертности разными методами...\n")

## 6.1 МЕТОД PAULY (1980) - температурно-зависимый
cat("  - Метод Pauly (температурная формула)\n")
M_Pauly <- tibble(
  source = "Pauly1980",
  T = T_grid,
  M = M_pauly(Linf, K, T_grid)
)

# Агрегирование: среднее по температурной сетке
M_Pauly_agg <- M_Pauly %>% 
  summarise(M = mean(M)) %>% 
  mutate(source = "Pauly1980_mean")

## 6.2 МЕТОД THEN (2015) - через максимальный возраст
cat("  - Метод Then (через максимальный возраст)\n")
Lq <- quantile(fish_m$CARAPACE, probs = qvec_tmax, na.rm = TRUE, names = FALSE)
tq <- t_from_L(Lq, Linf, K, t0)

M_Then <- tibble(
  source = "Then2015",
  q = qvec_tmax,
  Lq = as.numeric(Lq),
  tmax = as.numeric(tq),
  M = M_then_tmax(as.numeric(tq))
)

# Агрегирование: среднее по квантилям
M_Then_agg <- M_Then %>% 
  summarise(M = mean(M)) %>% 
  mutate(source = "Then2015_mean")

## 6.3 МЕТОД LORENZEN (2005) - размерно-специфичный
cat("  - Метод Lorenzen (размерно-специфичная смертность)\n")
M_Lor_agg <- tibble()

# Выбор коэффициентов длина-масса
LW_use <- LW %>% filter(GROUP == LW_group_use)
if (nrow(LW_use) == 0) {
  cat("  Предупреждение: выбранная группа не найдена, используем 'ALL'\n")
  LW_use <- LW %>% filter(GROUP == "ALL")
}

if (nrow(LW_use) > 0) {
  a_LW <- as.numeric(LW_use$a[1])
  b_LW <- as.numeric(LW_use$b[1])
  
  # Только особи с полной улавливаемостью
  L_sel <- fish_m$CARAPACE[fish_m$CARAPACE >= L_sel_min]
  L_sel <- L_sel[is.finite(L_sel)]
  
  if (length(L_sel) >= 50) {
    M_vals <- M_lorenzen_from_L(L_sel, a = a_LW, b = b_LW, c = c_Lor, d = d_Lor)
    M_Lor <- tibble(source = "Lorenzen", L = L_sel, M = M_vals)
    
    # Агрегирование: медиана по размерам
    M_Lor_agg <- M_Lor %>% 
      summarise(M = median(M)) %>% 
      mutate(source = "Lorenzen_median")
  }
}

## 6.4 МЕТОД Z-MINUS-F - через длинную кривую и улов/биомассу
cat("  - Метод Z-minus-F (через длинную кривую)\n")
M_ZF_agg <- tibble()

if (nrow(CATCH) > 0 && nrow(SIND) > 0) {
  # Оценка Z по годам
  Z_by_year <- survey_m %>%
    group_by(YEAR) %>%
    summarise(
      Z = Z_from_length_converted(CARAPACE, binw, Linf, K, t0, min_counts_bin),
      n = n(), 
      .groups = "drop"
    ) %>%
    filter(!is.na(Z), n >= 200)
  
  if (nrow(Z_by_year) > 0) {
    # Объединение с уловами и индексами
    df_cf <- Z_by_year %>%
      inner_join(CATCH, by = "YEAR") %>%
      inner_join(SIND, by = "YEAR") %>%
      filter(is.finite(Z), CATCH_T > 0, INDEX_T > 0)
    
    if (nrow(df_cf) > 0) {
      # Монте-Карло симуляция для q
      sdlog <- sqrt(log(1 + q_cv^2))
      meanlog <- log(q_median)
      set.seed(123)
      q_draws <- rlnorm(nmc_q, meanlog = meanlog, sdlog = sdlog)
      
      # Расчет M для каждого года
      M_mc <- map_df(seq_len(nrow(df_cf)), function(i) {
        yr  <- df_cf$YEAR[i]
        Z   <- df_cf$Z[i]
        C_y <- df_cf$CATCH_T[i]
        I_y <- df_cf$INDEX_T[i]
        
        B_draws <- I_y / q_draws
        U_draws <- pmin(C_y / B_draws, 0.95)
        F_draws <- -log(pmax(1 - U_draws, 1e-6))
        M_draws <- pmax(Z - F_draws, 0)
        
        tibble(YEAR = yr, M = M_draws)
      })
      
      # Агрегирование: медиана по всем симуляциям
      M_ZF_agg <- M_mc %>% 
        summarise(M = median(M)) %>% 
        mutate(source = "Z_minus_F_median")
    }
  }
}

# ------------------------ 7. ФОРМИРОВАНИЕ ПРИОРА ----------------------------

cat("Формирование априорного распределения...\n")

# 7.1 Объединение агрегированных оценок
M_aggregated <- bind_rows(
  M_Pauly_agg,
  M_Then_agg,
  M_Lor_agg,
  M_ZF_agg
) %>%
filter(is.finite(M), M > 0)

# 7.2 Сохранение сырых агрегированных оценок
write.csv(M_aggregated, "M_estimates_aggregated.csv", row.names = FALSE)

# 7.3 Построение логнормального приора для M
if (nrow(M_aggregated) >= 2) {
  logM <- log(M_aggregated$M)
  mu <- mean(logM, na.rm = TRUE)
  sd <- sd(logM, na.rm = TRUE)
  
  # Функция квантилей логнормального распределения
  qfun <- function(p) exp(mu + sd * qnorm(p))
  
  M_prior_summary <- tibble(
    M_med = qfun(0.5),     # медиана
    M_lo  = qfun(0.025),   # нижний 95% ДИ
    M_hi  = qfun(0.975)    # верхний 95% ДИ
  )
  
  # 7.4 Расчет приора для M/K
  set.seed(321)
  MK_draws <- tibble(M = exp(rnorm(10000, mu, sd))) %>%
    mutate(MK = M / K)
  
  MK_prior_summary <- MK_draws %>%
    summarise(
      MK_med = median(MK, na.rm = TRUE),
      MK_lo  = quantile(MK, 0.025, na.rm = TRUE),
      MK_hi  = quantile(MK, 0.975, na.rm = TRUE)
    )
  
  # 7.5 Сохранение результатов
  write.csv(M_prior_summary, "M_prior_summary.csv", row.names = FALSE)
  write.csv(MK_prior_summary, "MK_prior_summary.csv", row.names = FALSE)
  
  # ------------------------ 8. ВИЗУАЛИЗАЦИЯ И ОТЧЕТ --------------------------
  
  cat("Подготовка отчета...\n")
  
  # 8.1 Сводка по методам
  method_summary <- M_aggregated %>%
    mutate(Method = case_when(
      source == "Pauly1980_mean" ~ "Pauly (1980) - температурная формула",
      source == "Then2015_mean" ~ "Then et al. (2015) - через tmax",
      source == "Lorenzen_median" ~ "Lorenzen (2005) - размерно-специфичная",
      source == "Z_minus_F_median" ~ "Z-minus-F - через длинную кривую"
    )) %>%
    select(Method, M)
  
  # 8.2 Визуализация распределения оценок M
  library(ggplot2)
  
 p <- ggplot(M_aggregated, aes(x = source, y = M, fill = source)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = round(M, 3)), vjust = -0.5) +
    labs(title = "Сравнение оценок естественной смертности (M) разными методами",
         x = "Метод оценки",
         y = "Коэффициент естественной смертности (M, год??)") +
    theme_minimal() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave("M_methods_comparison.png", p, width = 8, height = 6, dpi = 300)
  
  # 8.3 Итоговый отчет
  cat("================================================================\n")
  cat("ИТОГОВЫЙ ОТЧЕТ ПО ОЦЕНКЕ ЕСТЕСТВЕННОЙ СМЕРТНОСТИ\n")
  cat("================================================================\n")
  cat(sprintf("Вид: краб-стригун (Chionoecetes opilio)\n"))
  cat(sprintf("Район: Баренцево море\n"))
  cat(sprintf("Параметры роста: Linf = %.1f мм, K = %.3f год??\n", Linf, K))
  cat("\n")
  cat("ОЦЕНКИ M ПО МЕТОДАМ:\n")
  print(method_summary)
  cat("\n")
  cat("ФИНАЛЬНЫЙ ПРИОР ДЛЯ M (логнормальное распределение):\n")
  print(M_prior_summary)
  cat("\n")
  cat("ФИНАЛЬНЫЙ ПРИОР ДЛЯ M/K (логнормальное распределение):\n")
  print(MK_prior_summary)
  cat("\n")
  cat("График сравнения методов сохранен в файле: M_methods_comparison.png\n")
  cat("Все данные сохранены в CSV-файлах в рабочей директории.\n")
  cat("================================================================\n")
  
} else {
  cat("ВНИМАНИЕ: Недостаточно данных для построения приора!\n")
  cat("Проверьте входные данные и настройки.\n")
}


p


