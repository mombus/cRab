#Вводная часть
#Цель: оценить параметры логистической селективности промысла по длине — SL50 (длина 50% улавливаемости) и SL95 (длина 95% улавливаемости) — для дальнейшего использования в length-based моделях (LBSPR, LBB, LIME).
#
#Идея: если считать, что распределение по длине в съёмке отражает доступность (availability) и наличие особей в популяции в заданный сезон/район, а промысловое распределение — это выборка тех же «доступных» особей, прошедших через «фильтр» селективности снасти, то:
#- p_survey(L) аппроксимирует истинную длиновую структуру наличности,
#- p_fish(L) пропорциональна p_survey(L) умножить на S(L),
#- S(L) — логистическая кривая селективности: S(L) = 1 / [1 + exp(-k (L ? SL50))],
#- удобно задавать её через SL50 и SR = SL95-SL50, где k = ln(19)/SR, а SL95 = SL50 + SR.
#
#Мы подбираем SL50 и SR так, чтобы вероятность наблюдать промысловые длины была максимальна (максимизация мультиномиального правдоподобия). Для устойчивости:
#- сравниваем «совпадающие» годы (пересечение лет съёмки и промысла),
#- задаём фиксированную ширину бина (например, 2 мм),
#- сглаживаем нули в p_survey маленькой добавкой.
#
#Результат: оценки SL50 и SL95 (с доверительными интервалами бутстрепом) и диагностические графики: кривая селективности на фоне «эмпирической» селективности, полученной из отношения промыслового и съёмочного распределений.
#
#Скрипт: оценка логистической селективности SL50 и SL95 (через SR)
#Скопируйте и запустите целиком. Комментарии закомментированы.
#
# estimate_selectivity_from_survey_vs_fishery.R
# Оценка логистической селективности промысла: SL50 и SL95 (= SL50 + SR)
# Метод: мультиномиальное ММП, p_fish ? p_survey * S(L; SL50, SR)
#
# ----------------------- НАСТРОЙКА СРЕДЫ И ПАКЕТОВ ----------------------------
 setwd("C:/LBM/")  # при необходимости укажите рабочую директорию

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
})

# ----------------------------- ПОЛЬЗОВАТЕЛЬСКИЕ ОПЦИИ -------------------------
sex_use     <- "M"         # пол
binw        <- 2           # ширина бина, мм
Lmin_mm     <- 40          # минимум длины для анализа
Lmax_mm     <- 280         # максимум длины для анализа
years_use   <- NULL        # например, 2016:2020; NULL = все пересечения лет
months_svy  <- NULL        # например, 8:10; NULL = все месяцы съёмки
months_fsh  <- NULL        # например, 10:12; NULL = все месяцы промысла
pool_years  <- FALSE       # если TRUE — агрегирует все годы в один срез

# Бутстреп для ДИ (необязательно; может быть долгим)
do_boot     <- TRUE
B_boot      <- 300         # число репликаций бутстрэпа
seed_boot   <- 123

# ---------------------------- ЗАГРУЗКА ДАННЫХ ---------------------------------
SURVEY  <- read.csv("SURVEYDATA.csv",  sep = ";", stringsAsFactors = FALSE)
FISHERY <- read.csv("FISHERYDATA.csv", sep = ";", stringsAsFactors = FALSE)

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

if (!is.null(months_svy)) survey  <- survey  %>% filter(MONTH %in% months_svy)
if (!is.null(months_fsh)) fishery <- fishery %>% filter(MONTH %in% months_fsh)

# ----------------------------- БИНЫ ПО ДЛИНЕ ----------------------------------
mk_bins <- function(Lmin, Lmax, bw) {
  brks <- seq(floor(Lmin/bw)*bw, ceiling(Lmax/bw)*bw, by = bw)
  mids <- brks[-1] - bw/2
  list(breaks = brks, mids = mids)
}
bins <- mk_bins(Lmin_mm, Lmax_mm, binw)

bin_counts <- function(x, breaks) {
  # Фильтруем данные по диапазону бинов
  x <- x[is.finite(x) & x >= min(breaks) & x <= max(breaks)]
  hist(x, breaks = breaks, plot = FALSE)$counts
}
# ------------------------- ГОДА ДЛЯ СРАВНЕНИЯ ---------------------------------
years_common <- intersect(unique(survey$YEAR), unique(fishery$YEAR))
if (!is.null(years_use)) years_common <- intersect(years_common, years_use)
years_common <- sort(as.integer(years_common))

if (length(years_common) == 0) stop("Нет пересечения лет между съёмкой и промыслом.")

# -------------------------- ДАННЫЕ ДЛЯ ММП ------------------------------------
# Сохраняем по каждому году: вектор p_svy (наличность) и n_fsh (мультиномиальные счётчики)
svy_list <- list()
fsh_list <- list()

if (pool_years) {
  svy_all <- survey %>% filter(YEAR %in% years_common) %>% pull(CARAPACE)
  fsh_all <- fishery %>% filter(YEAR %in% years_common) %>% pull(CARAPACE)
  svy_cnt <- bin_counts(svy_all, bins$breaks)
  fsh_cnt <- bin_counts(fsh_all, bins$breaks)
  # сглаживаем нули «+0.5»
  p_svy <- (svy_cnt + 0.5) / sum(svy_cnt + 0.5)
  svy_list[["ALL"]] <- p_svy
  fsh_list[["ALL"]] <- fsh_cnt
} else {
  for (yr in years_common) {
    svy_len <- survey %>% filter(YEAR == yr) %>% pull(CARAPACE)
    fsh_len <- fishery %>% filter(YEAR == yr) %>% pull(CARAPACE)
    svy_cnt <- bin_counts(svy_len, bins$breaks)
    fsh_cnt <- bin_counts(fsh_len, bins$breaks)
    if (sum(svy_cnt) >= 100 && sum(fsh_cnt) >= 100) {
      p_svy <- (svy_cnt + 0.5) / sum(svy_cnt + 0.5)
      svy_list[[as.character(yr)]] <- p_svy
      fsh_list[[as.character(yr)]] <- fsh_cnt
    }
  }
  if (length(svy_list) == 0) stop("Недостаточно данных по годам (мало наблюдений). Увеличьте период/объём.")
}

Lmid <- bins$mids

# ------------------------- МОДЕЛЬ СЕЛЕКТИВНОСТИ -------------------------------
# Логистическая кривая через SL50 и SR: S(L) = 1 / (1 + exp(-k(L-SL50))), где k = ln(19)/SR
logistic_sel <- function(L, SL50, SR) {
  k <- log(19) / pmax(SR, 1e-6)
  1 / (1 + exp(-k * (L - SL50)))
}

# Негативный логарифм правдоподобия (сумма по годам)
nll_selectivity <- function(par) {
  SL50 <- par[1]
  SR   <- exp(par[2])  # параметризация для SR>0
  nll <- 0
  for (nm in names(svy_list)) {
    pA <- svy_list[[nm]]                      # p_survey (наличность)
    S  <- logistic_sel(Lmid, SL50, SR)        # селективность по длине
    pF_unnorm <- pA * S
    pF <- pF_unnorm / pmax(sum(pF_unnorm), 1e-12)
    n  <- fsh_list[[nm]]
    # мультиномиальный лог-лик
    nll <- nll - sum(n * log(pmax(pF, 1e-12)))
  }
  nll
}

# ------------------------------ ИНИЦИАЛИЗАЦИЯ ---------------------------------
# Первичные догадки: SL50 — медиана промысловой длины; SR ~ 10–15 мм
if (pool_years) {
  L50_init <- median(rep(Lmid, times = fsh_list[["ALL"]]))
} else {
  all_counts <- Reduce("+", fsh_list)
  L50_init <- median(rep(Lmid, times = all_counts))
}
SR_init <- 12

par_init <- c(L50_init, log(SR_init))

# ------------------------------- ОПТИМИЗАЦИЯ ----------------------------------
fit <- optim(
  par = par_init,
  fn = nll_selectivity,
  method = "BFGS",
  hessian = TRUE,
  control = list(reltol = 1e-10, maxit = 1000)
)

SL50_est <- fit$par[1]
SR_est   <- exp(fit$par[2])
SL95_est <- SL50_est + SR_est

# Оценка стандартных ошибок из гессиана (асимптотика)
SE_SL50 <- NA_real_; SE_SR <- NA_real_
if (is.matrix(fit$hessian) && all(is.finite(fit$hessian))) {
  H <- fit$hessian
  V <- try(solve(H), silent = TRUE)
  if (!inherits(V, "try-error")) {
    # delta-method для SR = exp(theta2)
    var_SL50 <- V[1,1]
    var_t2   <- V[2,2]
    cov_1_2  <- V[1,2]
    SE_SL50  <- sqrt(max(var_SL50, 0))
    SE_SR    <- sqrt(max((exp(fit$par[2])^2) * var_t2, 0))
  }
}

# ----------------------------- БУТСТРЕП ДЛЯ ДИ --------------------------------
boot_tab <- tibble()
if (do_boot) {
  set.seed(seed_boot)
  resample_year <- function(n_vec) {
    # небараметрический бутстреп: выборка поособных длин из исходных в каждом году
    # для экономии: приближаем мультиномиалью по наблюдённым пропорциям
    N <- sum(n_vec)
    p <- (n_vec + 0.5) / sum(n_vec + 0.5)
    as.integer(rmultinom(1, size = N, prob = p)[,1])
  }
  for (b in 1:B_boot) {
    svy_list_b <- list(); fsh_list_b <- list()
    for (nm in names(svy_list)) {
      svy_list_b[[nm]] <- (resample_year((svy_list[[nm]] * 10000))) # масштаб, чтобы были целые
      svy_list_b[[nm]] <- (svy_list_b[[nm]] + 0.5) / sum(svy_list_b[[nm]] + 0.5)
      fsh_list_b[[nm]] <- resample_year(fsh_list[[nm]])
    }
    nll_b <- function(par) {
      SL50 <- par[1]; SR <- exp(par[2])
      nll <- 0
      for (nm in names(svy_list_b)) {
        pA <- svy_list_b[[nm]]
        S  <- logistic_sel(Lmid, SL50, SR)
        pF <- (pA * S)
        pF <- pF / pmax(sum(pF), 1e-12)
        n  <- fsh_list_b[[nm]]
        nll <- nll - sum(n * log(pF + 1e-12))
      }
      nll
    }
    fit_b <- try(optim(par = c(SL50_est, log(SR_est)), fn = nll_b, method = "BFGS", control = list(maxit = 500)), silent = TRUE)
    if (!inherits(fit_b, "try-error")) {
      SL50_b <- fit_b$par[1]; SR_b <- exp(fit_b$par[2])
      boot_tab <- bind_rows(boot_tab, tibble(SL50 = SL50_b, SR = SR_b, SL95 = SL50_b + SR_b))
    }
  }
}

# Доверительные интервалы
ci <- function(x) c(lo = quantile(x, 0.025, na.rm = TRUE), hi = quantile(x, 0.975, na.rm = TRUE))
SL50_CI <- if (nrow(boot_tab) > 10) ci(boot_tab$SL50) else c(lo = NA, hi = NA)
SR_CI   <- if (nrow(boot_tab) > 10) ci(boot_tab$SR)   else c(lo = NA, hi = NA)
SL95_CI <- if (nrow(boot_tab) > 10) ci(boot_tab$SL95) else c(lo = NA, hi = NA)

# ----------------------------- ВЫВОД РЕЗУЛЬТАТОВ ------------------------------
res <- tibble(
  SL50 = SL50_est,
  SL50_lo = SL50_CI["lo"], SL50_hi = SL50_CI["hi"],
  SR = SR_est,
  SR_lo = SR_CI["lo"], SR_hi = SR_CI["hi"],
  SL95 = SL95_est,
  SL95_lo = SL95_CI["lo"], SL95_hi = SL95_CI["hi"],
  SE_SL50 = SE_SL50,
  SE_SR = SE_SR,
  binw = binw,
  pool_years = pool_years
)
print(res)
write.csv(res, "selectivity_logistic_SL50_SL95.csv", row.names = FALSE)

# ----------------------------- ДИАГНОСТИКА/ГРАФИКИ ----------------------------
# ----------------------------- ДИАГНОСТИКА/ГРАФИКИ ----------------------------
# Создаем эмпирическую селективность с доверительными интервалами
emp_df <- tibble()
for (nm in names(svy_list)) {
  pA <- svy_list[[nm]]
  nF <- fsh_list[[nm]]
  pF <- nF / sum(nF)
  ratio <- pF / pmax(pA, 1e-12)
  # Рассчитываем доверительные интервалы для эмпирической селективности
  se_ratio <- sqrt((pF * (1 - pF) / sum(nF) + pA * (1 - pA) / sum(svy_cnt)) / (pA^2))
  lower <- ratio - 1.96 * se_ratio
  upper <- ratio + 1.96 * se_ratio
  # Нормируем к [0,1] для визуализации
  max_ratio <- max(ratio, na.rm = TRUE)
  ratio <- ratio / max_ratio
  lower <- lower / max_ratio
  upper <- upper / max_ratio
  
  emp_df <- bind_rows(emp_df, tibble(
    group = nm, 
    L = Lmid, 
    empS = ratio,
    lower = lower,
    upper = upper
  ))
}

# Создаем отдельный датафрейм для теоретической кривой
# Исправление: правильно вычисляем доверительные интервалы через бутстреп
if (nrow(boot_tab) > 10) {
  # Создаем более гладкую кривую для визуализации
  sel_df <- tibble(L = seq(Lmin_mm, Lmax_mm, by = 0.5)) %>%
    mutate(
      S = logistic_sel(L, SL50_est, SR_est),
      # Для каждого значения L вычисляем доверительные интервалы
      S_lo = sapply(L, function(x) {
        s_vals <- sapply(1:nrow(boot_tab), function(i) {
          logistic_sel(x, boot_tab$SL50[i], boot_tab$SR[i])
        })
        quantile(s_vals, 0.025)
      }),
      S_hi = sapply(L, function(x) {
        s_vals <- sapply(1:nrow(boot_tab), function(i) {
          logistic_sel(x, boot_tab$SL50[i], boot_tab$SR[i])
        })
        quantile(s_vals, 0.975)
      })
    )
} else {
  sel_df <- tibble(L = seq(Lmin_mm, Lmax_mm, by = 1)) %>%
    mutate(S = logistic_sel(L, SL50_est, SR_est))
}

# Создаем первый график: сравнение эмпирической и теоретической селективности
p1 <- ggplot() +
  # Эмпирические точки с доверительными интервалами
  geom_point(data = emp_df, aes(L, empS, color = group), alpha = 0.7, size = 3) +
  geom_linerange(data = emp_df, aes(L, ymin = lower, ymax = upper, color = group), alpha = 0.4, size = 1) +
  
  # Теоретическая кривая с доверительными интервалами
  geom_line(data = sel_df, aes(L, S), color = "#D81B60", linewidth = 1.8, alpha = 0.9) +
  geom_ribbon(data = sel_df, aes(L, ymin = S_lo, ymax = S_hi), fill = "#D81B60", alpha = 0.2) +
  
  # Вертикальные линии для SL50 и SL95
  geom_vline(xintercept = SL50_est, linetype = "dashed", color = "#1E88E5", linewidth = 1.2) +
  geom_vline(xintercept = SL95_est, linetype = "dotted", color = "#1E88E5", linewidth = 1.2) +
  
  # Подписи для SL50 и SL95
  annotate("text", x = SL50_est, y = 0.90, label = paste0("SL50 = ", round(SL50_est, 1), " мм"), 
           color = "black", hjust = 0, size = 4, fontface = "bold") +
  annotate("text", x = SL95_est, y = 1.03, label = paste0("SL95 = ", round(SL95_est, 1), " мм"), 
           color = "black", hjust = 0, size = 4, fontface = "bold") +
  
  # Настройки графика
  labs(
    x = "Длина карапакса, мм",
    y = "Нормированная селективность",
    title = "Сравнение эмпирической и теоретической селективности",
    subtitle = "Оценка SL50 и SL95 с доверительными интервалами (95%)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
    legend.position = "bottom",
    legend.title = element_blank(),
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 12)
  ) +
  scale_color_viridis_d(option = "plasma") +
  scale_y_continuous(limits = c(0, 1.1), breaks = seq(0, 1, by = 0.2)) +
  coord_cartesian(xlim = c(50, 170))  # Добавлено ограничение оси X до 180 мм

# Создаем второй график: кривая селективности с деталями
p2 <- ggplot(sel_df, aes(L, S)) +
  geom_line(color = "#1E88E5", linewidth = 1.8) +
  geom_ribbon(aes(ymin = S_lo, ymax = S_hi), fill = "#1E88E5", alpha = 0.2) +
  geom_vline(xintercept = SL50_est, linetype = "dashed", color = "#D81B60", linewidth = 1.2) +
  geom_vline(xintercept = SL95_est, linetype = "dotted", color = "#D81B60", linewidth = 1.2) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "#7570B3", linewidth = 1) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "#7570B3", linewidth = 1) +
  
  # Подписи для SL50 и SL95
  annotate("text", x = SL50_est + 1, y = 0.25, label = paste0("SL50 = ", round(SL50_est, 1), " мм"), 
           color = "#D81B60", size = 4, fontface = "bold") +
  annotate("text", x = SL95_est + 1, y = 0.75, label = paste0("SL95 = ", round(SL95_est, 1), " мм"), 
           color = "#D81B60", size = 4, fontface = "bold") +
  
  # Подписи для 50% и 95% уровней
  annotate("text", x = Lmin_mm + 5, y = 0.52, label = "50% улавливаемости", 
           color = "#7570B3", size = 3.5, fontface = "italic") +
  annotate("text", x = Lmin_mm + 5, y = 0.97, label = "95% улавливаемости", 
           color = "#7570B3", size = 3.5, fontface = "italic") +
  
  labs(
    x = "Длина карапакса, мм",
    y = "Вероятность улавливания (селективность)",
    title = "Логистическая кривая селективности",
    subtitle = "С доверительными интервалами (95%) и ключевыми параметрами"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 12)
  ) +
  scale_y_continuous(limits = c(0, 1.05), breaks = seq(0, 1, by = 0.2)) +
  coord_cartesian(xlim = c(50, 170))  # Добавлено ограничение оси X до 180 мм
# Сохраняем графики
ggsave("selectivity_comparison.png", p1, width = 12, height = 7, dpi = 300)
ggsave("selectivity_curve.png", p2, width = 10, height = 6, dpi = 300)

# Выводим графики в консоль
print(p1)
print(p2)
