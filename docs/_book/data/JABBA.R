# =============================================================
# ОЦЕНКА ЗАПАСА С ПОМОЩЬЮ JABBA (JUST ANOTHER BAYESIAN BIOMASS ASSESSMENT)
# =============================================================

# ------------------------- 1. ПОДГОТОВКА СРЕДЫ ---------------------------

# Загрузка необходимых пакетов
library(JABBA)   # Основной пакет для оценки запасов
library(reshape2) # Для преобразования данных (функция dcast)

# Установка рабочей директории (папки, где будут храниться результаты)
setwd("C:/BAKANEV/JABBA")

# Создание папки для результатов анализа
assessment <- "NEW JABBA"         # Название оценки
output.dir <- file.path(getwd(), assessment) # Создание пути к папке
dir.create(output.dir, showWarnings = FALSE) # Создание папки (если не существует)
setwd(output.dir) # Переход в созданную папку

# ------------------------- 2. ЗАГРУЗКА И ПОДГОТОВКА ДАННЫХ ---------------------------

# Создание вектора лет анализа
Year <- 2005:2024 # Последовательность лет от 2005 до 2024

# Вектор данных об уловах (catch)
Catch <- c(5,7,6,10,14,25,28,30,32,35,25,20,15,12,10,12,10,13,11,12)

# Вектор данных индекса обилия CPUE (catch per unit effort)
CPUE <- c(27.4,26.8,16.8,23.0,29.0,30.0,16.5,17.2,10.5,14.6,8.3,11.4,15.5,13.8,11.5,15.3,12.2,15.6,16.2,13.4)

# Вектор данных индекса обилия BESS
BESS <- c(NA,16.3,20.7,15.1,18.6,16.0,13.8,13.3,11.7,11.8,9.3,7.1,8.0,9.2,10.3,9.8,10.3,11.7,13.7,13.4)

# Форматирование данных в таблицы для JABBA
catch_data <- data.frame(year = Year, catch = Catch) # Таблица уловов
cpue_data <- data.frame(year = Year, CPUE = CPUE, BESS = BESS) # Таблица индексов

# Расчет стандартных ошибок (SE) для индексов
# Используем коэффициент вариации (CV) = 20% (0.2)
# ifelse() заменяет NA на NA, а для остальных - 0.2
se_data <- data.frame(
  year = Year,
  CPUE = ifelse(is.na(CPUE), NA, 0.2),
  BESS = ifelse(is.na(BESS), NA, 0.2)
)

# ------------------- 3. НАСТРОЙКА И ЗАПУСК МОДЕЛИ JABBA --------------------

# Создание входных данных для модели
jbinput <- build_jabba(
  catch = catch_data,    # Данные об уловах
  cpue = cpue_data,      # Данные индексов обилия
  se = se_data,          # Стандартные ошибки
  assessment = assessment, # Название оценки
  scenario = "SPiCT_adapted", # Сценарий модели
  model.type = "Schaefer", # Тип модели (Шефера)
  sigma.est = TRUE,      # Оценивать изменчивость процесса?
  r.prior = c(0.2, 0.5),    # Априорное распределение для r (среднее, SD)
  K.prior = c(189.6, 0.795),# Априорное для K (среднее, SD) 
  psi.prior = c(0.75, 0.25),# Априорное для начального заполнения
  igamma = c(0.001, 0.001), # Параметры для дисперсии наблюдений
  verbose = FALSE        # Отключить подробный вывод
)

# Запуск Байесовской модели (MCMC)
fit <- fit_jabba(
  jbinput,   # Входные данные
  ni = 50000, # Общее количество итераций
  nb = 10000, # Количество "выжигаемых" итераций (burn-in)
  nt = 5,     # Частота прореживания (thinning)
  nc = 2      # Количество цепей MCMC
)

# -------------------  вывод таблицы параметров --------------------
fit$pars


# ------------------- 4. ДИАГНОСТИКА И ВИЗУАЛИЗАЦИЯ РЕЗУЛЬТАТОВ --------------------

# Генерация стандартных диагностических графиков
jbplot_ensemble(fit)
jabba_plots(fit, output.dir = output.dir)

# Индивидуальные графики для детального анализа:
jbplot_catch(fit)       # График уловов
jbplot_cpuefits(fit)    # Сравнение модельных и наблюдаемых индексов
jbplot_ppdist(fit)      # Распределения априорных и апостериорных параметров
jbplot_residuals(fit)   # Остатки модели
jbplot_mcmc(fit)        # Диагностика сходимости MCMC
jbplot_trj(fit, type ) "B")      # Динамика биомассы
jbplot_trj(fit, type = "F")      # Динамика промысловой смертности
jbplot_trj(fit, type = "BBmsy")  # Отношение B/Bmsy
jbplot_trj(fit, type = "FFmsy")  # Отношение F/Fmsy
jbplot_spphase(fit)     # График продуктивности запаса с разметкой фаз Kobe
jbplot_kobe(fit)        # Коуб-график (B/Bmsy vs F/Fmsy)

# Дополнительные диагностики:
jbplot_runstest(fit)    # Тест на случайность остатков
jbplot_logfits(fit)     # Графики в логарифмической шкале
jbplot_procdev(fit)     # Отклонения процесса

# Сохранение временных рядов результатов
write.csv(fit$timeseries, file = "results.csv", row.names = FALSE)

# ------------------- 5. РЕТРОСПЕКТИВНЫЙ АНАЛИЗ --------------------

# Создание папки для результатов ретроспективы
retro.dir <- file.path(output.dir, "retro")
dir.create(retro.dir, showWarnings = FALSE)

# Запуск ретроспективного анализа (убираем по 1-5 лет)
hc <- hindcast_jabba(jbinput = jbinput, fit = fit, peels = 1:5)

# Визуализация ретроспективного анализа
mohnsrho <- jbplot_retro(
  hc, 
  as.png = FALSE,         # Чтобы сохранить как PNG-файл установите TRUE
  output.dir = retro.dir,
  xlim = c(2007, 2022)   # Ограничения по годам на графике
)

# Кросс-валидация
mase <- jbplot_hcxval(hc, as.png = FALSE, output.dir = retro.dir)

# ------------------- 6. ПРОГНОЗИРОВАНИЕ --------------------

# Прогноз на основе F (ловушечное усилие)
fw1 <- fw_jabba(
  fit,
  nyears = 10,       # Длина прогноза (лет)
  imp.yr = 1,        # Год внедрения новых правил
  imp.values = seq(10, 16, 2), # Варианты управления (уровни улова)
  quant = "Catch",   # Прогнозировать по уловам
  type = "abs",      # Абсолютные значения
  stochastic = TRUE  # Стохастический прогноз
)

# Графики ансамбля прогнозов
jbpar(mfrow = c(3, 2)) # Настройка макета графиков (3 строки, 2 столбца)
jbplot_ensemble(fw1)    # Основной график прогнозов

# График для B/Bmsy с кастомизацией
jbplot_ensemble(
  fw1,
  subplots = c(1),        # Только B/Bmsy
  add = TRUE,             # Добавить к текущему графику
  xlim = c(2020, 2035),   # Ограничение по годам
  legend.loc = "topleft"  # Позиция легенды
)

# ------------------- 7. ОБРАБОТКА ПРОГНОЗНЫХ ДАННЫХ --------------------

# Фильтрация данных прогноза (2025-2034) для выбранных сценариев
forecast_data <- subset(
  fw1, 
  year %in% 2025:2034 &    # Годы прогноза
    run %in% c("C10", "C12", "C14", "C16") & # Сценарии управления
    type == "prj"           # Только прогнозные значения
)

# Расчет медиан биомассы (B) по годам и сценариям
median_B <- aggregate(
  B ~ year + run,          # Формула: группировка по году и сценарию
  data = forecast_data, 
  FUN = median             # Функция агрегации
)

# Расчет медиан улова (Catch) по годам и сценариям
median_Catch <- aggregate(
  Catch ~ year + run, 
  data = forecast_data, 
  FUN = median
)

# Преобразование в широкий формат (годы по строкам, сценарии по столбцам)
b_table <- dcast(median_B, year ~ run, value.var = "B")
catch_table <- dcast(median_Catch, year ~ run, value.var = "Catch")

# Вывод таблиц
print("Медианная биомасса:")
print(b_table)

print("Медианные уловы:")
print(catch_table)

# Сохранение таблиц
write.csv(b_table, "biomass_forecast.csv", row.names = TRUE)
write.csv(catch_table, "catch_forecast.csv", row.names = TRUE)