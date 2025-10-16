# ==============================================================================
# НЕЧЁТКАЯ МОДЕЛЬ РЫБОПРОДУКТИВНОСТИ АЗОВСКОГО МОРЯ
# Основано на: Дяченко О.Ф.,2014,НЕЧЕТКАЯ МОДЕЛЬ ДЛЯ ОЦЕНКИ
# РЫБОПРОДУКТИВНОСТИ АЗОВСКОГО МОРЯ. Вестник ХНТУ, №1(48)
# 
# ОБУЧАЮЩИЙ МОДУЛЬ ДЛЯ ГИДРОБИОЛОГОВ:
# Этот скрипт демонстрирует применение нечёткой логики для оценки 
# рыбопродуктивности на основе двух ключевых факторов: количества планктона 
# и солёности воды. Модель позволяет работать как с числовыми данными, 
# так и с качественными описаниями.
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. ПОДГОТОВКА СРЕДЫ: УСТАНОВКА И ЗАГРУЗКА БИБЛИОТЕК
# ------------------------------------------------------------------------------
# Библиотеки для визуализации и работы с данными
if (!require("ggplot2")) install.packages("ggplot2")  # для создания графиков
if (!require("reshape2")) install.packages("reshape2") # для преобразования данных
if (!require("gridExtra")) install.packages("gridExtra") # для компоновки графиков

library(ggplot2)
library(reshape2)
library(gridExtra)

# ------------------------------------------------------------------------------
# 2. ОСНОВНОЙ МАТЕМАТИЧЕСКИЙ АППАРАТ: ФУНКЦИИ ПРИНАДЛЕЖНОСТИ
# ------------------------------------------------------------------------------

# ТРЕУГОЛЬНАЯ ФУНКЦИЯ ПРИНАДЛЕЖНОСТИ:
# Определяет, насколько значение x принадлежит к нечёткому множеству
# Параметры: a - начало, b - вершина, c - конец треугольника
triangular_mf <- function(x, a, b, c) {
  result <- numeric(length(x))
  for (i in 1:length(x)) {
    if (x[i] <= a | x[i] >= c) {
      result[i] <- 0
    } else if (x[i] > a & x[i] <= b) {
      if (b == a) result[i] <- 1
      else result[i] <- (x[i] - a) / (b - a)
    } else if (x[i] > b & x[i] < c) {
      if (c == b) result[i] <- 1
      else result[i] <- (c - x[i]) / (c - b)
    }
  }
  return(pmax(pmin(result, 1), 0))
}

# ФУНКЦИИ ПРИНАДЛЕЖНОСТИ ДЛЯ ПЛАНКТОНА (г/м?):
# Преобразует количественные данные в качественные оценки
plankton_mf <- function(x, term) {
  switch(term,
         "Мало" = triangular_mf(x, 0, 50, 100),      # 0-100 г/м? - мало
         "Средне" = triangular_mf(x, 50, 100, 150),  # 50-150 г/м? - средне  
         "Много" = triangular_mf(x, 100, 150, 200))  # 100-200 г/м? - много
}

# ФУНКЦИИ ПРИНАДЛЕЖНОСТИ ДЛЯ СОЛЁНОСТИ (промилле):
# Отражает толерантность рыб к солёности воды
salinity_mf <- function(x, term) {
  switch(term,
         "Низкая" = triangular_mf(x, 0, 2, 10),     # 0-10 ‰ - низкая (оптимально для азовских рыб)
         "Средняя" = triangular_mf(x, 5, 10, 15),   # 5-15 ‰ - средняя
         "Высокая" = triangular_mf(x, 10, 15, 20))  # 10-20 ‰ - высокая (неблагоприятно)
}

# ФУНКЦИИ ПРИНАДЛЕЖНОСТИ ДЛЯ ПРОДУКТИВНОСТИ (кг/га):
# Конечный результат - оценка рыбопродуктивности
productivity_mf <- function(x, term) {
  switch(term,
         "Очень_Низкая" = triangular_mf(x, 0, 10, 20),  # 0-20 кг/га
         "Низкая" = triangular_mf(x, 10, 20, 30),       # 10-30 кг/га
         "Средняя" = triangular_mf(x, 20, 35, 50),      # 20-50 кг/га
         "Высокая" = triangular_mf(x, 40, 50, 60),      # 40-60 кг/га
         "Очень_Высокая" = triangular_mf(x, 50, 60, 70)) # 50-70 кг/га
}

# ------------------------------------------------------------------------------
# 3. БАЗА ЗНАНИЙ: НЕЧЁТКИЕ ПРАВИЛА "ЕСЛИ-ТО"
# ------------------------------------------------------------------------------

# НЕЧЁТКИЕ ПРАВИЛА (алгоритм Мамдани):
# На основе биологических знаний о взаимосвязи факторов
fuzzy_rules_silent <- function(plankton, salinity) {
  # Вычисляем степени принадлежности для входных переменных
  plankton_low <- plankton_mf(plankton, "Мало")
  plankton_medium <- plankton_mf(plankton, "Средне")
  plankton_high <- plankton_mf(plankton, "Много")
  
  salinity_low <- salinity_mf(salinity, "Низкая")
  salinity_medium <- salinity_mf(salinity, "Средняя")
  salinity_high <- salinity_mf(salinity, "Высокая")
  
  # ПРАВИЛО 1: Много корма + Низкая солёность = Очень высокая продуктивность
  rule1 <- min(plankton_high, salinity_low)
  
  # ПРАВИЛО 2: Средний корм + Средняя солёность = Средняя продуктивность  
  rule2 <- min(plankton_medium, salinity_medium)
  
  # ПРАВИЛО 3: Мало корма + Высокая солёность = Очень низкая продуктивность
  rule3 <- min(plankton_low, salinity_high)
  
  return(c(rule1, rule2, rule3))
}

# ------------------------------------------------------------------------------
# 4. ДЕФАЗЗИФИКАЦИЯ: ПРЕОБРАЗОВАНИЕ НЕЧЁТКОГО РЕЗУЛЬТАТА В ЧИСЛО
# ------------------------------------------------------------------------------

# МЕТОД ЦЕНТРА ТЯЖЕСТИ (Center of Gravity):
# Находит "средневзвешенное" значение продуктивности
defuzzify <- function(activations) {
  productivity_values <- seq(0, 70, length.out = 100)
  
  # Агрегируем все активированные правила
  aggregated_mf <- pmax(
    pmin(activations[1], productivity_mf(productivity_values, "Очень_Высокая")),
    pmin(activations[2], productivity_mf(productivity_values, "Средняя")),
    pmin(activations[3], productivity_mf(productivity_values, "Очень_Низкая")),
    na.rm = TRUE
  )
  
  # Вычисляем центр тяжести полученной фигуры
  numerator <- sum(productivity_values * aggregated_mf, na.rm = TRUE)
  denominator <- sum(aggregated_mf, na.rm = TRUE)
  
  if (is.na(denominator) || denominator == 0) return(0)
  return(numerator / denominator)
}

# =============================================================================
# ВИЗУАЛИЗАЦИЯ 1: КАК МОДЕЛЬ "ВИДИТ" ВХОДНЫЕ ДАННЫЕ
# =============================================================================

# СОЗДАЁМ ДАННЫЕ ДЛЯ ВИЗУАЛИЗАЦИИ ФУНКЦИЙ ПРИНАДЛЕЖНОСТИ:

# Для солёности: диапазон 0-20 промилле
salinity_seq <- seq(0, 20, length.out = 100)
salinity_data <- data.frame(
  salinity = salinity_seq,
  Низкая = salinity_mf(salinity_seq, "Низкая"),
  Средняя = salinity_mf(salinity_seq, "Средняя"),
  Высокая = salinity_mf(salinity_seq, "Высокая")
)

# Для планктона: диапазон 0-200 г/м?
plankton_seq <- seq(0, 200, length.out = 100)
plankton_data <- data.frame(
  plankton = plankton_seq,
  Мало = plankton_mf(plankton_seq, "Мало"),
  Средне = plankton_mf(plankton_seq, "Средне"),
  Много = plankton_mf(plankton_seq, "Много")
)

# ПРЕОБРАЗУЕМ ДАННЫЕ ДЛЯ ПОСТРОЕНИЯ ГРАФИКОВ:
salinity_melted <- melt(salinity_data, id.vars = "salinity", 
                       variable.name = "Терм", value.name = "Принадлежность")
plankton_melted <- melt(plankton_data, id.vars = "plankton", 
                       variable.name = "Терм", value.name = "Принадлежность")

# СОЗДАЁМ СДВОЕННЫЙ ГРАФИК:
p1 <- ggplot(salinity_melted, aes(x = salinity, y = Принадлежность, color = Терм)) +
  geom_line(linewidth = 1.2) +
  labs(title = "Функции принадлежности для солёности воды",
       subtitle = "Как модель интерпретирует разные уровни солёности",
       x = "Солёность (промилле)", y = "Степень принадлежности") +
  theme_minimal() +
  scale_color_manual(values = c("blue", "green", "red"))

p2 <- ggplot(plankton_melted, aes(x = plankton, y = Принадлежность, color = Терм)) +
  geom_line(linewidth = 1.2) +
  labs(title = "Функции принадлежности для количества планктона", 
       subtitle = "Как модель интерпретирует кормовую базу",
       x = "Планктон (г/м?)", y = "Степень принадлежности") +
  theme_minimal() +
  scale_color_manual(values = c("red", "orange", "darkgreen"))

# ВЫВОДИМ ОБА ГРАФИКА РЯДОМ:
grid.arrange(p1, p2, ncol = 2, 
             top = "ФУНКЦИИ ПРИНАДЛЕЖНОСТИ: Как модель 'понимает' входные данные")

# =============================================================================
# ВИЗУАЛИЗАЦИЯ 2: КАК ФАКТОРЫ ВЛИЯЮТ НА РЕЗУЛЬТАТ
# =============================================================================

# СОЗДАЁМ ВСЕ ВОЗМОЖНЫЕ КОМБИНАЦИИ ПЛАНКТОНА И СОЛЁНОСТИ:
plankton_range <- seq(0, 200, length.out = 30)   # 30 значений планктона
salinity_range <- seq(0, 20, length.out = 30)    # 30 значений солёности

results <- matrix(0, nrow = length(plankton_range), ncol = length(salinity_range))

# РАСЧЁТ ПРОДУКТИВНОСТИ ДЛЯ КАЖДОЙ КОМБИНАЦИИ:
cat("Вычисление продуктивности для 900 комбинаций факторов...\n")
for (i in 1:length(plankton_range)) {
  for (j in 1:length(salinity_range)) {
    activations <- fuzzy_rules_silent(plankton_range[i], salinity_range[j])
    results[i, j] <- defuzzify(activations)
  }
}

# ПОДГОТАВЛИВАЕМ ДАННЫЕ ДЛЯ ТЕПЛОВОЙ КАРТЫ:
results_df <- expand.grid(Планктон = plankton_range, Солёность = salinity_range)
results_df$Продуктивность <- as.vector(results)

# СОЗДАЁМ ТЕПЛОВУЮ КАРТУ:
ggplot(results_df, aes(x = Планктон, y = Солёность, fill = Продуктивность)) +
  geom_tile() +
  scale_fill_gradient2(low = "red", mid = "yellow", high = "green", 
                      midpoint = 35, name = "Продуктивность\n(кг/га)") +
  labs(title = "ЗАВИСИМОСТЬ РЫБОПРОДУКТИВНОСТИ ОТ ПЛАНКТОНА И СОЛЁНОСТИ",
       subtitle = "Тёплые цвета = высокая продуктивность, холодные = низкая",
       x = "Планктон (г/м?)", y = "Солёность (промилле)") +
  theme_minimal()

# =============================================================================
# ПРИМЕР 1: ТОЧНЫЙ РАСЧЁТ ПО ЧИСЛОВЫМ ДАННЫМ
# =============================================================================

# ФУНКЦИЯ ДЛЯ ПОЛНОГО РАСЧЁТА:
calculate_productivity <- function(plankton_val, salinity_val) {
  # Применяем нечёткие правила
  activations <- fuzzy_rules_silent(plankton_val, salinity_val)
  productivity <- defuzzify(activations)
  
  # Вычисляем степени принадлежности для детального анализа
  plankton_deg <- c(
    Мало = plankton_mf(plankton_val, "Мало"),
    Средне = plankton_mf(plankton_val, "Средне"), 
    Много = plankton_mf(plankton_val, "Много")
  )
  
  salinity_deg <- c(
    Низкая = salinity_mf(salinity_val, "Низкая"),
    Средняя = salinity_mf(salinity_val, "Средняя"),
    Высокая = salinity_mf(salinity_val, "Высокая")
  )
  
  # Возвращаем все результаты
  list(
    plankton = plankton_val,
    salinity = salinity_val,
    productivity = productivity,
    plankton_degrees = plankton_deg,
    salinity_degrees = salinity_deg
  )
}

# ПРИМЕР 1: БЛАГОПРИЯТНЫЕ УСЛОВИЯ (много корма, низкая солёность)
example1 <- calculate_productivity(180, 8)

# ПРИМЕР 2: НЕБЛАГОПРИЯТНЫЕ УСЛОВИЯ (мало корма, высокая солёность)  
example2 <- calculate_productivity(30, 18)

# ВЫВОД РЕЗУЛЬТАТОВ:
cat("\n", strrep("=", 70), "\n")
cat("ПРИМЕР 1 - БЛАГОПРИЯТНЫЕ УСЛОВИЯ (оптимально для Азовского моря):\n")
cat("Планктон:", example1$plankton, "г/м?\n")
cat("Солёность:", example1$salinity, "промилле\n")
cat("Рассчитанная продуктивность:", round(example1$productivity, 2), "кг/га\n")
cat("Степени принадлежности для планктона:\n")
print(round(example1$plankton_degrees, 3))
cat("Степени принадлежности для солёности:\n") 
print(round(example1$salinity_degrees, 3))

cat("\n", strrep("-", 70), "\n")

cat("ПРИМЕР 2 - НЕБЛАГОПРИЯТНЫЕ УСЛОВИЯ (антропогенное воздействие):\n")
cat("Планктон:", example2$plankton, "г/м?\n")
cat("Солёность:", example2$salinity, "промилле\n")
cat("Рассчитанная продуктивность:", round(example2$productivity, 2), "кг/га\n")
cat("Степени принадлежности для планктона:\n")
print(round(example2$plankton_degrees, 3))
cat("Степени принадлежности для солёности:\n")
print(round(example2$salinity_degrees, 3))

# =============================================================================
# ПРИМЕР 2: БЫСТРАЯ ОЦЕНКА ПО КАЧЕСТВЕННЫМ ОПИСАНИЯМ
# =============================================================================

# ФУНКЦИЯ ДЛЯ ЛИНГВИСТИЧЕСКОГО ВВОДА (удобно при отсутствии точных данных):
estimate_by_linguistic <- function(plankton_term = "Средне", salinity_term = "Средняя") {
  # Проверяем корректность введённых терминов
  valid_plankton <- c("Мало", "Средне", "Много")
  valid_salinity <- c("Низкая", "Средняя", "Высокая")
  
  if (!(plankton_term %in% valid_plankton))
    stop("ОШИБКА: Неверный терм для планктона. Используйте: 'Мало', 'Средне', 'Много'")
  if (!(salinity_term %in% valid_salinity))
    stop("ОШИБКА: Неверный терм для солёности. Используйте: 'Низкая', 'Средняя', 'Высокая'")
  
  # Преобразуем слова в характерные числовые значения
  plankton_val <- switch(plankton_term,
    "Мало" = 50,    # характерное значение для "Мало"
    "Средне" = 100, # характерное значение для "Средне"
    "Много" = 150   # характерное значение для "Много"
  )
  
  salinity_val <- switch(salinity_term,
    "Низкая" = 5,   # характерное значение для "Низкая"  
    "Средняя" = 10, # характерное значение для "Средняя"
    "Высокая" = 15  # характерное значение для "Высокая"
  )
  
  # Вычисляем продуктивность
  res <- calculate_productivity(plankton_val, salinity_val)
  
  # Выводим понятный результат
  cat("\n", strrep("=", 70), "\n")
  cat("ЛИНГВИСТИЧЕСКАЯ ОЦЕНКА (быстрая, по качественным описаниям):\n")
  cat("Входные данные:\n")
  cat("  • Планктон:", plankton_term, "\n")
  cat("  • Солёность:", salinity_term, "\n")
  cat("Прогнозируемая продуктивность:", round(res$productivity, 2), "кг/га\n")
  
  # Возвращаем результат для дальнейшего использования
  invisible(res)
}

# ДЕМОНСТРАЦИЯ ЛИНГВИСТИЧЕСКОГО ВВОДА:
cat("\n", strrep("=", 70), "\n")
cat("ДЕМОНСТРАЦИЯ РАБОТЫ С КАЧЕСТВЕННЫМИ ОПИСАНИЯМИ:\n")
cat("(используйте эти функции когда нет точных измерений)\n\n")

# Пример A: Оптимальные условия
estimate_by_linguistic("Много", "Низкая")

# Пример B: Критические условия  
estimate_by_linguistic("Мало", "Высокая")

# Пример C: Средние условия
estimate_by_linguistic("Средне", "Средняя")

cat("\n", strrep("=", 70), "\n")
cat("ОБУЧАЮЩИЙ МОДУЛЬ ЗАВЕРШЁН!\n")
cat("Вы научились:\n")
cat("• Понимать принципы нечёткой логики в гидробиологии\n")  
cat("• Использовать модель для численных расчётов\n")
cat("• Проводить быстрые оценки по качественным описаниям\n")
cat("• Интерпретировать визуализации взаимосвязей факторов\n")
cat(strrep("=", 70), "\n")