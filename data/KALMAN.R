# ПРАКТИКУМ: ФИЛЬТР КАЛМАНА ДЛЯ ОЦЕНКИ ЗАПАСА ГИДРОБИОНТОВ
# =====================================================================================================================
# Автор: Баканев С.В. Для курса "Оценка водных биоресурсов при недостатке данных в среде R (для начинающих)"
# Цель: Показать преимущества фильтра Калмана перед классическими методами
# =====================================================================================================================

library(ggplot2)
library(dplyr)
library(tidyr)

# ЧАСТЬ 1: ГЕНЕРАЦИЯ И ВИЗУАЛИЗАЦИЯ СИНТЕТИЧЕСКИХ ДАННЫХ
# ======================================================
# В реальных условиях мы не знаем истинного состояния запаса,
# но для демонстрации методов создадим синтетические данные

cat("=== ЧАСТЬ 1: ГЕНЕРАЦИЯ И ВИЗУАЛИЗАЦИЯ ДАННЫХ ===\n\n")

# Устанавливаем seed для воспроизводимости результатов
set.seed(123)
days <- 1:20  # 20 дней наблюдений

# СОЗДАЕМ РЕАЛИСТИЧНУЮ ДИНАМИКУ ЗАПАСА:
# 1. Экспоненциальное снижение (естественная смертность + вылов)
# 2. Миграционные всплески (приток/отток особей)
# 3. Случайные колебания (неучтенные факторы)

cat("Создание модели истинного запаса...\n")
migration_pattern <- c(rep(0,5), 50, 80, 30, rep(0,5), -40, -60, -20, rep(0,3))
true_stock <- 1000 * exp(-0.08 * days) + migration_pattern + rnorm(20, sd = 25)

cat("Генерация данных промысла...\n")
# Создаем dataframe с данными наблюдений
data <- data.frame(
  Day = days,
  # Суточный вылов (2-5% от текущего запаса)
  Catch = pmax(round(true_stock * runif(20, 0.02, 0.05)), 0),
  # CPUE (Catch Per Unit Effort) - основной наблюдаемый показатель
  # В реальности зависит от множества факторов, здесь - линейная зависимость + шум
  CPUE = true_stock * 0.002 + rnorm(20, sd = 0.1)
)

# БАЗОВАЯ СТАТИСТИКА ДАННЫХ
cat("Средний улов на усилие (CPUE):", round(mean(data$CPUE), 3), "\n")
cat("Общий вылов за период:", sum(data$Catch), "условных единиц\n")
cat("Максимальный суточный вылов:", max(data$Catch), "\n")
cat("Минимальный суточный вылов:", min(data$Catch), "\n\n")

# ВИЗУАЛИЗАЦИЯ ИСХОДНЫХ ДАННЫХ
cat("Построение графика исходных данных...\n")
plot_raw_data <- ggplot(data, aes(x = Day)) +
  # Столбцы - суточный вылов (масштабируем для совмещения с CPUE)
  geom_col(aes(y = Catch/50), fill = "lightblue", alpha = 0.7, width = 0.7) +
  # Точки и линия - наблюдения CPUE
  geom_point(aes(y = CPUE), color = "darkred", size = 3) +
  geom_line(aes(y = CPUE), color = "darkred", linetype = "dashed", alpha = 0.7) +
  # Две оси Y для разных масштабов
  scale_y_continuous(
    name = "Улов на усилие (CPUE)",
    sec.axis = sec_axis(~ . * 50, name = "Суточный улов (особи)")
  ) +
  labs(title = "Исходные данные промысла",
       subtitle = "Красные точки/линия - CPUE (индекс обилия), синие столбцы - уловы",
       caption = "CPUE используется как прокси-переменная для оценки запаса") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

print(plot_raw_data)

# ЧАСТЬ 2: КЛАССИЧЕСКАЯ МОДЕЛЬ ЛЕСЛИ
# ===================================
# Модель Лесли - традиционный подход для оценки динамики популяций
# Предполагает линейное изменение CPUE во времени

cat("\n=== ЧАСТЬ 2: КЛАССИЧЕСКАЯ МОДЕЛЬ ЛЕСЛИ ===\n\n")

# Линейная регрессия CPUE от времени
cat("Построение линейной регрессии CPUE ~ Day...\n")
lm_model <- lm(CPUE ~ Day, data = data)

# Добавляем предсказания модели в dataframe
data$LM_Prediction <- predict(lm_model)

# ВЫВОД РЕЗУЛЬТАТОВ РЕГРЕССИИ
cat("РЕЗУЛЬТАТЫ ЛИНЕЙНОЙ РЕГРЕССИИ:\n")
cat("Наклон (темп снижения CPUE):", round(coef(lm_model)[2], 5), "\n")
cat("R? (доля объясненной дисперсии):", round(summary(lm_model)$r.squared, 3), "\n")
cat("p-значение:", round(summary(lm_model)$coefficients[2,4], 5), "\n\n")

# ВИЗУАЛИЗАЦИЯ МОДЕЛИ ЛЕСЛИ
cat("Визуализация результатов линейной модели...\n")
plot_leslie <- ggplot(data, aes(x = Day, y = CPUE)) +
  geom_point(color = "darkred", size = 3, alpha = 0.7) +
  geom_line(aes(y = LM_Prediction), color = "darkblue", size = 1.2, 
            linetype = "solid") +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.2, color = "blue", 
              fill = "lightblue") +
  labs(title = "Классическая модель Лесли: линейная регрессия CPUE",
       subtitle = paste("Наклон =", round(coef(lm_model)[2], 4), 
                       ", R? =", round(summary(lm_model)$r.squared, 3)),
       x = "День промысла", 
       y = "Улов на усилие (CPUE)") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

print(plot_leslie)

# ЧАСТЬ 3: ФИЛЬТР КАЛМАНА ДЛЯ ОЦЕНКИ ЗАПАСА И МИГРАЦИИ
# ====================================================
# Фильтр Калмана позволяет оценивать ненаблюдаемые состояния системы
# (запас, миграцию) на основе зашумленных наблюдений (CPUE)

cat("\n=== ЧАСТЬ 3: ФИЛЬТР КАЛМАНА ===\n\n")

# ПАРАМЕТРЫ ФИЛЬТРА КАЛМАНА
cat("Инициализация параметров фильтра Калмана...\n")
q <- 0.002      # Коэффициент улавливаемости (связь между запасом и CPUE)
sigma_N <- 30   # Шум процесса для запаса (неопределенность в динамике запаса)
sigma_m <- 20   # Шум процесса для миграции (неопределенность в миграции)
sigma_y <- 0.08 # Шум измерений для CPUE (ошибка наблюдений)

# НАЧАЛЬНЫЕ УСЛОВИЯ
N_est <- data$CPUE[1] / q  # Начальная оценка запаса из первого наблюдения CPUE
m_est <- 0                 # Начальная оценка миграции (предполагаем 0)
P <- matrix(c(100, 0, 0, 100), nrow = 2)  # Начальная матрица ковариации ошибок

# ДЛЯ ХРАНЕНИЯ РЕЗУЛЬТАТОВ
results <- data.frame(
  Day = integer(),
  CPUE_Observed = numeric(),
  Stock_Estimate = numeric(),
  Migration_Estimate = numeric(),
  CPUE_Estimated = numeric()
)

cat("Запуск основного цикла фильтра Калмана...\n")

# ОСНОВНОЙ ЦИКЛ ФИЛЬТРА КАЛМАНА (для каждого дня наблюдений)
for (t in 1:nrow(data)) {
  
  # ----- ШАГ 1: ПРОГНОЗ (PREDICTION STEP) -----
  # Прогнозируем состояние системы на следующий день
  
  # Уравнения состояния:
  N_pred <- N_est - data$Catch[t] + m_est  # Запас: предыдущий - вылов + миграция
  m_pred <- 0.9 * m_est                    # Миграция: авторегрессия 1-го порядка
  
  # Матрица перехода состояний (описывает динамику системы)
  F_matrix <- matrix(c(1, 1,     # Влияние миграции на запас
                       0, 0.9),  # Авторегрессия для миграции
                     nrow = 2, byrow = TRUE)
  
  # Матрица ковариации шума процесса
  Q <- matrix(c(sigma_N^2, 0, 
                0, sigma_m^2), nrow = 2)
  
  # Прогноз ковариации ошибок
  P_pred <- F_matrix %*% P %*% t(F_matrix) + Q
  
  # ----- ШАГ 2: КОРРЕКЦИЯ (UPDATE STEP) -----
  # Обновляем прогноз на основе новых наблюдений
  
  # Матрица измерений (связь состояний с наблюдениями)
  H <- matrix(c(q, 0), nrow = 1)  # Только запас влияет на CPUE
  
  # Прогнозируемое значение CPUE
  y_pred <- q * N_pred
  
  # Инновация (невязка) - разница между наблюдением и прогнозом
  innovation <- data$CPUE[t] - y_pred
  
  # Ковариация инноваций
  S <- H %*% P_pred %*% t(H) + sigma_y^2
  
  # Коэффициент усиления Калмана (оптимальные веса для обновления)
  K <- P_pred %*% t(H) %*% solve(S)
  
  # ОБНОВЛЕНИЕ ОЦЕНОК СОСТОЯНИЙ
  state_update <- K %*% innovation
  N_est <- N_pred + state_update[1,1]  # Обновляем оценку запаса
  m_est <- m_pred + state_update[2,1]  # Обновляем оценку миграции
  
  # ОБНОВЛЕНИЕ КОВАРИАЦИИ ОШИБОК
  P <- (diag(2) - K %*% H) %*% P_pred
  
  # СОХРАНЕНИЕ РЕЗУЛЬТАТОВ ДЛЯ ЭТОГО ДНЯ
  results <- rbind(results, data.frame(
    Day = t,
    CPUE_Observed = data$CPUE[t],
    Stock_Estimate = N_est,
    Migration_Estimate = m_est,
    CPUE_Estimated = y_pred
  ))
}

cat("Фильтр Калмана успешно завершил работу!\n\n")

# ОБЪЕДИНЕНИЕ ВСЕХ ДАННЫХ ДЛЯ АНАЛИЗА
final_data <- cbind(data, results[, c("Stock_Estimate", "Migration_Estimate", "CPUE_Estimated")])

# ЧАСТЬ 4: ЗАКЛЮЧИТЕЛЬНАЯ АНАЛИТИКА И ВИЗУАЛИЗАЦИЯ
# ================================================

cat("=== ЧАСТЬ 4: СРАВНИТЕЛЬНЫЙ АНАЛИЗ РЕЗУЛЬТАТОВ ===\n\n")

# ОСНОВНЫЕ РЕЗУЛЬТАТЫ ФИЛЬТРА КАЛМАНА
cat("ОЦЕНКИ ФИЛЬТРА КАЛМАНА:\n")
cat("• Начальный запас:", round(results$Stock_Estimate[1], 1), "особей\n")
cat("• Конечный запас:", round(results$Stock_Estimate[20], 1), "особей\n")
cat("• Снижение запаса за период:", 
    round(results$Stock_Estimate[1] - results$Stock_Estimate[20], 1), 
    "особей (", 
    round((results$Stock_Estimate[1] - results$Stock_Estimate[20]) / results$Stock_Estimate[1] * 100, 1),
    "%)\n")
cat("• Средняя миграция:", round(mean(results$Migration_Estimate), 1), "особей/день\n")
cat("• Максимальный приток:", round(max(results$Migration_Estimate), 1), "особей/день\n")
cat("• Максимальный отток:", round(min(results$Migration_Estimate), 1), "особей/день\n\n")

# ГРАФИК 1: СРАВНЕНИЕ МЕТОДОВ ОЦЕНКИ CPUE
cat("Построение сравнительных графиков...\n")
plot_comparison <- ggplot(final_data) +
  # Наблюдения
  geom_point(aes(x = Day, y = CPUE, color = "Наблюдения CPUE"), size = 3, alpha = 0.7) +
  # Модель Лесли
  geom_line(aes(x = Day, y = LM_Prediction, color = "Модель Лесли"), 
            linetype = "dashed", size = 1.2, alpha = 0.8) +
  # Фильтр Калмана
  geom_line(aes(x = Day, y = CPUE_Estimated, color = "Фильтр Калмана"), 
            size = 1.5, alpha = 0.9) +
  # Настройки графика
  labs(title = "Сравнение методов оценки динамики запаса",
       subtitle = "Фильтр Калмана vs Классическая модель Лесли",
       x = "День промысла", 
       y = "Улов на усилие (CPUE)",
       color = "Метод") +
  scale_color_manual(values = c(
    "Наблюдения CPUE" = "darkred",
    "Модель Лесли" = "darkblue", 
    "Фильтр Калмана" = "darkgreen"
  )) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        legend.position = "bottom")

print(plot_comparison)

# ГРАФИК 2: ДИНАМИКА ЗАПАСА И МИГРАЦИИ
plot_stock_migration <- ggplot(final_data) +
  # Оценка запаса
  geom_line(aes(x = Day, y = Stock_Estimate, color = "Оценка запаса"), 
            size = 1.3) +
  # Миграционный поток (масштабируем для визуализации)
  geom_line(aes(x = Day, y = Migration_Estimate * 10 + 400, 
                color = "Миграция (x10)"), size = 1.1, alpha = 0.8) +
  # Две оси Y
  scale_y_continuous(
    name = "Оценка запаса (особи)",
    sec.axis = sec_axis(~ (. - 400) / 10, name = "Миграционный поток (особи/день)")
  ) +
  labs(title = "Динамика запаса и миграционного потока",
       subtitle = "Результаты работы фильтра Калмана",
       x = "День промысла",
       color = "Переменные") +
  scale_color_manual(values = c(
    "Оценка запаса" = "steelblue",
    "Миграция (x10)" = "orange"
  )) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        legend.position = "bottom")

print(plot_stock_migration)

# ГРАФИК 3: ДЕТАЛЬНЫЙ АНАЛИЗ МИГРАЦИИ
plot_migration_detail <- ggplot(final_data, aes(x = Day, y = Migration_Estimate)) +
  geom_line(color = "firebrick", size = 1.3) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5, size = 1) +
  geom_point(aes(color = ifelse(Migration_Estimate > 0, "Иммиграция", "Эмиграция")), 
             size = 3) +
  labs(title = "Оценка миграционного потока фильтром Калмана",
       subtitle = "Положительные значения - приток, отрицательные - отток",
       x = "День промысла", 
       y = "Чистый миграционный поток (особи/день)",
       color = "Тип миграции") +
  scale_color_manual(values = c("Иммиграция" = "green4", "Эмиграция" = "red3")) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        legend.position = "bottom")

print(plot_migration_detail)

# ТАБЛИЦА РЕЗУЛЬТАТОВ
cat("=== ТАБЛИЦА РЕЗУЛЬТАТОВ (первые 10 дней) ===\n")
results_table <- head(final_data[, c("Day", "Catch", "CPUE", "Stock_Estimate", "Migration_Estimate")], 10) %>%
  mutate(across(where(is.numeric), round, 1))

print(results_table)

# СТАТИСТИКА КАЧЕСТВА ОЦЕНОК
cat("\n=== СТАТИСТИКА КАЧЕСТВА ОЦЕНОК ===\n")

# Среднеквадратичные ошибки
mse_leslie <- mean((final_data$CPUE - final_data$LM_Prediction)^2)
mse_kalman <- mean((final_data$CPUE - final_data$CPUE_Estimated)^2)

cat("СРЕДНЕКВАДРАТИЧНЫЕ ОШИБКИ (MSE):\n")
cat("• Модель Лесли:", round(mse_leslie, 4), "\n")
cat("• Фильтр Калмана:", round(mse_kalman, 4), "\n")
cat("• Улучшение:", round((mse_leslie - mse_kalman) / mse_leslie * 100, 1), "%\n\n")

# ЗАКЛЮЧЕНИЕ
cat("=== ВЫВОДЫ И РЕКОМЕНДАЦИИ ===\n")
cat("1. Фильтр Калмана позволяет оценивать ненаблюдаемые параметры (запас, миграцию)\n")
cat("2. Учитывает не только наблюдения, но и знания о динамике системы\n")
cat("3. Обеспечивает более точные оценки по сравнению с классическими методами\n")
cat("4. Особенно полезен при наличии миграционных процессов\n")
cat("5. Требует тщательной настройки параметров шумов\n\n")

cat("Анализ завершен! Для углубленного изучения рекомендуется:\n")
cat("• Экспериментировать с параметрами шумов (sigma_N, sigma_m, sigma_y)\n")
cat("• Добавить сезонные компоненты в модель\n")
cat("• Рассмотреть нелинейные версии фильтра Калмана\n")