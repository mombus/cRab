# ==============================================================================
# ПРАКТИЧЕСКОЕ ЗАНЯТИЕ: АНАЛИЗ ФАКТОРОВ И ПРОГНОЗ ПОПОЛНЕНИЯ ЗАПАСА
# Курс: "Оценка водных биоресурсов в среде R (для начинающих)"
# Структура:
# 1) Подготовка данных и выюор предикторов
# 2) Базовое сравнение моделей (5-fold CV + holdout)
# 3) Выбор лучшей прогностической модели (time-slice CV на 3 года + хронологический тест)
# 4) Прогноз 2022–2024 (ансамбль CUBIST+LM) и график 1990–2024 с ДИ
# ------------------------------------------------------------------------------
# Пояснения к занятию (для начинающих):
# - Мы работаем с временным рядом пополнения запаса R3haddock и набором факторов
#   среды/биомассы. Цель — построить понятные и проверяемые модели прогноза.
# - Сначала отберём информативные предикторы (Boruta и LASSO), затем сравним
#   разные модели машинного обучения на кросс-валидации (CV), после чего выберем
#   лучшую схему по time-slice CV (учитывая хронологию), и сделаем прогноз.
# - Ничего в коде ниже не изменено по сути — добавлены только комментарии,
#   логические подписи и подсказки для понимания каждого шага.
# ==============================================================================


# ==============================================================================
# 1) ВЫБОР ПРЕДИКТОРОВ
# ------------------------------------------------------------------------------
# Цель блока: привести данные к числовому виду, обработать пропуски, сократить
# мультиколлинеарность (сильные корреляции), а затем автоматически выделить
# кандидатов-предикторов двумя методами (Boruta, LASSO). В конце сформируем
# финальный пул признаков и проверим их значимость в простой LM.
# ==============================================================================

# Установка и подключение необходимых библиотек
# Для автоматического отбора предикторов нам понадобятся дополнительные пакеты
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  readxl, tidyverse, caret, corrplot, mgcv, randomForest, xgboost,
  Boruta,GGally, FactoMineR, glmnet, recipes, rsample  # Новые библиотеки для автоматического отбора
)

# Очистка среды и установка рабочей директории
# Совет: rm(list=ls()) очищает все объекты в памяти R; setwd задаёт папку,
# где искать/сохранять файлы. Убедитесь, что путь корректен на вашей машине.
rm(list = ls())
setwd("C:/RECRUITMENT/")

# Пакеты для расширенного отбора предикторов
# Boruta — обёртка над Random Forest для отбора признаков;
# glmnet — регуляризация (LASSO/ElasticNet) для отбора/усиления обобщающей способности;
# FactoMineR — PCA и другие многомерные методы (используем как утилиту).
library(Boruta)   # Алгоритм обертки для отбора признаков
library(glmnet)   # LASSO-регрессия
library(FactoMineR) # PCA анализ


# Загрузка и первичная обработка данных
# Шаги: фильтруем годы, приводим типы к числовому, заменяем строковые "NA" на NA.
DATA <- readxl::read_excel("RECRUITMENT.xlsx", sheet = "RECRUITMENT") %>%
  filter(YEAR > 1989 & YEAR < 2022) %>%
  # Преобразуем необходимые столбцы в числовой формат
  mutate(
    across(starts_with("T"), as.numeric),
    across(starts_with("I"), as.numeric),
    across(starts_with("O"), as.numeric),
  ) %>%
  # Обработка пропущенных значений (заменяем строку "NA" на NA)
  mutate(across(where(is.character), ~na_if(., "NA")))

# 1. Подготовка данных -------------------------------------------------------
# Выделим все возможные предикторы, включая географию и индексы трески
# Примечание: оставляем только числовые переменные, т.к. большинство моделей
# требует числовой вход без категориальных уровней.
predictors <- DATA %>% 
  select(-YEAR, -R3haddock) %>% 
  select_if(is.numeric) # Только числовые переменные

# Целевая переменная
response <- DATA$R3haddock

# 2. Обработка пропусков -----------------------------------------------------
# Заполнение медианными значениями — простой и устойчивый способ справиться с NA.
# Альтернативы: множественная иммутация (mice), KNN-impute и др.
predictors_filled <- predictors %>%
  mutate(across(everything(), ~ifelse(is.na(.), median(., na.rm = TRUE), .)))

# 3. Предварительный анализ корреляций ---------------------------------------
# Зачем: высокие корреляции затрудняют интерпретацию и могут вредить ряду моделей.
cor_matrix <- cor(predictors_filled, use = "complete.obs")
corrplot(cor_matrix, method = "circle", type = "upper", tl.cex = 0.7)

# Удаляем высокоскоррелированные предикторы (r > 0.8)
# Это механическое сокращение мультиколлинеарности до этапа отбора.
high_cor <- findCorrelation(cor_matrix, cutoff = 0.8)
predictors_filtered <- predictors_filled[, -high_cor]

# 4. Автоматизированный отбор Boruta (обертка Random Forest) -----------------
# Идея: определить признаки, которые важнее, чем случайный шум (shadow features).
set.seed(123)
boruta_output <- Boruta(
  x = predictors_filtered, 
  y = response,
  maxRuns = 100,
  doTrace = 2
)

# Визуализация результатов
plot(boruta_output, cex.axis = 0.7, las = 2)
boruta_stats <- attStats(boruta_output)
selected_vars <- getSelectedAttributes(boruta_output, withTentative = TRUE)

# 5. LASSO с более строгим критерием ------------------------------------------
# Идея: L1-регуляризация зануляет коэффициенты «слабых» предикторов.
# Выбор lambda.1se вместо lambda.min — более консервативный (простая модель).
x <- as.matrix(predictors_filtered)
y <- response

# Кросс-валидация
cv_fit <- cv.glmnet(x, y, alpha = 1, nfolds = 10)
plot(cv_fit)

# ?? ИСПОЛЬЗУЕМ lambda.1se вместо lambda.min — СТРОЖЕ!
lasso_coef <- coef(cv_fit, s = "lambda.1se")  # <-- Ключевое изменение!
lasso_vars <- rownames(lasso_coef)[lasso_coef[,1] != 0][-1]  # исключаем (Intercept)


# 6. Сравнение отобранных предикторов ----------------------------------------
# Полезно видеть, какие признаки отмечают оба метода (устойчивые кандидаты).
cat("Boruta selected:", length(selected_vars), "variables\n")
print(selected_vars)

cat("\nLASSO selected:", length(lasso_vars), "variables\n")
print(lasso_vars)

# 7. Финальный набор предикторов (объединение результатов) -------------------
# Логика: объединяем списки, добавляем биологически важные переменные вручную.
final_vars <- union(selected_vars, lasso_vars) 

# Добавляем обязательные переменные по биологической логике
mandatory <- c("haddock68")
final_vars <- union(final_vars, mandatory) %>% unique()

# 8. Проверка значимости -----------------------------------------------------
# Быстрая оценка значимости с LM: не как окончательный вывод, а как sanity-check.
final_model <- lm(response ~ as.matrix(predictors_filled[, final_vars]))
summary(final_model)

# 9. Формирование финального датасета ----------------------------------------
# Собираем набор с откликом и выбранными предикторами; удалим строки с NA.
model_data <- DATA %>%
  select(R3haddock, all_of(final_vars)) %>%
  drop_na()

# Просмотр структуры финальных данных
glimpse(model_data)

# Визуализация важности переменных
# Внимание: важности от RF — относительные; сопоставляйте с предметной логикой.
var_importance <- randomForest(R3haddock ~ ., data = model_data, importance = TRUE)
varImpPlot(var_importance, main = "Важность предикторов")


# ==============================================================================
#  ПОДГОТОВКА ДАННЫХ
# Создаём NAOspring, фиксируем финальный набор признаков, сохраняем CSV.
# ------------------------------------------------------------------------------
# Цель блока: стандартизировать набор признаков для дальнейшего сравнения
# моделей и обеспечить воспроизводимость (фиксированный CSV с нужными полями).
# ==============================================================================

# 1.1 Пакеты и окружение
# Примечание: блок повторяет базовую инициализацию для автономного запуска.
if (!require("pacman")) install.packages("pacman")
pacman::p_load(readxl, tidyverse, caret, corrplot)

rm(list = ls())
set.seed(123)
setwd("C:/RECRUITMENT/")

# 1.2 Загрузка исходных данных и приведение типов
DATA <- readxl::read_excel("RECRUITMENT.xlsx", sheet = "RECRUITMENT") %>%
  filter(YEAR > 1989 & YEAR < 2022) %>%
  mutate(
    across(starts_with("T"), as.numeric),
    across(starts_with("I"), as.numeric),
    across(starts_with("O"), as.numeric),
    across(where(is.character), ~na_if(., "NA"))
  )

# 1.3 Создаём NAOspring (если есть NAO3, NAO4, NAO5)
# Идея: агрегируем весенний индекс NAO как среднее за месяцы 3–5.
if (all(c("NAO3","NAO4","NAO5") %in% names(DATA))) {
  DATA <- DATA %>%
    mutate(NAOspring = rowMeans(pick(NAO3, NAO4, NAO5), na.rm = TRUE)) %>%
    select(-NAO3, -NAO4, -NAO5)
}

# 1.4 Финальный учебный набор предикторов (фиксируем)
# Важно: проверяем присутствие нужных колонок и формируем компактный датасет.
needed <- c("codTSB", "T12", "I5", "NAOspring", "haddock68")
stopifnot(all(needed %in% names(DATA)))

# Сохраняем YEAR в CSV (ниже он будет отброшен при обучении, но нужен для графика)
model_data <- DATA %>%
  select(YEAR, all_of(needed), R3haddock) %>%
  drop_na()

write.csv(model_data, "selected_predictors_dataset.csv", row.names = FALSE)
glimpse(model_data)

# (необязательно) Глянуть попарные связи и корреляции
# ggpairs может быть медленным, оставим по желанию
 ggpairs(model_data, columns = 2:7,
         lower = list(continuous = wrap("smooth", alpha = 0.3, size = 0.5)),
         upper = list(cor = wrap("cor", size = 3)))



# ==============================================================================
# 2) БАЗОВОЕ СРАВНЕНИЕ МОДЕЛЕЙ (5-FOLD CV + HOLDOUT)
# Единые фолды CV, тренировочно-тестовое разбиение, сводка метрик.
# ------------------------------------------------------------------------------
# Идея блока: быстрая «панель» сравнения разных семейств моделей на одинаковых
# условиях (одинаковые фолды CV) и внешний тест (holdout). Это помогает увидеть
# уровни ошибок и выбрать несколько лидеров для более строгой проверки далее.
# ==============================================================================

# 2.1 Пакеты и данные
pacman::p_load(mgcv, randomForest, xgboost, nnet, earth, kernlab, pls, Cubist, ranger, gbm, lattice)

model_data <- read.csv("selected_predictors_dataset.csv", header = TRUE, stringsAsFactors = FALSE)
# Если YEAR отсутствует (на всякий случай), создадим
if (!"YEAR" %in% names(model_data)) {
  model_data$YEAR <- seq(1990, by = 1, length.out = nrow(model_data))
}

# Используем только предикторы и отклик (YEAR исключаем)
model_data <- model_data %>%
  select(codTSB, T12, I5, NAOspring, haddock68, R3haddock) %>%
  na.omit()

# 2.2 Holdout и CV-контроллер
# Пропорция 80/20 обеспечивает внешний тест; внутри train — 5-fold CV для
# корректной настройки моделей и оценки средней ошибки.
train_idx <- caret::createDataPartition(model_data$R3haddock, p = 0.8, list = FALSE)
train <- model_data[train_idx, ]
test  <- model_data[-train_idx, ]

ctrl <- caret::trainControl(method = "cv", number = 5, savePredictions = "final")

# 2.3 Кастомный GAM (mgcv) для caret (bs="tp", REML, select=TRUE)
# GAM даёт гладкие нелинейности по каждому признаку; REML стабилизирует оценку.
gam_spec <- list(
  type = "Regression", library = "mgcv", loop = NULL,
  parameters = data.frame(parameter = "none", class = "character", label = "none"),
  grid = function(x,y,len=NULL,search="grid") data.frame(none = NA),
  fit = function(x,y,...) {
    df <- x; df$R3haddock <- y
    mgcv::gam(
      R3haddock ~ s(codTSB,bs="tp") + s(T12,bs="tp") + s(I5,bs="tp") +
                  s(NAOspring,bs="tp") + s(haddock68,bs="tp"),
      data=df, method="REML", select=TRUE, ...
    )
  },
  predict = function(modelFit, newdata, submodels = NULL) {
    predict(modelFit, newdata = newdata, type = "response")
  },
  prob = NULL, sort = function(x) x
)

# 2.4 Обучение моделей
# Подсказка: разные методы по-разному чувствительны к масштабу, числу признаков
# и мультиколлинеарности. Мы применяем одинаковые фолды CV для честного сравнения.
# --- 4. Линейная регрессия (LM)
# Учебный смысл: базовая линейная модель; ориентир для сравнения.
lm_model    <- caret::train(R3haddock ~ ., data = train, method = "lm", trControl = ctrl)
# --- 5. Обобщённая линейная модель (GLM: Gamma с лог-ссылкой)
# Учебный смысл: модель для положительных откликов; допускает нелинейность в шкале log.
glm_model   <- caret::train(R3haddock ~ ., data = train, method = "glm",
                            family = Gamma(link = "log"), trControl = ctrl)
# --- 6. Обобщённая аддитивная модель (GAM, mgcv: bs="tp", REML, select=TRUE)
# Учебный смысл: гибкие гладкие нелинейности по каждому предиктору.
gam_model   <- caret::train(x = train[, -which(names(train)=="R3haddock")],
                            y = train$R3haddock, method = gam_spec, trControl = ctrl)
# --- 7. Random Forest (rf: ntree=1000, mtry=1)
# Учебный смысл: ансамбль деревьев; устойчив к шуму; нелинейности/взаимодействия “из коробки”.
rf_model    <- caret::train(R3haddock ~ ., data = train, method = "rf", trControl = ctrl,
                            ntree = 1000, tuneGrid = data.frame(mtry = 1), importance = TRUE)
# --- 8. XGBoost (xgbTree) — точно как в “первом варианте”
# Учебный смысл: бустинг деревьев; сильная ML?модель, легко переобучается без валидации.
xgb_grid    <- expand.grid(nrounds=100, max_depth=4, eta=0.1, gamma=0,
                           colsample_bytree=0.8, min_child_weight=1, subsample=0.8)
xgb_model   <- caret::train(R3haddock ~ ., data = train, method = "xgbTree",
                            trControl = ctrl, tuneGrid = xgb_grid, verbose = 0)
# --- 9. Нейросеть (MLP, nnet: линейный выход, стандартизация)
# Учебный смысл: универсальный аппроксиматор; чувствителен к масштабу; требует регуляризации.
nnet_model  <- caret::train(R3haddock ~ ., data = train, method = "nnet",
                            trControl = ctrl, preProcess = c("center","scale"),
                            tuneGrid = expand.grid(size = 5, decay = 0.1),
                            linout = TRUE, trace = FALSE, MaxNWts = 5000)
# --- 10. Elastic Net (glmnet)
# Учебный смысл: регуляризация (L1/L2), борьба с мультиколлинеарностью, частичный отбор признаков.
glmnet_model<- caret::train(R3haddock ~ ., data = train, method = "glmnet",
                            trControl = ctrl, preProcess = c("center","scale"), tuneLength = 10)
# --- 11. MARS (earth)
# Учебный смысл: кусочно-линейные сплайны + простые взаимодействия; гибкая интерпретация.
earth_model <- caret::train(R3haddock ~ ., data = train, method = "earth",
                            trControl = ctrl, tuneLength = 10)
# --- 12. SVM с радиальным ядром (svmRadial)
# Учебный смысл: ядровой метод; улавливает сложные нелинейности; важна стандартизация.
svm_model   <- caret::train(R3haddock ~ ., data = train, method = "svmRadial",
                            trControl = ctrl, preProcess = c("center","scale"), tuneLength = 8)
# --- 13. k-ближайших соседей (kNN)
# Учебный смысл: простая интуитивная нелинейная модель на расстояниях; чувствительна к масштабу.
knn_model   <- caret::train(R3haddock ~ ., data = train, method = "knn",
                            trControl = ctrl, preProcess = c("center","scale"), tuneLength = 15)
# --- 14. Ranger (быстрый Random Forest)
# Учебный смысл: альтернативная/быстрая реализация леса; сравнить с randomForest.
ranger_model<- caret::train(R3haddock ~ ., data = train, method = "ranger",
                            trControl = ctrl, tuneLength = 3, importance = "impurity")
# --- 15. GBM (классический градиентный бустинг)
# Учебный смысл: другой бустинг деревьев; полезно сравнить с XGBoost.
gbm_model   <- caret::train(R3haddock ~ ., data = train, method = "gbm",
                            trControl = ctrl,
                            tuneGrid = expand.grid(n.trees=100, interaction.depth=1,
                                                   shrinkage=0.1, n.minobsinnode=2),
                            distribution = "gaussian", bag.fraction = 1, verbose = FALSE)
# --- 16. PLS (Partial Least Squares)
# Учебный смысл: проекция на скрытые компоненты с учетом отклика; решает мультиколлинеарность.
pls_model   <- caret::train(R3haddock ~ ., data = train, method = "pls",
                            trControl = ctrl, preProcess = c("center","scale"), tuneLength = 10)
# --- 17. Cubist (правила + деревья)
# Учебный смысл: интерпретируемые правила с комитетами; часто силен на табличных данных.
cubist_model<- caret::train(R3haddock ~ ., data = train, method = "cubist",
                            trControl = ctrl, tuneLength = 5)

# 2.5 Метрики и оценка на тесте
# Замечание: RMSE/MAE — абсолютные ошибки; R2 — доля объяснённой вариации;
# MAPE/sMAPE — относительные ошибки (осторожно при малых значениях отклика).
rmse  <- function(a, p) sqrt(mean((a - p)^2, na.rm = TRUE))
mae   <- function(a, p) mean(abs(a - p), na.rm = TRUE)
r2    <- function(a, p) 1 - sum((a - p)^2, na.rm = TRUE) / sum((a - mean(a))^2, na.rm = TRUE)
mape  <- function(a, p) mean(abs((a - p) / a), na.rm = TRUE) * 100
smape <- function(a, p) mean(2 * abs(p - a) / (abs(a) + abs(p)), na.rm = TRUE) * 100
metrics_vec <- function(y, pred) c(RMSE=rmse(y,pred), MAE=mae(y,pred), R2=r2(y,pred),
                                   MAPE=mape(y,pred), sMAPE=smape(y,pred))

y_test <- test$R3haddock
preds_test <- list(
  LM=predict(lm_model,test), GLM=predict(glm_model,test), GAM=predict(gam_model,test),
  RF=predict(rf_model,test), XGB=predict(xgb_model,test), NNET=predict(nnet_model,test),
  ENet=predict(glmnet_model,test), MARS=predict(earth_model,test), SVM=predict(svm_model,test),
  kNN=predict(knn_model,test), RANGER=predict(ranger_model,test), GBM=predict(gbm_model,test),
  PLS=predict(pls_model,test), CUBIST=predict(cubist_model,test)
)
metrics_table <- do.call(rbind, lapply(names(preds_test), function(nm){
  data.frame(Model = nm, t(metrics_vec(y_test, preds_test[[nm]])), row.names = NULL)
})) %>% arrange(RMSE, MAE)
print(round(metrics_table, 2))

# 2.6 CV-резюме
# Сводим результаты CV по всем моделям и смотрим распределения ошибок.
results <- caret::resamples(list(
  LM=lm_model, GLM=glm_model, GAM=gam_model, RF=rf_model, XGB=xgb_model, NNET=nnet_model,
  ENet=glmnet_model, MARS=earth_model, SVM=svm_model, kNN=knn_model, RANGER=ranger_model,
  GBM=gbm_model, PLS=pls_model, CUBIST=cubist_model
))
summary(results)
lattice::dotplot(results, metric = "RMSE")



# ==============================================================================
# 3) ВЫБОР ЛУЧШЕЙ ПРОГНОСТИЧЕСКОЙ МОДЕЛИ (TIME-SLICE CV НА 3 ГОДА + ХРОНО-ТЕСТ)
# Делим последние годы в тест, внутри train — скользящее окно, h=3.
# ------------------------------------------------------------------------------
# Почему time-slice: временные данные нельзя случайно перемешивать, иначе мы
# «подсматриваем в будущее». Создаём серии обучающих/валидационных окон,
# увеличивая тренировочный период, и тестируем на ближайшем горизонте (3 года).
# ==============================================================================

# 3.1 Данные для time-slice (с YEAR)
model_data <- read.csv("selected_predictors_dataset.csv", header = TRUE, stringsAsFactors = FALSE)
if (!"YEAR" %in% names(model_data)) {
  model_data$YEAR <- seq(1990, by = 1, length.out = nrow(model_data))
}
# Хронологический порядок
model_data <- model_data %>% arrange(YEAR)

# Исходные фичи (исключаем YEAR)
md_for_fit <- model_data %>% select(codTSB, T12, I5, NAOspring, haddock68, R3haddock)

# 3.2 Хронологический holdout (последние годы)
# Идея: отложим ~20% последних лет как полностью внешний тест будущего качества.
n <- nrow(md_for_fit)
holdout_frac <- 0.2
n_test <- max(4, ceiling(n * holdout_frac))
train_ts <- head(md_for_fit, n - n_test)
test_ts  <- tail(md_for_fit, n_test)

# 3.3 Time-slice CV (h=3, expanding window рекомендован: fixedWindow=FALSE)
# initialWindow — размер первого «обучающего» фрагмента; horizon — горизонт
# валидации (здесь 3 года). Далее окно расширяется.
n_train <- nrow(train_ts)
initial_frac <- 0.6
horizon      <- 3
initialWindow <- max(10, floor(initial_frac * n_train))
if (initialWindow + horizon > n_train) initialWindow <- n_train - horizon

slices <- caret::createTimeSlices(1:n_train, initialWindow = initialWindow,
                                  horizon = horizon, fixedWindow = FALSE)
ctrl_ts <- caret::trainControl(method = "cv", index = slices$train, indexOut = slices$test,
                               savePredictions = "final")

# 3.4 Обучение (ядро набора, без GBM — он нестабилен на малом n в timeslice)
# Примечание: используем ту же рецептуру, что и в базовом сравнении, но с
# хронологическими срезами.
fit_ts <- function(method, form, data, ctrl, ...) {
  out <- try(caret::train(form, data = data, method = method, trControl = ctrl, ...), TRUE)
  if (inherits(out,"try-error")) NULL else out
}
lm_ts   <- fit_ts("lm",        R3haddock ~ ., train_ts, ctrl_ts)
glm_ts  <- fit_ts("glm",       R3haddock ~ ., train_ts, ctrl_ts, family = Gamma(link="log"))
gam_ts  <- caret::train(x = train_ts[, -which(names(train_ts)=="R3haddock")],
                        y = train_ts$R3haddock, method = gam_spec, trControl = ctrl_ts)
rf_ts   <- fit_ts("rf",        R3haddock ~ ., train_ts, ctrl_ts, ntree=1000, tuneGrid=data.frame(mtry=1))
xgb_ts  <- fit_ts("xgbTree",   R3haddock ~ ., train_ts, ctrl_ts, tuneGrid = xgb_grid, verbose = 0)
rgr_ts  <- fit_ts("ranger",    R3haddock ~ ., train_ts, ctrl_ts, tuneLength=3)
nnet_ts <- fit_ts("nnet",      R3haddock ~ ., train_ts, ctrl_ts,
                  preProcess=c("center","scale"),
                  tuneGrid=expand.grid(size=5,decay=0.1), linout=TRUE, trace=FALSE, MaxNWts=5000)
svm_ts  <- fit_ts("svmRadial", R3haddock ~ ., train_ts, ctrl_ts, preProcess=c("center","scale"), tuneLength=8)
knn_ts  <- fit_ts("knn",       R3haddock ~ ., train_ts, ctrl_ts, preProcess=c("center","scale"), tuneLength=15)
enet_ts <- fit_ts("glmnet",    R3haddock ~ ., train_ts, ctrl_ts, preProcess=c("center","scale"), tuneLength=10)
mars_ts <- fit_ts("earth",     R3haddock ~ ., train_ts, ctrl_ts, tuneLength=10)
pls_ts  <- fit_ts("pls",       R3haddock ~ ., train_ts, ctrl_ts, preProcess=c("center","scale"), tuneLength=10)
cub_ts  <- fit_ts("cubist",    R3haddock ~ ., train_ts, ctrl_ts, tuneLength=5)

models_ts <- list(LM=lm_ts, GLM=glm_ts, GAM=gam_ts, RF=rf_ts, XGB=xgb_ts, RANGER=rgr_ts,
                  NNET=nnet_ts, SVM=svm_ts, kNN=knn_ts, ENet=enet_ts, MARS=mars_ts, PLS=pls_ts, CUBIST=cub_ts)
models_ts <- models_ts[!vapply(models_ts, is.null, logical(1))]

# 3.5 Ранжирование по time-slice CV и по хронологическому тесту
# Сначала ранжируем по средним ошибкам на валидационных срезах, затем — по внешнему тесту.
cv_metrics <- function(m) {
  if (is.null(m$pred) || !"Resample" %in% names(m$pred)) return(c(RMSE=NA, MAE=NA))
  by_slice <- m$pred %>% group_by(Resample) %>%
    summarise(RMSE=rmse(obs,pred), MAE=mae(obs,pred), .groups="drop")
  c(RMSE = mean(by_slice$RMSE, na.rm = TRUE), MAE = mean(by_slice$MAE, na.rm = TRUE))
}
cv_rank <- do.call(rbind, lapply(models_ts, cv_metrics)) %>% as.data.frame()
cv_rank$Model <- rownames(cv_rank)
cv_rank <- cv_rank[is.finite(cv_rank$RMSE), ] %>% relocate(Model) %>% arrange(RMSE, MAE)
cat("\nTime-slice CV (h=3), средние RMSE/MAE:\n"); print(cv_rank)

preds_ts <- lapply(models_ts, function(m) try(predict(m, newdata = test_ts), TRUE))
keep <- vapply(preds_ts, function(p) is.numeric(p) && length(p)==nrow(test_ts) && all(is.finite(p)), logical(1))
preds_ts <- preds_ts[keep]
test_rank <- do.call(rbind, lapply(names(preds_ts), function(nm){
  data.frame(Model=nm, t(metrics_vec(test_ts$R3haddock, preds_ts[[nm]])), row.names = NULL)
})) %>% arrange(RMSE, MAE)
cat("\nХронологический тест (последние годы), RMSE/MAE/R2:\n"); print(test_rank)



# ==============================================================================
# 4) ПРОГНОЗ 2022–2024 (АНСАМБЛЬ CUBIST+LM) И ГРАФИК 1990–2024 С ДИ
# Прогнозные линии (медиана и ДИ) — пунктир; исторические — сплошные.
# Можно задать свои сценарии предикторов (user_future); по умолчанию — средние.
# ------------------------------------------------------------------------------
# Логика ансамбля: комбинируем сильную нелинейную модель (Cubist) с простой и
# устойчивой линейной (LM). Веса можно настраивать. Доверительные интервалы
# получаем эмпирически из распределения остатков (простая и наглядная эвристика).
# ==============================================================================

# 4.1 Полные модели для прогноза (на всех данных) и вес ансамбля
model_data <- read.csv("selected_predictors_dataset.csv", header = TRUE, stringsAsFactors = FALSE)
if (!"YEAR" %in% names(model_data)) {
  model_data$YEAR <- seq(1990, by = 1, length.out = nrow(model_data))
}
model_data <- model_data %>% arrange(YEAR)

cubist_full <- caret::train(R3haddock ~ codTSB + T12 + I5 + NAOspring + haddock68,
                            data = model_data, method = "cubist",
                            trControl = caret::trainControl(method="none"),
                            tuneGrid = if (exists("cubist_model")) cubist_model$bestTune else NULL,
                            tuneLength = if (exists("cubist_model")) 1 else 5)

lm_full <- caret::train(R3haddock ~ codTSB + T12 + I5 + NAOspring + haddock68,
                        data = model_data, method = "lm",
                        trControl = caret::trainControl(method="none"))

alpha_opt <- if (exists("alpha_opt")) alpha_opt else 0.75
predict_ensemble <- function(newdata, alpha = alpha_opt) {
  alpha * predict(cubist_full, newdata) + (1 - alpha) * predict(lm_full, newdata)
}

# 4.2 Остатки для ДИ (из CV, если есть; иначе — по фитам)
# Эмпирические квантилы остатков дают «практические» интервалы прогноза без
# предположения нормальности ошибок (хотя строгий PI требует аккуратности).
get_residuals_for_pi <- function() {
  if (exists("lm_model") && exists("cubist_model") &&
      !is.null(lm_model$pred) && !is.null(cubist_model$pred)) {
    pl <- lm_model$pred %>% select(Resample,rowIndex,obs,p_lm=pred)
    pc <- cubist_model$pred %>% select(Resample,rowIndex,p_cu=pred)
    inner_join(pl, pc, by=c("Resample","rowIndex")) %>%
      mutate(p_ens = alpha_opt * p_cu + (1 - alpha_opt) * p_lm,
             resid = obs - p_ens) %>%
      pull(resid) %>% .[is.finite(.)]
  } else {
    model_data$R3haddock - predict_ensemble(model_data)
  }
}
resids <- get_residuals_for_pi()
q025 <- as.numeric(quantile(resids, 0.025, na.rm = TRUE))
q250 <- as.numeric(quantile(resids, 0.250, na.rm = TRUE))
q750 <- as.numeric(quantile(resids, 0.750, na.rm = TRUE))
q975 <- as.numeric(quantile(resids, 0.975, na.rm = TRUE))

# 4.3 Сценарии будущего (по умолчанию — средние; можно переопределить user_future)
fc_start <- 2022
pred_cols <- c("codTSB","T12","I5","NAOspring","haddock68")
train_period <- model_data %>% filter(YEAR > 1989 & YEAR < fc_start)
mu <- train_period %>% summarise(across(all_of(pred_cols), ~mean(.x, na.rm = TRUE))) %>% as.list()

# Пример пользовательского сценария:
# user_future <- tibble::tribble(
#   ~YEAR, ~codTSB, ~T12, ~I5, ~NAOspring, ~haddock68,
#   2022, 2100000, 5.1, 48,  0.3, 120000,
#   2023, 2050000, 4.8, 50, -0.1, 115000,
#   2024, 2150000, 5.0, 47,  0.2, 118000
# )
if (!exists("user_future")) user_future <- NULL

build_future <- function(years, mu, user_df=NULL) {
  df <- tibble::tibble(YEAR = years)
  for (v in pred_cols) df[[v]] <- mu[[v]]
  if (!is.null(user_df)) {
    for (i in seq_len(nrow(user_df))) {
      yr <- user_df$YEAR[i]
      if (yr %in% years) {
        idx <- which(df$YEAR == yr)
        for (v in intersect(pred_cols, names(user_df))) {
          val <- user_df[[v]][i]
          if (!is.na(val)) df[[v]][idx] <- val
        }
      }
    }
  }
  df
}
future_years <- fc_start:2024
scenario_future <- build_future(future_years, mu, user_future)

# 4.4 Прогноз и таблица ДИ
pred_future <- predict_ensemble(scenario_future)
forecast_tbl <- tibble::tibble(
  YEAR      = scenario_future$YEAR,
  pred_mean = as.numeric(pred_future),
  PI50_low  = pred_future + q250, PI50_high = pred_future + q750,
  PI95_low  = pred_future + q025, PI95_high = pred_future + q975
)
print(forecast_tbl)

# 4.5 Непрерывный ряд 1990–2024 и график: ленты сплошные; линии медианы/ДИ — сплошные до 2021, пунктир с 2022
pred_df <- bind_rows(
  model_data %>% select(YEAR, all_of(pred_cols)),
  scenario_future
) %>% distinct(YEAR, .keep_all = TRUE) %>% arrange(YEAR)

pred_df$Pred      <- as.numeric(predict_ensemble(pred_df))
pred_df$PI50_low  <- pred_df$Pred + q250
pred_df$PI50_high <- pred_df$Pred + q750
pred_df$PI95_low  <- pred_df$Pred + q025
pred_df$PI95_high <- pred_df$Pred + q975

hist_df <- model_data %>% select(YEAR, R3haddock)

ggplot() +
  geom_ribbon(data = pred_df, aes(x = YEAR, ymin = PI95_low, ymax = PI95_high),
              fill = "grey80", alpha = 0.25) +
  geom_ribbon(data = pred_df, aes(x = YEAR, ymin = PI50_low, ymax = PI50_high),
              fill = "grey60", alpha = 0.35) +
  geom_line(data = subset(pred_df, YEAR < fc_start), aes(x = YEAR, y = PI95_low),
            color = "grey45", linewidth = 0.6) +
  geom_line(data = subset(pred_df, YEAR < fc_start), aes(x = YEAR, y = PI95_high),
            color = "grey45", linewidth = 0.6) +
  geom_line(data = subset(pred_df, YEAR < fc_start), aes(x = YEAR, y = PI50_low),
            color = "grey35", linewidth = 0.6) +
  geom_line(data = subset(pred_df, YEAR < fc_start), aes(x = YEAR, y = PI50_high),
            color = "grey35", linewidth = 0.6) +
  geom_line(data = subset(pred_df, YEAR >= fc_start-1), aes(x = YEAR, y = PI95_low),
            color = "grey45", linewidth = 0.6, linetype = "dashed") +
  geom_line(data = subset(pred_df, YEAR >= fc_start-1), aes(x = YEAR, y = PI95_high),
            color = "grey45", linewidth = 0.6, linetype = "dashed") +
  geom_line(data = subset(pred_df, YEAR >= fc_start-1), aes(x = YEAR, y = PI50_low),
            color = "grey35", linewidth = 0.6, linetype = "dashed") +
  geom_line(data = subset(pred_df, YEAR >= fc_start-1), aes(x = YEAR, y = PI50_high),
            color = "grey35", linewidth = 0.6, linetype = "dashed") +
  geom_line(data = subset(pred_df, YEAR < fc_start), aes(x = YEAR, y = Pred),
            color = "steelblue4", linewidth = 1) +
  geom_line(data = subset(pred_df, YEAR >= fc_start-1), aes(x = YEAR, y = Pred),
            color = "steelblue4", linewidth = 1, linetype = "dashed") +
  geom_point(data = hist_df, aes(x = YEAR, y = R3haddock),
             color = "black", size = 2, alpha = 0.9) +
  scale_x_continuous(expand = expansion(mult = c(0, 0))) +
  labs(
    title = "Пополнение R3haddock: факт (1990–2021) и прогноз (2022–2024)\nАнсамбль CUBIST+LM; непрерывные ДИ, прогноз — пунктир",
    x = "Год", y = "R3haddock"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none")

# 9) Таблица прогноза 2022–2024
print(forecast_tbl)


