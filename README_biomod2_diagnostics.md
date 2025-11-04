## Диагностика надежности моделей BIOMOD2

Этот репозиторий содержит скрипт `biomod2_diagnostics.R`, дополняющий ваш основной прогон biomod2 статистической диагностикой как для текущего состояния, так и для прогноза (будущих сценариев).

### Что делает скрипт
- Вычисляет и сохраняет:
  - Внутренние метрики моделей и ансамблей (копии результатов `get_evaluations`, `get_variables_importance`).
  - Калибровку (reliability): таблица по бинам, Brier score, ECE, PNG-график.
  - ROC и PR кривые (AUC и AUPRC) с PNG-графиками.
  - Оптимальные пороги (по TSS) для бинаризации прогнозов.
  - Непрерывный индекс Бойса (CBI) и кривую F-ratio.
  - Неопределенность по алгоритмам (SD, CV) для текущего и будущего.
  - Δ-карту (Future − Current) для EMmean (в шкале [0,1]).
  - MESS (экстраполяционный риск) для будущих условий относительно текущих.

Все результаты сохраняются в директорию `diagnostics/run_YYYYMMDD_HHMMSS/` в виде CSV/PNG/RDS.

### Требуемые пакеты
Установите при необходимости:
```r
install.packages(c(
  "biomod2", "dplyr", "tidyr", "purrr", "ggplot2", "readr",
  "pROC", "precrec", "ecospat", "yardstick", "dismo", "ggpubr"
))
```

### Ожидаемые объекты в памяти (после вашего основного скрипта)
- ТЕКУЩЕЕ: `myBiomodModelOut`, `myBiomodProj`, `myBiomodEM`, (опц.) `myBiomodEMProj`, `myExpl`, `myResp`, `DATA`.
- БУДУЩЕЕ (опц.): `myBiomodProj1`, `myBiomodEMProj1`, `myExplP1`, `DATA`.

Скрипт сам попытается собрать недостающий `myBiomodEMProj` из `myBiomodProj` и `myBiomodEM`.

### Как запустить
1) Выполните ваш основной скрипт biomod2 (моделирование, проекция и, при наличии, ансамбль). Убедитесь, что рабочая директория установлена корректно (`getwd()`).
2) Запустите диагностику:
```r
source("biomod2_diagnostics.R")
```
3) Путь к папке с результатами будет выведен в консоль, например:
```
Diagnostics outputs will be written to: /.../diagnostics/run_20250826_153210
```

### Независимая валидация (опционально)
Если в рабочей директории есть файл `independent_validation.csv` с колонками `x, y, occ` и теми же предикторами, что в `myExpl`, скрипт автоматически:
- Спроецирует модели на эти данные,
- Посчитает калибровку, ROC/PR, оптимальные пороги, CBI,
- Сохранит результаты как `*_independent.*`.

### Список ключевых выходов
- Внутренние метрики:
  - `eval_single_models.csv`, `varimp_single_models.csv`, `eval_ensemble.csv`, `varimp_ensemble.csv` (+ соответствующие RDS).
- Калибровка:
  - `calibration_current_table.csv`, `calibration_current.png` (и `*_independent.*`, если есть внешний датасет).
- ROC/PR:
  - `roc_current.png`, `pr_current.png`.
- Порог(и):
  - `optimal_thresholds_current.csv` (и `*_independent.csv`).
- Бойс:
  - `boyce_current.csv`, `boyce_current_curve.png`.
- Неопределенность:
  - `mapdata_current_with_uncertainty.csv`, `uncertainty_sd_current.png`.
  - (Будущее) `mapdata_future_with_uncertainty.csv`, `uncertainty_sd_future.png`.
- Сравнение сценариев:
  - `delta_future_minus_current.csv`, `delta_future_current.png`.
- Экстраполяция:
  - `mapdata_future_with_MESS.csv`, `mess_future_vs_current.png`.

### Примечания
- Масштаб прогнозов: если максимум EMmean > 1.5, прогноз масштабируется в [0,1] делением на 1000 (типично для biomod2). Это влияет на калибровку, ROC/PR, пороги и Δ-карты.
- `MESS` может не посчитаться, если предикторы не совпадают по именам/типам — тогда соответствующие файлы не создаются.
- Количество бинов для калибровки (`num_bins`), шаг порога (`step`), и пр. можно изменить в функции в `biomod2_diagnostics.R`.