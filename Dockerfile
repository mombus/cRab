FROM r-base:4.3.0

# Установка системных зависимостей
RUN apt-get update && apt-get install -y \
    libjags-dev \
    libgsl-dev \
    && rm -rf /var/lib/apt/lists/*

# Установка R пакетов
RUN R -e "install.packages(c('rjags', 'coda', 'ggplot2', 'dplyr', 'tidyr', 'gridExtra', 'viridis', 'scales', 'bayesplot', 'MCMCvis'), repos='https://cloud.r-project.org')"

# Копирование скрипта
COPY salmon_bayesian_model.R /app/

WORKDIR /app

CMD ["Rscript", "salmon_bayesian_model.R"]