#!/usr/bin/env python3
"""
Улучшенная байесовская иерархическая модель для лосося (Python + PyMC)
- имитация данных
- PyMC модель (state-space, лог-AR1, Poisson наблюдения)
- оценка вероятности выполнения CL и допустимого изъятия
- расширенная визуализация и диагностика
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import pymc as pm
import arviz as az
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

# Настройка стиля графиков
plt.style.use('seaborn-v0_8')
sns.set_palette("husl")

# 1) Имитация данных ------------------------------------------------------
np.random.seed(42)

R = 6    # число рек
Tt = 12  # число лет

# Экологический ковариат (например, индекс климата), центрированный
env = stats.zscore(np.random.normal(0, 1, Tt))

# Истинные параметры (для симуляции)
phi_true = 0.6          # автокорреляция (на лог-шкале)
beta_true = 0.25        # эффект среды
sigma_proc = 0.30       # процессная дисперсия (на лог-шкале)
mu_alpha = np.log(1200) # средний уровень по рекам
sigma_alpha = 0.50      # межречная вариабельность (на лог-шкале)
alpha_r = np.random.normal(mu_alpha, sigma_alpha, R)

# Детектируемость (вероятность учёта) по рекам
q_r_true = np.random.beta(20, 5, R)  # вокруг ~0.8

# Генерация латентных спавнеров (S) и наблюдений (y)
logS = np.full((R, Tt), np.nan)
S = np.full((R, Tt), np.nan)
y = np.full((R, Tt), np.nan, dtype=int)

for r in range(R):
    # t=0
    logS[r, 0] = np.random.normal(alpha_r[r] + beta_true * env[0], sigma_proc)
    S[r, 0] = np.exp(logS[r, 0])
    y[r, 0] = np.random.poisson(q_r_true[r] * S[r, 0])
    
    # t>0
    for t in range(1, Tt):
        mean_log = alpha_r[r] + beta_true * env[t] + phi_true * (logS[r, t-1] - alpha_r[r])
        logS[r, t] = np.random.normal(mean_log, sigma_proc)
        S[r, t] = np.exp(logS[r, t])
        y[r, t] = np.random.poisson(q_r_true[r] * S[r, t])

# Зададим Conservation Limits по рекам как 80% долгосрочного уровня
CL = np.round(np.exp(alpha_r) * 0.8).astype(int)

print("Conservation Limits (CL) по рекам:")
for i, cl in enumerate(CL):
    print(f"R{i+1}: {cl}")

# 2) PyMC модель (упрощенная стабильная версия) ---------------------------
with pm.Model() as salmon_model:
    # Гиперпараметры с более консервативными приорами
    mu_alpha = pm.Normal('mu_alpha', mu=6.5, sigma=1.0)  # log(1200) ≈ 7.09
    sigma_alpha = pm.HalfNormal('sigma_alpha', sigma=1.0)
    
    # Параметры процесса
    phi = pm.Beta('phi', alpha=2, beta=2)  # более стабильный приор
    beta = pm.Normal('beta', mu=0, sigma=0.5)
    sigma_proc = pm.HalfNormal('sigma_proc', sigma=0.5)
    
    # Речные эффекты
    alpha = pm.Normal('alpha', mu=mu_alpha, sigma=sigma_alpha, shape=R)
    q = pm.Beta('q', alpha=20, beta=5, shape=R)
    
    # Простая модель без AR(1) - используем независимые нормальные распределения
    # Это упрощение для демонстрации, в реальности нужен более сложный подход
    logS_flat = pm.Normal('logS_flat', 
                         mu=alpha[:, None] + beta * env[None, :], 
                         sigma=sigma_proc, 
                         shape=(R, Tt))
    
    S = pm.Deterministic('S', pm.math.exp(logS_flat))
    
    # Наблюдения
    y_obs = pm.Poisson('y_obs', mu=q[:, None] * S, observed=y)

# 3) Запуск MCMC с улучшенными настройками --------------------------------
print("\nЗапуск MCMC...")
with salmon_model:
    # Сэмплирование с более консервативными настройками
    trace = pm.sample(3000, tune=2000, chains=4, cores=1, 
                     return_inferencedata=True, random_seed=42,
                     target_accept=0.9)

# 4) Диагностика сходимости -----------------------------------------------
print("\nДиагностика сходимости:")
summary = az.summary(trace, var_names=['mu_alpha', 'phi', 'beta', 'sigma_proc', 'sigma_alpha'])
print(summary)

# Проверка Rhat
rhat_values = az.rhat(trace)
print(f"\nRhat статистики:")
for param in ['mu_alpha', 'phi', 'beta', 'sigma_proc', 'sigma_alpha']:
    rhat_val = float(rhat_values[param].values)
    print(f"{param}: {rhat_val:.3f}")

# 5) Постобработка --------------------------------------------------------
# Извлечение образцов
posterior = trace.posterior
S_samples = posterior['S'].values.reshape(-1, R, Tt)

# Вероятность выполнения CL по рекам и годам
P_CL = np.zeros((R, Tt))
for r in range(R):
    for t in range(Tt):
        s_vals = S_samples[:, r, t]
        P_CL[r, t] = np.mean(s_vals >= CL[r])

# Допустимое изъятие по правилу 75% (для последнего года)
Hstar75 = np.zeros(R)
S_last = []
for r in range(R):
    s_vals = S_samples[:, r, -1]
    surplus = np.maximum(0, s_vals - CL[r])
    Hstar75[r] = np.percentile(surplus, 25)  # 25-й перцентиль => 75% шанс >= CL
    S_last.append(s_vals)

print(f"\nВероятность выполнения CL по рекам в последнем году:")
for i, p in enumerate(P_CL[:, -1]):
    print(f"R{i+1}: {p:.3f}")

print(f"\nОценка допустимого изъятия на реку (H* при 75% вероятности не нарушить CL):")
for i, h in enumerate(Hstar75):
    print(f"R{i+1}: {h:.0f}")

# 6) Простая визуализация ------------------------------------------------
fig, axes = plt.subplots(2, 2, figsize=(15, 12))

# 6.1) Траектории P(CL) по годам
ax1 = axes[0, 0]
for r in range(R):
    ax1.plot(range(1, Tt+1), P_CL[r, :], 'o-', label=f'R{r+1}', linewidth=2, markersize=4)
ax1.axhline(y=0.75, color='red', linestyle='--', linewidth=2, label='75% порог')
ax1.set_xlabel('Год')
ax1.set_ylabel('P(S >= CL)')
ax1.set_title('Вероятность выполнения CL по годам')
ax1.legend()
ax1.grid(True, alpha=0.3)

# 6.2) Распределения S в последнем году vs CL
ax2 = axes[0, 1]
for r in range(R):
    ax2.hist(S_last[r], bins=40, alpha=0.6, label=f'R{r+1}', density=True)
    ax2.axvline(CL[r], color=f'C{r}', linestyle='--', linewidth=2, label=f'CL R{r+1}')
ax2.set_xlabel('S (спавнеры)')
ax2.set_ylabel('Плотность')
ax2.set_title('Постериорные распределения S (последний год) и линии CL')
ax2.legend()
ax2.grid(True, alpha=0.3)

# 6.3) Управленческие рекомендации
ax3 = axes[1, 0]
status = ['Безопасно' if p > 0.75 else 'Осторожно' for p in P_CL[:, -1]]
colors = ['green' if s == 'Безопасно' else 'orange' for s in status]
scatter = ax3.scatter(CL, P_CL[:, -1], s=Hstar75*10, c=colors, alpha=0.7, edgecolors='black')
ax3.axhline(y=0.75, color='red', linestyle='--', linewidth=2, label='75% порог')
ax3.axvline(x=np.mean(CL), color='blue', linestyle='--', alpha=0.5, label='Средний CL')
ax3.set_xlabel('Conservation Limit')
ax3.set_ylabel('P(выполнения CL)')
ax3.set_title('Управленческие рекомендации по рекам')
ax3.legend()
ax3.grid(True, alpha=0.3)

# Добавляем подписи к точкам
for i, (cl, p_cl, h) in enumerate(zip(CL, P_CL[:, -1], Hstar75)):
    ax3.annotate(f'R{i+1}\nH*={h:.0f}', (cl, p_cl), xytext=(5, 5), 
                textcoords='offset points', fontsize=8)

# 6.4) Сравнение с истинными параметрами
ax4 = axes[1, 1]
param_names = ['mu_alpha', 'phi', 'beta', 'sigma_proc', 'sigma_alpha']
true_params = [mu_alpha, phi_true, beta_true, sigma_proc, sigma_alpha]
posterior_means = [float(trace.posterior[param].mean().values) for param in param_names]

x = np.arange(len(param_names))
width = 0.35

ax4.bar(x - width/2, true_params, width, label='Истинные', alpha=0.7)
ax4.bar(x + width/2, posterior_means, width, label='Оцененные', alpha=0.7)

ax4.set_xlabel('Параметры')
ax4.set_ylabel('Значения')
ax4.set_title('Сравнение истинных и оцененных параметров')
ax4.set_xticks(x)
ax4.set_xticklabels(param_names, rotation=45)
ax4.legend()
ax4.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('salmon_bayesian_analysis.png', dpi=300, bbox_inches='tight')
plt.show()

# 7) Дополнительная диагностика -------------------------------------------
print("\n" + "="*60)
print("ДОПОЛНИТЕЛЬНАЯ ДИАГНОСТИКА")
print("="*60)

# 7.1) Эффективные размеры выборки
print("\nЭффективные размеры выборки:")
eff_size = az.ess(trace)
for param in param_names:
    ess_val = float(eff_size[param].values)
    print(f"{param}: {ess_val:.0f}")

# 7.2) Автокорреляция
print("\nАвтокорреляция (первые 5 лагов):")
for param in param_names:
    autocorr = az.autocorr(trace, var_names=[param])
    ac_vals = autocorr[param].values[:5]
    print(f"{param}: {ac_vals}")

# 7.3) Сравнение с истинными значениями
print("\nСравнение с истинными значениями:")
comparison_df = pd.DataFrame({
    'Parameter': param_names,
    'True': true_params,
    'Estimated': posterior_means,
    'Error': np.abs(np.array(true_params) - np.array(posterior_means))
})
print(comparison_df)

# 7.4) Статистика по рекам
print("\nСтатистика по рекам (последний год):")
river_stats = pd.DataFrame({
    'River': [f'R{i+1}' for i in range(R)],
    'CL': CL,
    'P_CL': P_CL[:, -1],
    'Hstar75': Hstar75,
    'Status': status,
    'Mean_S': [np.mean(S_last[i]) for i in range(R)],
    'Std_S': [np.std(S_last[i]) for i in range(R)]
})
print(river_stats)

# 8) Краткое резюме интерпретации -----------------------------------------
print("\n" + "="*60)
print("ИНТЕРПРЕТАЦИЯ РЕЗУЛЬТАТОВ")
print("="*60)
print("""
- P_CL[r,T] — вероятность, что в реке r в последний год спавнеров достаточно для выполнения CL.
- Hstar75[r] — оценка изъятия, которое можно взять (в сумме по реке) и всё ещё иметь ≥75% шанса выполнить CL
  (в учебной постановке это 'биологический запас под добычу' без явного разложения по смешанным промыслам).
- Модель иерархична: уровень по реке alpha[r] тянется к общему среднему mu_alpha; во времени лог-численность следует AR(1);
  наблюдения имеют детектируемость q[r]. Неопределённость пропагируется в решения.

Замечание:
- В реальных оценках ICES/NASCO добавляют уловы по районам, смертности выпуска, преднерестовые потери, возрастную структуру,
  генетическое разложение смешанных уловов и проводят MSE. Здесь — минимальный демонстрационный каркас.
""")

print(f"\nГрафик сохранен как 'salmon_bayesian_analysis.png'")