# Загрузка необходимых библиотек
library(neuralnet)
library(ggplot2)
###########DOWNLOAD DATA########################################
setwd("C:/COURSES/ЛЕКЦИИ")


###1. Линейная регрессия: основы экологического моделирования###

# Данные: масса (w) и длина тела (lt) гадюк
w <- c(85, 90, 85, 95, 95, 135, 165, 135, 140)
lt <- c(51, 51, 52, 54, 54, 59, 59, 60, 62)

# Построение линейной модели
lreg <- lm(w ~ lt)  # w = a0 + a1*lt
summary(lreg)

# Визуализация
plot(lt, w, main = "Зависимость массы от длины", 
     xlab = "Длина тела (см)", ylab = "Масса (г)", pch = 19)
abline(lreg, col = "red", lwd = 2)

###2. Численная оптимизация параметров моделей###

nls_model <- nls(w ~ a0 + a1 * lt, start = list(a0 = 1, a1 = 1))
summary(nls_model)

# Сравнение с МНК
cat("МНК (lm):", coef(lreg), "\n")
cat("Подгонка (nls):", coef(nls_model))

###3. Множественная регрессия: учет комплексных факторов###
# Расширенный набор данных
w <- c(40, 156, 105, 85, 80, 50, 75, 48, 75, 67)
lt <- c(44, 59, 49, 50, 54, 43, 49, 42, 47, 47)
lc <- c(70, 78, 66, 90, 83, 70, 62, 75, 40, 80)

# Множественная регрессия
multi_reg <- glm(w ~ lt + lc)
summary(multi_reg)

###4. Моделирование нелинейных экологических зависимостей###

# Линеаризация через логарифмирование
log_model <- lm(log(w) ~ log(lt))

# Преобразование коэффициентов
a0 <- exp(coef(log_model)[1])
a1 <- coef(log_model)[2]

# Визуализация
plot(lt, w, main = "Степенная зависимость", 
     xlab = "Длина тела (см)", ylab = "Масса (г)")
curve(a0 * x^a1, add = TRUE, col = "blue", lwd = 2)

###5. Логистическая регрессия: пороговые эффекты в экологии###
# Данные по смертности дафний
K <- c(100, 126, 158, 200, 251, 316, 398, 501, 631, 794, 1000)
p <- c(0, 0, 0, 0, 0, 0.5, 0.5, 1, 1, 1, 1)
d <- data.frame(K, p)

# Логистическая модель
logit_model <- glm(p ~ K, family = binomial(), data = d)

# Визуализация S-образной кривой
plot(d$K, d$p, xlab = "Концентрация токсиканта", 
     ylab = "Доля погибших", main = "Логистическая регрессия")
lines(d$K, predict(logit_model, type = "response"), 
      col = "darkgreen", lwd = 2)

###6. Нейронные сети: от регрессии к нейронным конструкциям###
# Простейшая нейронная сеть (аналог линейной регрессии)
nn_simple <- neuralnet(p ~ K, data = d, hidden = 0)

# Визуализация
plot(nn_simple, rep = "best")

###7. Нейроны как преобразователи нелинейностей###

# Сеть с одним скрытым нейроном
nn_1hidden <- neuralnet(p ~ K, data = d, hidden = 1)

# Сравнение с логистической регрессией
plot(d$K, predict(logit_model, type = "response"), 
     type = "l", col = "darkgreen", lwd = 2,
     xlab = "Концентрация", ylab = "Смертность")
lines(d$K, predict(nn_1hidden, d), col = "blue", lty = 2)
legend("topleft", legend = c("Логит-регрессия", "Нейронная сеть"),
       col = c("darkgreen", "blue"), lty = 1:2)

###8. Классификация в экологии: определение пола животных###
######8.1 Расчеты для модели диагностики пола гадюк с нулевым нейроном###
v<-read.csv( "vipkar.csv")

head(v,3)

plot(v$lc,v$ns,ylim=c(-0.1,1.1))
nv<-neuralnet(formula = ns~lc, data = v, hidden = 0)
wn<-(as.vector(nv$net.result[[1]]))
points(v$lc,wn,pch=16)
legend('left',legend=c(1:2),pch=c(1,16))
plot(nv)
plot(v$ns,wn)
abline(h=0.5,lty=2)
sum(v$ns==round(wn,0))/nrow(v)

######8.2 Расчеты для модели диагностики пола гадюк с одним нейроном###
v<-read.csv( "vipkar.csv")

head(v,3)

plot(v$lc,v$ns,ylim=c(-0.1,1.1))
nv<-neuralnet(formula = ns~lc, data = v, hidden = 1)
wn<-(as.vector(nv$net.result[[1]]))
points(v$lc,wn,pch=16)
legend('left',legend=c(1:2),pch=c(1,16))
plot(nv)
plot(v$ns,wn)
abline(h=0.5,lty=2)
sum(v$ns==round(wn,0))/nrow(v)

######8.3 Расчеты для модели диагностики пола гадюк с тремя нейроном###
v<-read.csv( "vipkar.csv")

head(v,3)

plot(v$lc,v$ns,ylim=c(-0.1,1.1))
nv<-neuralnet(formula = ns~lc+lt+w, data = v, hidden = 3)
wn<-(as.vector(nv$net.result[[1]]))
points(v$lc,wn,pch=16)
legend('left',legend=c(1:2),pch=c(1,16))
plot(nv)
plot(v$ns,wn)
abline(h=0.5,lty=2)
sum(v$ns==round(wn,0))/nrow(v)

###9. Прогнозирование в пространственной экологии###

v<-read.csv( "kihzsdat.csv")
head(v,3)

nro<-sample(1:24,12) ; d<-v[nro,]
nv<-neuralnet(formula = vb ~ fo+me+bo, data = d, hidden = c(5))
nr<-compute(x=nv, covariate=d)$net.result
rnr<-round(nr,0)
data.frame(d$name,d$vb,rnr)

sum(d$vb==rnr)/nrow(d)

nro2<-sample(1:24,12) ; d2<-v[nro2,]
nr2<-compute(x=nv, covariate=d2)$net.result
rnr2<-round(nr2,0)
data.frame(d2$name,d2$vb,rnr2)

sum(d2$vb==rnr2)/nrow(d)


# Прогноз для новых условий
new_conditions <- data.frame(
  fo = c(57.9, 35.3, 83.0),
  me = c(4.1, 0.0, 7.3),
  bo = c(3.4, 7.9, 11.5)
)
predict(nv, new_conditions)






