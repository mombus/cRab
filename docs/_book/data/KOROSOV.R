# �������� ����������� ���������
library(neuralnet)
library(ggplot2)
###########DOWNLOAD DATA########################################
setwd("C:/COURSES/������")


###1. �������� ���������: ������ �������������� �������������###

# ������: ����� (w) � ����� ���� (lt) �����
w <- c(85, 90, 85, 95, 95, 135, 165, 135, 140)
lt <- c(51, 51, 52, 54, 54, 59, 59, 60, 62)

# ���������� �������� ������
lreg <- lm(w ~ lt)  # w = a0 + a1*lt
summary(lreg)

# ������������
plot(lt, w, main = "����������� ����� �� �����", 
     xlab = "����� ���� (��)", ylab = "����� (�)", pch = 19)
abline(lreg, col = "red", lwd = 2)

###2. ��������� ����������� ���������� �������###

nls_model <- nls(w ~ a0 + a1 * lt, start = list(a0 = 1, a1 = 1))
summary(nls_model)

# ��������� � ���
cat("��� (lm):", coef(lreg), "\n")
cat("�������� (nls):", coef(nls_model))

###3. ������������� ���������: ���� ����������� ��������###
# ����������� ����� ������
w <- c(40, 156, 105, 85, 80, 50, 75, 48, 75, 67)
lt <- c(44, 59, 49, 50, 54, 43, 49, 42, 47, 47)
lc <- c(70, 78, 66, 90, 83, 70, 62, 75, 40, 80)

# ������������� ���������
multi_reg <- glm(w ~ lt + lc)
summary(multi_reg)

###4. ������������� ���������� ������������� ������������###

# ������������ ����� ����������������
log_model <- lm(log(w) ~ log(lt))

# �������������� �������������
a0 <- exp(coef(log_model)[1])
a1 <- coef(log_model)[2]

# ������������
plot(lt, w, main = "��������� �����������", 
     xlab = "����� ���� (��)", ylab = "����� (�)")
curve(a0 * x^a1, add = TRUE, col = "blue", lwd = 2)

###5. ������������� ���������: ��������� ������� � ��������###
# ������ �� ���������� ������
K <- c(100, 126, 158, 200, 251, 316, 398, 501, 631, 794, 1000)
p <- c(0, 0, 0, 0, 0, 0.5, 0.5, 1, 1, 1, 1)
d <- data.frame(K, p)

# ������������� ������
logit_model <- glm(p ~ K, family = binomial(), data = d)

# ������������ S-�������� ������
plot(d$K, d$p, xlab = "������������ ����������", 
     ylab = "���� ��������", main = "������������� ���������")
lines(d$K, predict(logit_model, type = "response"), 
      col = "darkgreen", lwd = 2)

###6. ��������� ����: �� ��������� � ��������� ������������###
# ���������� ��������� ���� (������ �������� ���������)
nn_simple <- neuralnet(p ~ K, data = d, hidden = 0)

# ������������
plot(nn_simple, rep = "best")

###7. ������� ��� ��������������� �������������###

# ���� � ����� ������� ��������
nn_1hidden <- neuralnet(p ~ K, data = d, hidden = 1)

# ��������� � ������������� ����������
plot(d$K, predict(logit_model, type = "response"), 
     type = "l", col = "darkgreen", lwd = 2,
     xlab = "������������", ylab = "����������")
lines(d$K, predict(nn_1hidden, d), col = "blue", lty = 2)
legend("topleft", legend = c("�����-���������", "��������� ����"),
       col = c("darkgreen", "blue"), lty = 1:2)

###8. ������������� � ��������: ����������� ���� ��������###
######8.1 ������� ��� ������ ����������� ���� ����� � ������� ��������###
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

######8.2 ������� ��� ������ ����������� ���� ����� � ����� ��������###
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

######8.3 ������� ��� ������ ����������� ���� ����� � ����� ��������###
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

###9. ��������������� � ���������������� ��������###

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


# ������� ��� ����� �������
new_conditions <- data.frame(
  fo = c(57.9, 35.3, 83.0),
  me = c(4.1, 0.0, 7.3),
  bo = c(3.4, 7.9, 11.5)
)
predict(nv, new_conditions)






