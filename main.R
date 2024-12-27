#############################################################
####### code for doubly robust augmented model transfer inference with completely missing covariates
#############################################################
library(BB)
source("functions.R")
str1 = "setting 1: two nuisance models are correctly specified."
str2 = "setting 2: only density ratio model is correct, imputation model is mispecified."
str3 = "setting 3: only imputation model is correct, density ratio model is mispecified."
strs = c(str1,str2,str3)
## choose setting 1 or 2 or 3
setting = 1 
# 1. generate data --------
M = 2000 # total sample size of source and target data
## parameter values
gamma = c(-1, 1, -2) 
eta = c(1, -0.6, -0.5) 
## generate dataset
set.seed(123)
Y = rnorm(M, mean = 0, sd = 1)
Z = rnorm(M, mean = 0, sd = 2) # q = 2
eps = rnorm(M, 0, 0.2)
if(setting == 1){
  X = gamma[1] + gamma[2]*Y + gamma[3]*Z + eps 
  Ps = expit(eta[1] + eta[2]*Y + eta[3]*Z)
}
if(setting == 2){
  X = gamma[1] + gamma[2]*Y + 2*gamma[3]*Z + 0.5*Y*Z + eps 
  Ps = expit(eta[1] + eta[2]*Y + eta[3]*Z)
}
if(setting == 3){
  X = gamma[1] + gamma[2]*Y + gamma[3]*Z + eps
  Ps = expit(2.2*eta[1] + eta[2]*Y + eta[3]*Z - Y*Z)
}
S = rbinom(M,size=1,prob=Ps) 
ns = sum(S) # sample size for source data
nt = M - ns # sample size for target data
data = cbind(X,Y,Z,S) 
# observed source data
Xs = data[which(data[,4]==1),1]
Ys = data[which(data[,4]==1),2] 
Zs = data[which(data[,4]==1),3] 
# observed target data
Yt = data[which(data[,4]==0),2]
Zt = data[which(data[,4]==0),3]

# 2. estimate nuisance models --------
## estimated correct density ratio model: omega(Ys,Zs)
fit = glm(S ~ Y + Z, family = "binomial")
eta_hat = fit$coefficients
omega_hat = exp(-eta_hat[1] - eta_hat[2]*Ys - eta_hat[3]*Zs)

## estimated correct imputation model: m_i(Ys,Zs), m_i(Yt,Zt), i = 1,2
fit = lm(Xs ~ Ys + Zs)
s_eps_hat = sd(fit$residuals)
gamma_hat = fit$coefficients
m1.s = gamma_hat[1] + gamma_hat[2]*Ys + gamma_hat[3]*Zs # m_1(Ys,Zs)
m1.t = gamma_hat[1] + gamma_hat[2]*Yt + gamma_hat[3]*Zt # m_1(Yt,Zt)
m2.s = s_eps_hat^2 + m1.s^2
m2.t = s_eps_hat^2 + m1.t^2

# 3. estimate parameters using the proposed methods --------
## input
data = data # data frame of (X,Y,Z,S): source data (S=1) and target data (S=0)
omega_hat = omega_hat # estimated density ratio
m1.s = m1.s; m2.s = m2.s # imputed X and X^2 for source data
m1.t = m1.t; m2.t = m2.t # imputed X and X^2 for target data
## run our method 
result = fn_dr(data,omega_hat,m1.s,m1.t,m2.s,m2.t)
## output
print(strs[setting])
round(result,3) # point estimators and standard errors for three methods: IW, IMP, and DR (proposed doubly robust method)




