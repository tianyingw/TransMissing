## functions
library(BB) # dfsane: Derivative-Free Spectral Approach for solving nonlinear systems of equations
expit<-function(x){
  1/(1+exp(-x))
}
# Comparable method 1: Importance Weighting
ImportanceWeight_Only<-function(omega, beta_theta_init=NULL){
  # weighted estimating equations
  augeq<-function(beta_theta){
    temp = t(cbind(Xs,1,Zs)) %*% (omega*(Ys - beta_theta[1] - Xs*beta_theta[2] - Zs*beta_theta[3]))/ns
    temp[,1]
  }
  if(is.null(beta_theta_init)){
    beta_theta_init = rep(0,3) # p+q
  }
  beta_theta_iw = dfsane(beta_theta_init,augeq,control=list(maxit=2500,trace = FALSE))$par
  names(beta_theta_iw) = c("Intercept", "Xt", "Zt")
  return(beta_theta_iw)
}
# Comparable method 2: Imputation 
Imputation_Only<-function(m1.t, m2.t, beta_theta_init=NULL){
  augeq<-function(beta_theta){
    temp1 = t(m1.t) %*% (Yt - beta_theta[1] - Zt*beta_theta[3])/nt - mean(m2.t)*beta_theta[2]
    temp2 = t(cbind(1,Zt)) %*% (Yt - beta_theta[1] - m1.t*beta_theta[2] - Zt*beta_theta[3])/nt
    temp = c(temp1, temp2)
    return(temp)
  }
  if(is.null(beta_theta_init)){
    beta_theta_init = rep(0,3) # p+q
  }
  beta_theta_imp = dfsane(beta_theta_init,augeq,control=list(maxit=2500,trace = FALSE))$par
  names(beta_theta_imp) = c("Intercept", "Xt", "Zt")
  return(beta_theta_imp)
}
# Proposed Doubly Robust method
Doubly_Robust <- function(omega, m1.s, m2.s, m1.t, m2.t, beta_theta_init=NULL){
  augeq<-function(beta_theta){
    temp1 = t(Xs - m1.s) %*% (omega*(Ys - beta_theta[1] - Zs*beta_theta[3])) / ns + t(m2.s - Xs^2) %*% omega * beta_theta[2] / ns 
    temp1 = temp1 + t(m1.t) %*% (Yt - beta_theta[1] - Zt*beta_theta[3]) / nt - mean(m2.t)*beta_theta[2]
    temp2 = t(cbind(1,Zs)) %*% (omega*(m1.s - Xs)) * beta_theta[2] / ns
    temp2 = temp2 + t(cbind(1,Zt)) %*% (Yt - beta_theta[1] - m1.t*beta_theta[2]- Zt*beta_theta[3]) / nt
    temp = c(temp1, temp2)
    return(temp)
  }
  if(is.null(beta_theta_init)){
    beta_theta_init = rep(0,3) # p+q
  }
  beta_theta_dr = dfsane(beta_theta_init,augeq,control=list(maxit=2500,trace = FALSE))$par
  names(beta_theta_dr) = c("Intercept", "Xt", "Zt")
  return(beta_theta_dr)
}

# functions for estimating parameters using above three proposed methods
fn_dr = function(data,omega_hat,m1.s,m1.t,m2.s,m2.t){
  B = 500  # 500 Bootstrap for variance estimation
  beta_theta_naive = beta_theta_iw = beta_theta_im = beta_theta_dr = rep(NA,3)
  b.se.naive = b.se.iw = b.se.im = b.se.dr = rep(NA,3) # intercept, Xt, Zt
  b.beta.theta.naive = b.beta.theta.iw= b.beta.theta.im = b.beta.theta.dr = matrix(NA,B,3) # estimator from bootstrap
  
  X = data[,1]; Y = data[,2]; Z = data[,3]; S = data[,4]
  ## source data
  Xs = data[which(data[,4]==1),1]
  Ys = data[which(data[,4]==1),2]
  Zs = data[which(data[,4]==1),3]
  ## target data
  Xt = data[which(data[,4]==0),1] # artificially missing 
  Yt = data[which(data[,4]==0),2]
  Zt = data[which(data[,4]==0),3]
  
  X0 = X; Y0 = Y; Z0 = Z; S0 = S
  Xs0 = Xs; Ys0 = Ys; Zs0 = Zs
  Xt0 = Xt; Yt0 = Yt; Zt0 = Zt
  
  ## naive method: using source data
  fit <- lm(Ys ~ Xs + Zs)
  beta_theta_naive = fit$coefficients
  
  ## IW method
  beta_theta_iw <- tryCatch({
    ImportanceWeight_Only(omega_hat, beta_theta_init=beta_theta_naive)
  }, warning = function(w) {
    message("Warning：", conditionMessage(w))
    return(NA)
  })
  ## IMP method
  beta_theta_im <- tryCatch({
    Imputation_Only(m1.t, m2.t, beta_theta_init=beta_theta_naive)
  }, warning = function(w) {
    message("Warning：", conditionMessage(w))
    return(NA)
  })
  ## Proposed method
  beta_theta_dr <- tryCatch({
    Doubly_Robust(omega_hat, m1.s, m2.s, m1.t, m2.t, beta_theta_init=beta_theta_naive) 
  }, warning = function(w) {
    message("Warning：", conditionMessage(w))
    return(NA)
  })
  
  ## bootstrap for variance estimation------------
  for(b in 1:B){
    b.s.data = matrix(0,ns,3) # source: X,Y,Z,S
    b.s.data[,1] = Xs0
    b.s.data[,2] = Ys0
    b.s.data[,3] = Zs0
    index = sample(1:ns,size = ns,replace = TRUE)
    b.s.data = b.s.data[index,]
    Xs = b.s.data[,1]
    Ys = b.s.data[,2]
    Zs = b.s.data[,3]
    
    b.t.data = matrix(0,nt,3) # target: X (artificially missing),Y,Z,S
    b.t.data[,1] = Xt0
    b.t.data[,2] = Yt0
    b.t.data[,3] = Zt0
    index = sample(1:nt,size = nt,replace = TRUE)
    b.t.data = b.t.data[index,]
    Xt = b.t.data[,1] # artificially missing
    Yt = b.t.data[,2]
    Zt = b.t.data[,3]
    
    X = c(Xs,Xt); Y = c(Ys,Yt); Z = c(Zs,Zt); S = c(rep(1,ns),rep(0,nt))
    
    ## nuisance models--------------
    # estimated correct density ratio model: omega(Ys,Zs)
    fit = glm(S ~ Y + Z, family = "binomial")
    eta_hat = fit$coefficients
    omega_hat = exp(-eta_hat[1] - eta_hat[2]*Ys - eta_hat[3]*Zs)
    
    # estimated correct imputation model: m_i(Ys,Zs), m_i(Yt,Zt), i = 1,2
    fit = lm(Xs ~ Ys + Zs)
    s_eps_hat = sd(fit$residuals)
    gamma_hat = fit$coefficients
    m1.s = gamma_hat[1] + gamma_hat[2]*Ys + gamma_hat[3]*Zs # m_1(Ys,Zs)
    m1.t = gamma_hat[1] + gamma_hat[2]*Yt + gamma_hat[3]*Zt # m_1(Yt,Zt)
    m2.s = s_eps_hat^2 + m1.s^2
    m2.t = s_eps_hat^2 + m1.t^2
    
    ## different methods--------------
    fit <- lm(Ys ~ Xs + Zs)
    b.beta.theta.naive[b,] = fit$coefficients
    
    b.beta.theta.iw[b,] <- tryCatch({
      ImportanceWeight_Only(omega_hat, beta_theta_init=b.beta.theta.naive[b,])
    }, warning = function(w) {
      #message("Warning：", conditionMessage(w))
      return(NA)
    })
    b.beta.theta.im[b,] <- tryCatch({
      Imputation_Only(m1.t, m2.t, beta_theta_init=b.beta.theta.naive[b,]) # might be bad if m2.t has a large bias
    }, warning = function(w) {
      #message("Warning：", conditionMessage(w))
      return(NA)
    })
    b.beta.theta.dr[b,] <- tryCatch({
      Doubly_Robust(omega_hat, m1.s, m2.s, m1.t, m2.t, beta_theta_init=b.beta.theta.naive[b,]) 
    }, warning = function(w) {
      #message("Warning：", conditionMessage(w))
      return(NA)
    })
  }
  for(j in 1:3){
    b.se.naive[j] = sd(na.omit(b.beta.theta.naive[,j]))
    b.se.iw[j]    = sd(na.omit(b.beta.theta.iw[,j]))
    b.se.im[j]    = sd(na.omit(b.beta.theta.im[,j]))
    b.se.dr[j]    = sd(na.omit(b.beta.theta.dr[,j]))
  }
  
  names(beta_theta_iw) = c("Intercept", "Xt", "Zt")
  result= data.frame(estimate.iw = beta_theta_iw, estimate.imp = beta_theta_im, estimate.dr = beta_theta_dr, 
                     se.iw = b.se.iw, se.imp = b.se.im, se.dr = b.se.dr)
  return(result)
}
