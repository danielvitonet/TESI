# ledermann bound (maximum number of factors)
ledermann <- function(m) {
  as.integer(floor((2*m+1)/2 - sqrt((2*m+1)^2/4 - m^2 + m)))
}

# functions for calculating KL divergence
my.ecdf  <-  function(x)   {
  x   <-   sort(x)
  x.u <-   unique(x)
  n  <-  length(x) 
  x.rle  <-  rle(x)$lengths
  y  <-  (cumsum(x.rle)-0.5) / n
  FUN  <-  approxfun(x.u, y, method="linear", yleft=0, yright=1,
                     rule=2)
  return(FUN)
}   

KL_est  <-  function(x, y)   {
  dx  <-  diff(sort(unique(x)))
  dy  <-  diff(sort(unique(y)))
  ex  <-  min(dx) ; ey  <-  min(dy)
  e   <-  min(ex, ey)/2
  n   <-  length(x)    
  P  <-   my.ecdf(x) ; Q  <-  my.ecdf(y)
  KL  <-  sum( log( (P(x)-P(x-e))/(Q(x)-Q(x-e)))) / n
  return(KL)              
}

# function for quantile weighted CRPS
qwCRPS_sample <-function(true,mcmc,tau,weighting="none"){
  require(pracma)
  
  tau_len <- length(tau)
  Q.tau <- quantile(mcmc,probs=tau)
  true_rep <- rep(true,tau_len)
  QS.vec <- (true_rep-Q.tau)*(tau-((true_rep<=Q.tau)*1))
  
  weights <- switch(tolower(weighting),
                    "none" = 1,
                    "tails" = (2*tau-1)^2,
                    "right" = tau^2,
                    "left" = (1-tau)^2,
                    "center" = tau*(1-tau))
  wghs <- QS.vec*weights
  return(pracma::trapz(tau,wghs))
}

# This file contains several auxiliary functions.
ncind <- function(y,mu,sig,q){
  sample(1:length(q),
         size = 1,
         prob = q*dnorm(y,mu,sig))
}

# get posteriors for the horseshoe prior (see Makalic & Schmidt, 2015)
get.hs <- function(bdraw,lambda.hs,nu.hs,tau.hs,zeta.hs){
  k <- length(bdraw)
  if (is.na(tau.hs)){
    tau.hs <- 1   
  }else{
    tau.hs <- invgamma::rinvgamma(1,shape=(k+1)/2,rate=1/zeta.hs+sum(bdraw^2/lambda.hs)/2) 
  }
  
  lambda.hs <- invgamma::rinvgamma(k,shape=1,rate=1/nu.hs+bdraw^2/(2*tau.hs))
  
  nu.hs <- invgamma::rinvgamma(k,shape=1,rate=1+1/lambda.hs)
  zeta.hs <- invgamma::rinvgamma(1,shape=1,rate=1+1/tau.hs)
  
  ret <- list("psi"=(lambda.hs*tau.hs),"lambda"=lambda.hs,"tau"=tau.hs,"nu"=nu.hs,"zeta"=zeta.hs)
  return(ret)
}

# lag variables
mlag <- function(X,lag){
  p <- lag
  X <- as.matrix(X)
  Traw <- nrow(X)
  N <- ncol(X)
  Xlag <- matrix(NA,Traw,p*N)
  for (ii in 1:p){
    Xlag[(p+1):Traw,(N*(ii-1)+1):(N*ii)]=X[(p+1-ii):(Traw-ii),(1:N)]
  }
  return(Xlag)
}

# companion matrix
get_companion <- function(Beta_,varndxv){
  nn <- varndxv[[1]]
  nd <- varndxv[[2]]
  nl <- varndxv[[3]]
  
  nkk <- nn*nl+nd
  
  Jm <- matrix(0,nkk,nn)
  Jm[1:nn,1:nn] <- diag(nn)
  
  MM <- rbind(t(Beta_),cbind(diag((nl-1)*nn), matrix(0, (nl-1)*nn, nn)))
  
  return(list(MM=MM,Jm=Jm))
}
remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
  H <- 15 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  y
}

# prior moments for Minnesota prior
get_V <- function(a_bar_1,a_bar_2,ind1,sigma_sq1,p1,set.cons,M1,K1){
  V_i <- matrix(0,K1,M1)
  #this double loop fills the prior covariance matrix
  for (i in 1:M1){ #for each i equation
    for (j in 1:K1){ #for each variable on the rhs
      if (set.cons==1) {
        if (j==K1) {
          V_i[j,i] <- 10^3*sigma_sq1[i,1] #variance on constant, trend, dummies and on ex variables         
        }    
        else if (any(j==ind1[i,])) {
          ll <- which(ind1[i,]==j)
          V_i[j,i] <- a_bar_1/(ll^2) #variance on own lags
        }else{
          ll <- which(ind1==j,arr.ind=TRUE)[2]
          kj <- ind1[which(ind1==j,arr.ind=TRUE)[1],1]
          V_i[j,i] <- (a_bar_2*sigma_sq1[i,1])/((ll^2)*sigma_sq1[kj,1])
        }
      }else{
        if (any(j==ind1[i,])) {
          ll <- which(ind1[i,]==j)
          V_i[j,i] <- a_bar_1/(ll^2) #variance on own lags
        } 
        else{
          ll <- which(ind1==j,arr.ind=TRUE)[2]
          kj <- ind1[which(ind1==j,arr.ind=TRUE)[1],1]
          V_i[j,i] <- (a_bar_2*sigma_sq1[i,1])/((ll^2)*sigma_sq1[kj,1])
        }
      }
    }
  }  
  return(V_i)
}

# -----------------------------------------------------------------------------------------------
# function to draw the factor loadings (basic linear regression)
get.facload <- function(yy,xx,l_sd){
  V_prinv <- diag(NCOL(xx))/l_sd
  V_lambda <- solve(crossprod(xx) + V_prinv)
  lambda_mean <- V_lambda %*% (crossprod(xx,yy))
  
  lambda_draw <- lambda_mean + t(chol(V_lambda)) %*% rnorm(NCOL(xx))
  return(lambda_draw)
}

# factor loadings draw
get.Lambda <- function(eps,fac,S,pr,m,q,id.fac){
  L <- matrix(0,m,q)
  if(id.fac){
    for(jj in 1:m){
      if (jj<=q){
        normalizer <- exp(0.5*S[,jj])
        yy0 <- (eps[,jj]-fac[,jj])/normalizer
        xx0 <- fac[,1:(jj-1),drop=FALSE]/normalizer
        if (jj>1){
          l_sd <- pr[jj,1:(jj-1)]
          lambda0 <- get.facload(yy0,xx0,l_sd=l_sd)
        }else{
          lambda0 <- 1
        }
        
        if (jj>1){
          L[jj,1:(jj-1)] <- lambda0
          L[jj,jj] <- 1
        }else if (jj==1){
          L[jj,jj] <- 1
        }
      }else{
        normalizer <- exp(0.5*S[,jj])
        yy0 <- (eps[,jj])/normalizer
        xx0 <- fac[,,drop=FALSE]/normalizer
        l_sd <- pr[jj,]
        lambda0 <- get.facload(yy0,xx0,l_sd=l_sd)
        L[jj,] <- lambda0
      }
    }
  }else{
    for(jj in 1:m){
      normalizer <- exp(0.5*S[,jj])
      yy0 <- (eps[,jj])/normalizer
      xx0 <- fac[,,drop=FALSE]/normalizer
      l_sd <- pr[jj,]
      lambda0 <- get.facload(yy0,xx0,l_sd=l_sd)
      L[jj,] <- lambda0
    }
  }
  return(L)
}

# sample the latent factors
get.factors <- function(e,S,H,L,q,t){
  F_raw <- matrix(0,t,q)
  for (tt in 1:t){
    normalizer <- exp(-S[tt,]/2)
    Lt <- L*normalizer
    yt <- e[tt,]*normalizer
    
    if (q==1) fac.varinv <-  1/H[tt] else fac.varinv <- diag(q)/H[tt]
    fac.Sigma <-  solve(crossprod(Lt)+fac.varinv)
    fac.mean <- fac.Sigma%*%crossprod(Lt,yt)
    
    F_temp <- try(fac.mean + t(chol(fac.Sigma)) %*% rnorm(q),silent=TRUE)
    if (is(F_temp,"try-error")) F_temp <- fac.mean + t(chol(fac.Sigma+diag(q)*1e-6)) %*% rnorm(q)
    F_raw[tt,] <- F_temp
  }
  return(F_raw)
}