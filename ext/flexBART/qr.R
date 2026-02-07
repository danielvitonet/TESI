get.hs <- function(bdraw,lambda.hs,nu.hs,tau.hs,zeta.hs){
  k <- length(bdraw)
  
  lambda.hs <- invgamma::rinvgamma(k,shape=1,rate=1/nu.hs+bdraw^2/(2*tau.hs))
  tau.hs <- invgamma::rinvgamma(1,shape=(k+1)/2,rate=1/zeta.hs+sum(bdraw^2/lambda.hs)/2)
  nu.hs <- invgamma::rinvgamma(k,shape=1,rate=1+1/lambda.hs)
  zeta.hs <- invgamma::rinvgamma(1,shape=1,rate=1+1/tau.hs)
  
  ret <- list("psi"=(lambda.hs*tau.hs),"lambda"=lambda.hs,"tau"=tau.hs,"nu"=nu.hs,"zeta"=zeta.hs)
  return(ret)
}


# function to estimate a quantile regression with several sparsity/homogeneity inducing priors
qr <- function(Y,X,X.out,nburn=1000,nsave=1000,thinfac=1,prior="flat",p=0.5,quiet=FALSE){
  require(MASS)
  ntot <- nburn + nsave
  
  nthin <- round(nsave/thinfac)
  thin.set <- floor(seq(nburn+1,ntot,length.out=nthin))
  in.thin <- 0
  
  K <- ncol(X)
  T <- NROW(Y)
  
  b_draw <- b_prior <- rep(0,K)
  omega_p <- rep(1,K)
  sig2_draw <- rep(1)
  var_draw <- rep(1,T)
  
  offset_hs <- 1e-8 # offsetting threshold for HS prior, for stability
  if(prior=="flat"){
    omega_p <- omega_p*10
  }else if(prior=="hs"){
    lambda <- nu <- matrix(1,K)
    zeta <- tau <- 1
  }
  
  # prior for scale parameter
  n0 <- 3/2
  s0 <- 0.3/2
  
  # quantile regression
  theta <- (1-2*p)/(p*(1-p))
  tau2 <- 2/(p*(1-p))
  v <- matrix(1,T)
  
  # sampling
  beta_store <- array(NA,dim=c(nthin,K))
  sig2_store <- array(NA,dim=c(nthin))
  omega_store <- array(NA,dim=c(nthin,K))
  fcst_store <- array(NA,dim=c(nthin,1))
  
  if(!quiet){
    pb <- txtProgressBar(min = 0, max = ntot, style = 3) #start progress bar
  }
  
  for(irep in 1:ntot){
    # sample regression coefficients
    norm <- 1/sqrt(var_draw)
    X_new <- X*norm
    y_new <- (Y-theta*v)*norm
    
    Vpo <- solve(crossprod(X_new) + diag(1/omega_p))
    bpo <- Vpo %*% (crossprod(X_new,y_new) + diag(1/omega_p)%*%b_prior)
    b_draw <- try(as.numeric(bpo + t(chol(Vpo)) %*% rnorm(K)),silent=TRUE)
    if(is(b_draw,"try-error")){
      b_draw <- try(mvrnorm(1,bpo,Vpo),silent=TRUE)
      if(is(b_draw,"try-error")){
        b_draw <- as.numeric(bpo)
      }
    }
    
    # sample auxiliary quantity v_p
    delta2 <- ((Y - X%*%as.numeric(b_draw))^2)/(tau2 * sig2_draw)
    gamma2 <- (2/sig2_draw) + ((theta^2) / (tau2 * sig2_draw))
    for(tt in 1:T){
      v[tt] <- GIGrvg::rgig(n=1,1/2,delta2[tt],gamma2)
    }
    
    # sample scale parameter sigma_p
    eps <- Y - X%*%as.numeric(b_draw) - theta*v
    ntilde <- n0 + 3*T
    stilde <- s0 + 2*sum(v) + sum(eps^2/(tau2 * v))
    sig2_draw <- 1/rgamma(1,ntilde/2,stilde/2)
    var_draw <- as.numeric(tau2 * sig2_draw * v)
    
    # common shrinkage
    if(prior=="hs"){
      hs_draw <- get.hs(bdraw=as.numeric(b_draw),
                        lambda.hs=as.numeric(lambda),nu.hs=as.numeric(nu),
                        tau.hs=as.numeric(tau),zeta.hs=as.numeric(zeta))
      omega_draw <- hs_draw$psi
      
      omega_draw[omega_draw<offset_hs] <- offset_hs
      omega_draw[omega_draw>10] <- 10
      
      omega_p[] <- omega_draw
      omega_p[K] <- 10
      lambda[] <- hs_draw$lambda
      nu[] <- hs_draw$nu
      tau <- hs_draw$tau
      zeta <- hs_draw$zeta
    }
    
    # storage
    if(irep %in% thin.set){
      in.thin <- in.thin+1
      beta_store[in.thin,] <- b_draw
      sig2_store[in.thin] <- sig2_draw
      fcst_store[in.thin] <- sum(X.out*b_draw)
    }
    
    if(!quiet) setTxtProgressBar(pb, irep)
  }
  
  # ret_list <- fcst_store
  ret_list <- list("b"=beta_store,"sig2"=sig2_store,"omega"=omega_store,"fcst"=fcst_store)
  return(ret_list)
}

