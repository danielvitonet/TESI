# ------------------------------------------------------------------------------------
# MCMC setup
flexBART <- function(Yraw, nburn=5000, nsave = 5000, thinfac = 1, prior="HS", prior.sig, model ="mixBART", sv = "SV", fc.approx="exact", restr.var = NULL, fhorz = 3, quiet = FALSE, pr.mean = prior.mean){
  # auxiliary functions
  #source("aux_func.R")

  
  # required packages
  require(mvtnorm)            # L   for mix$mixBART
  require(invgamma)           # L   for mix$mixBART
  require(factorstochvol)     # L   for mix$mixBART
  require(stochvol)           # R\L sampling stochastic volatilities

  Z.sv <- NULL
  
  # standardize data
  Ymu <- apply(Yraw, 2, mean,na.rm=T)
  Ysd <- apply(Yraw, 2, sd,na.rm=T)
  Yraw <- apply(Yraw, 2, function(x){(x-mean(x,na.rm=T))/sd(x,na.rm=T)})
  
  iter.update <- 250
  
  nthin <- round(thinfac * nsave)
  ntot <- nburn + nsave
  thin.set <- floor(seq(nburn+1,ntot,length.out=nthin))
  in.thin <- 0
  
  p <- 12 # lag order
  cons <- 0 # whether constant should be included
  id.L <- FALSE # whether factor model for the errors should be identified
  
  # create design matrices Y/X
  if(cons){
    X <- cbind(mlag(Yraw,p),1)[(p+1):nrow(Yraw),]
  }else{
    X <- cbind(mlag(Yraw,p))[(p+1):nrow(Yraw),]
  }
  Y <- Yraw[(p+1):nrow(Yraw), ] # leave original input matrix unchanged
  
  T <- nrow(Y)
  M <- ncol(Y)
  K <- ncol(X)
  
  
  # initialize draws
  PHI_draw  <- A_draw <- A_approx <- A_OLS <- ginv(crossprod(X))%*%crossprod(X,Y)
  A_prior <- matrix(0, K, M); A_prior[1:M, 1:M] <- pr.mean
  Sig_OLS <- crossprod(Y-X%*%A_OLS)/T
  
  # factor structure in the errors
  eta <- E.hat <- Y - X%*%PHI_draw
  Q <- ledermann(M)
  if (Q==1 || Q==0) Q  <- 2
  
  
  Lambda <- matrix(0,M,Q)
  Ft <- prcomp(eta)$x[,1:Q]
  Omega <- matrix(1,T,Q)
  theta.Lambda <- Lambda^0
  if(id.L){
    id.Lambda <- rbind(lower.tri(diag(Q)),matrix(TRUE,M-Q,Q))
    Lambda[1:Q,1:Q] <- diag(Q)
  }else{
    id.Lambda <- matrix(TRUE,M,Q)
  }
  ASV_approx <- matrix(0,K+1,M)
  AFSV_approx <- matrix(0,K+1,Q)
  
  # HS prior for factor loadings (by column)
  lambda.Lambda <- matrix(1,M,Q)
  nu.Lambda <- matrix(1,M,Q)
  zeta.Lambda <- tau.Lambda <- rep(1,Q)
  
  # covariance related objects
  eta <- matrix(NA,T,M)
  H <- matrix(-3,T,M)
  
  # objects for SV
  sv_draw <- list()
  sv_latent <- list()
  for (mm in seq_len(M)){
    sv_draw[[mm]] <- list(mu = 0, phi = 0.99, sigma = 0.01, nu = Inf, rho = 0, beta = NA, latent0 = 0)
    sv_latent[[mm]] <- rep(0,T)
  }
  
  svfac_draw <- list()
  svfac_latent <- list()
  for (qq in seq_len(Q)){
    svfac_draw[[qq]] <- list(mu = 0, phi = 0.99, sigma = 0.01, nu = Inf, rho = 0, beta = NA, latent0 = 0)
    svfac_latent[[qq]] <- rep(0,T)
  }
  
  sv_priors <- list()
  if(sv == "SV"){
    for(mm in 1:M){
      sv_priors[[mm]] <- specify_priors(
        mu = sv_normal(mean = 0, sd = 1),
        phi = sv_beta(shape1 = 5, shape2 = 1.5),
        sigma2 = sv_gamma(shape = 0.5, rate = 1/(2*0.01)),
        nu = sv_infinity(),
        rho = sv_constant(0)
      ) 
    }
  }
  
  svfac_priors <- list()
  if(sv == "SV"){
    for(qq in 1:Q){
      svfac_priors[[qq]] <- specify_priors(
        mu = sv_normal(mean = 0, sd = 1),
        phi = sv_beta(shape1 = 5, shape2 = 1.5),
        sigma2 = sv_gamma(shape = 0.5, rate = 1/(2*0.01)),
        nu = sv_infinity(),
        rho = sv_constant(0)
      ) 
    }
  }
  
  # priors
  theta_A <- matrix(10,K,M)
  a0sig <- 3
  b0sig <- 0.3
  
  if(prior=="Minn"){
    # MH preliminaries
    accept <- 0
    accept2 <- 0
    scale1 <- .43
    scale2 <- .45
    
    shrink.1 <- 0.5
    shrink.2 <- 0.5
    
    ind <- matrix(0,M,p)
    for (i in 1:M){
      ind[i,] <- seq(i,K-cons, by=M)
    }
    
    sigma_sq  <- matrix(0,M,1) #vector which stores the residual variance
    for (i in 1:M){
      Ylag_i <- mlag(Yraw[,i],p)
      Ylag_i <- Ylag_i[(p+1):nrow(Ylag_i),]
      Y_i <- Yraw[(p+1):nrow(Yraw),i]
      alpha_i <- solve(t(Ylag_i)%*%Ylag_i)%*%t(Ylag_i)%*%Y_i
      sigma_sq[i,1] <- (1/(nrow(Yraw)-p))*t(Y_i-Ylag_i%*%alpha_i)%*%(Y_i-Ylag_i%*%alpha_i) 
    }
    
    V_i <- get_V(shrink.1,shrink.2,ind1=ind,sigma_sq1=sigma_sq,p1=p,set.cons=cons,M1=M,K1=K)
    theta_PHI <- theta_A <- V_i
  }
  
  # BART initialization
  cgm.level <- 0.95
  cgm.exp <- 2
  sd.mu <- 2
  num.trees <- 250
  
  control <- dbartsControl(verbose = FALSE, keepTrainingFits = TRUE, useQuantiles = FALSE,
                           keepTrees = FALSE, n.samples = ntot,
                           n.cuts = 100L, n.burn = nburn, n.trees = num.trees, n.chains = 1,
                           n.threads = 1, n.thin = 1L, printEvery = 1,
                           printCutoffs = 0L, rngKind = "default", rngNormalKind = "default",
                           updateState = FALSE)
  
  sampler.list <- list()
  sampler.hetero.list <- list()
  sampler.fac.list <- list()
  svdraw.list <- list()
  
  for (jj in seq_len(M)){
    if (sv == "heteroBART" || sv == "SV"){ 
      prior.sig = c(10000^50, 0.5)
      sigma.init <- 1
    }else{
      sigma.init <- sqrt(Sig_OLS[jj,jj]) 
    }
    
    sampler.list[[jj]] <- dbarts(Y[,jj]~X, control = control,tree.prior = cgm(cgm.exp, cgm.level), node.prior = normal(sd.mu), n.samples = nsave, weights=rep(1,T), sigma=sigma.init, resid.prior = chisq(prior.sig[[1]], prior.sig[[2]]))
    if (sv == "heteroBART"){
      sampler.hetero.list[[jj]] <-  dbarts(E.hat[, jj] ~ cbind(seq(1, T), X), control = control,tree.prior = cgm(cgm.exp, cgm.level), node.prior = normal(sd.mu), n.samples = nsave, weights=rep(1,T), sigma=3)
    }
  }
  
  for (qq in seq_len(Q)){
    if (sv == "heteroBART" || sv == "SV"){ 
      prior.sig = c(10000^50, 0.5)
      sigma.init <- 1
    }else{
      sigma.init <- 1
    }
    
    if (sv == "heteroBART"){
      sampler.fac.list[[qq]] <-  dbarts(Ft[,qq] ~ cbind(seq(1, T), X), control = control,tree.prior = cgm(cgm.exp, cgm.level), node.prior = normal(sd.mu), n.samples = nsave, weights=rep(1,T), sigma=3)
    }
  }
  
  sampler.run <- list()
  sampler.sv.run <- list()
  sampler.fsv.run <- list()
  if (sv == "heteroBART"){
    # Stuff necessary for the 10-component mixture)
    m_st  = c(1.92677, 1.34744, 0.73504, 0.02266, -0.85173, -1.97278, -3.46788, -5.55246, -8.68384, -14.65000)
    v_st2 = c(0.11265, 0.17788, 0.26768, 0.40611,  0.62699,  0.98583,  1.57469,  2.54498,  4.16591,   7.33342)
    q_st  = c(0.00609, 0.04775, 0.13057, 0.20674,  0.22715,  0.18842,  0.12047,  0.05591,  0.01575,   0.00115)
  }
  
  sigma.fac.mat <- matrix(1,Q,1)
  sigma.mat <- matrix(1, M, 1)
  count.mat <- matrix(0, M*p, M)
  
  # initialize HS prior on linear VAR coefs
  lambda.VAR <- 1
  nu.VAR <- 1
  tau.VAR <- 1
  zeta.VAR <- 1
  sv_parms_mat <- matrix(NA, M, 4)
  svfac_parms_mat <- matrix(NA, Q, 4)
  Sig_t <- array(0,dim=c(T,M,M))
  
  Y.fit.BART <- Y.fit.VAR <- Y.fit.full <- Y*0
  
  # storage objects
  if (!is.null(restr.var)){
    sl.quants <- c(0, 0.005, 0.01, seq(0.05, 1, by=0.05))
    range.conditional <- quantile(Y[,restr.var], sl.quants)
    R <- length(range.conditional)
    H_store <- array(NA, dim=c(nthin, T, M)) #Store stochastic volatilities
    fcst_store <- Hfcst_store <- array(NA,dim=c(nthin,fhorz,M,R))
  }else{
    H_store <- array(NA, dim=c(nthin, T, M)) #Store stochastic volatilities
    fcst_store <- Hfcst_store <- array(NA,dim=c(nthin,fhorz,M))
  }
  
  # -----------------------------------------------------------------------------
  # start Gibbs sampler
  if(!quiet){
    pb <- txtProgressBar(min = 0, max = ntot, style = 3)
    start  <- Sys.time()
  }
  
  error.full <- eta - X%*%PHI_draw
  X.ginv <- MASS::ginv(X)
  F.ginv <- MASS::ginv(cbind(seq(1, T), X))
  count.dirt <- 0
  for (irep in 1:ntot){
    check <- TRUE # algorithm stability
    
    # 1) sample model coefficients (either linear VAR or BART)
    Y_ <- Y - tcrossprod(Ft,Lambda)
    while (check){
      for (mm in seq_len(M)){
        if (model == "mixBART" || model == "BART"){
          sampler.list[[mm]]$setResponse(Y_[,mm])
          
          # This part estimates the eq-specific BART model
          rep_mm <- sampler.list[[mm]]$run(0L, 1L)
          sampler.run[[mm]] <- rep_mm
          
          sig.draw <- rep_mm$sigma
          sigma.mat[mm,] <- sig.draw
          A_approx[,mm] <- X.ginv%*%rep_mm$train
          
          if (any(is.na(rep_mm$train))) break
          Y.fit.BART[,mm] <- rep_mm$train
          eta[,mm] <- Y_[,mm] - rep_mm$train
          A_draw[,mm] <- X.ginv%*%rep_mm$train
          count.mat[,mm] <- rep_mm$varcount
        }else{
          A_draw <- matrix(0, K, M)
          eta <- Y_
        }
        
        # Here we simulate the regression part of the model
        norm_mm <- as.numeric(exp(-.5*H[, mm]))
        eta_mm <- eta[,mm] * norm_mm
        
        if (model=="mixBART" || model=="BVAR"){
          #This part allows for estimating a VAR model with BART shocks that depend on the lagged values of Y as well
          u_mm <- X*norm_mm
          theta_PHI <- theta_A[,mm]
          v0.inv <- diag(1/c(theta_PHI)) 
          
          V.cov <- solve(crossprod(u_mm) + v0.inv)
          mu.cov <- V.cov %*% (crossprod(u_mm, eta_mm) + v0.inv%*%c(A_prior[, mm]))
          
          mu.phi.a0.draw <- try(mu.cov + t(chol(V.cov)) %*% rnorm(ncol(V.cov)), silent=TRUE)
          if (is(mu.phi.a0.draw, "try-error")){
            mu.phi.a0.draw <-  as.numeric(mvtnorm::rmvnorm(1, mu.cov, as.matrix(Matrix::forceSymmetric((V.cov)))))
          } 
          
          error.full[,mm] <- eta[,mm]
          PHI_draw[,mm] <- mu.phi.a0.draw[1:K]
        }else{
          PHI_draw <- matrix(0, K, M)
        }
      }
      
      max.eigen <- (max(abs(eigen(get_companion(PHI_draw, c(M, 0, p))$MM)$values)))
      if (irep > 0.05*nburn){
        if (max.eigen < 3) check <- FALSE
      }else{
        check <- FALSE
      }
    }
    Y.fit.VAR <- X%*%PHI_draw
    Y.fit.full <- Y.fit.BART + Y.fit.VAR
    
    shock_norm <- Y - Y.fit.full - tcrossprod(Ft,Lambda)
    for (mm in seq_len(M)){
      if (sv == "SV"){
        svdraw_mm <- svsample_general_cpp(shock_norm[,mm]/sigma.mat[mm], startpara = sv_draw[[mm]], startlatent = sv_latent[[mm]], priorspec = sv_priors[[mm]])
        sv_draw[[mm]][c("mu", "phi", "sigma")] <- as.list(svdraw_mm$para[, c("mu", "phi", "sigma")])
        sv_latent[[mm]] <- svdraw_mm$latent
        sv_parms_mat[mm, ] <- c(svdraw_mm$para[, c("mu", "phi", "sigma")], svdraw_mm$latent[T])
        
        normalizer <- as.numeric(exp(-.5*svdraw_mm$latent))
        weights.new <- as.numeric(exp(-svdraw_mm$latent))
        sampler.list[[mm]]$setWeights(weights.new)
        H[,mm] <- log(sigma.mat[mm]^2) + svdraw_mm$latent
      }else if (sv =="heteroBART"){
        eta.star <- shock_norm[ , mm]/sigma.mat[mm]
        eta.off <- 1e-15 
        eta.star <- log((eta.star - eta.off)^2)
        
        # Sample the mixture components (previous draw of ht)
        z <- sapply(eta.star-H[,mm], ncind, m_st, sqrt(v_st2), q_st)
        
        # Subset mean and variances to the sampled mixture components; (n x p) matrices
        m_st_all = matrix(m_st[z],T, 1)
        v_st2_all = matrix(v_st2[z], T, 1)
        eta.star <- matrix((eta.star - m_st_all))
        
        # Do BART for conditional variances
        sampler.hetero.list[[mm]]$setResponse(eta.star)
        sampler.hetero.list[[mm]]$setWeights(1/v_st2_all)
        
        rep.hetero.i <- sampler.hetero.list[[mm]]$run(0L, 1L)
        sampler.sv.run[[mm]] <- rep.hetero.i
        
        ht.i <- rep.hetero.i$train
        ht.i[ht.i < -10] <- -10
        ASV_approx[,mm] <- F.ginv%*%ht.i
        
        weights.new <- as.numeric(exp(-ht.i))
        sampler.list[[mm]]$setWeights(weights.new)
        H[ , mm] <- ht.i + log(as.numeric(sigma.mat[mm])^2)
      }else{
        H[,mm] <- log(sigma.mat[mm]^2)
        if(model=="BVAR"){
          H[,mm] <- log(1/rgamma(1,0.01+T/2,0.01+sum(shock_norm[,mm]^2)))
        }
      }
    }
    H[H<log(1e-6)] <- log(1e-6)

    # draw factors and loadings
    eps <- Y - Y.fit.full
    Ft <- get.factors(eps,S=(H),H=exp(Omega),L=Lambda,q=Q,t=T)
    if(!id.L) Ft <- apply(Ft,2,function(x) (x-mean(x))/sd(x)) # normalize factor draw
    Lambda <- get.Lambda(eps,fac=Ft,S=(H),pr=theta.Lambda,m=M,q=Q,id.fac=!id.Lambda[1,Q])
    
    # shrinkage on loadings by columns
    for(qq in 1:Q){
      hs_draw <- get.hs(bdraw=as.numeric(Lambda[id.Lambda[,qq],qq]),
                        lambda.hs=as.numeric(lambda.Lambda[id.Lambda[,qq],qq]),
                        nu.hs=nu.Lambda[id.Lambda[,qq],qq],tau.hs=tau.Lambda[qq],zeta.hs=zeta.Lambda[qq])
      theta.Lambda[id.Lambda[,qq],qq] <- hs_draw$psi
      lambda.Lambda[id.Lambda[,qq],qq] <- hs_draw$lambda
      nu.Lambda[id.Lambda[,qq],qq] <- hs_draw$nu
      tau.Lambda[qq] <- hs_draw$tau
      zeta.Lambda[qq] <- hs_draw$zeta
    }
    theta.Lambda[theta.Lambda<1e-5] <- 1e-5
    
    # draw variances of factors
    if(sv=="SV"){
      for(qq in 1:Q){
        svdraw_qq <- svsample_general_cpp(Ft[,qq]/sigma.fac.mat[qq], startpara = svfac_draw[[qq]], 
                                          startlatent = svfac_latent[[qq]], priorspec = svfac_priors[[qq]])
        svfac_draw[[qq]][c("mu", "phi", "sigma")] <- as.list(svdraw_qq$para[, c("mu", "phi", "sigma")])
        svfac_latent[[qq]] <- svdraw_qq$latent
        svfac_parms_mat[qq, ] <- c(svdraw_qq$para[, c("mu", "phi", "sigma")], svdraw_qq$latent[T])
        Omega[,qq] <- log(sigma.fac.mat[qq]^2) + svdraw_qq$latent
      }
    }else if(sv=="heteroBART"){
      for(qq in 1:Q){
        eta.star <- Ft[,qq]/sigma.fac.mat[qq]
        eta.off <- 1e-15 
        eta.star <- log((eta.star - eta.off)^2)
        
        # Sample the mixture components (previous draw of ht)
        z <- sapply(eta.star-Omega[,qq],ncind,m_st,sqrt(v_st2),q_st)
        
        # Subset mean and variances to the sampled mixture components; (n x p) matrices
        m_st_all = matrix(m_st[z],T,1)
        v_st2_all = matrix(v_st2[z],T,1)
        eta.star <- matrix((eta.star - m_st_all))
        
        # Do BART for conditional variances
        sampler.fac.list[[qq]]$setResponse(eta.star)
        sampler.fac.list[[qq]]$setWeights(1/v_st2_all)
        
        rep.fac.q <- sampler.fac.list[[qq]]$run(0L, 1L)
        sampler.fsv.run[[qq]] <- rep.fac.q
        
        omega.q <- rep.fac.q$train
        omega.q[omega.q < -8] <- -8
        AFSV_approx[,qq] <- F.ginv%*%omega.q
        Omega[,qq] <- omega.q + log(as.numeric(sigma.fac.mat[qq])^2)
      }
    }else{
      for(qq in 1:Q){
        Omega[,qq] <- 0 # log(1/rgamma(1,3+T/2,0.3+sum(Ft[,qq]^2)))
      }
    }
    
    for(tt in seq_len(T)){
      if (Q>1) S.t <- Lambda%*%tcrossprod(diag(exp(Omega[tt,])),Lambda) + diag(exp(H[tt,])) else S.t <- exp(Omega[tt,])*Lambda%*%t(Lambda) + diag(exp(H[tt,]))
      Sig_t[tt,,] <- S.t
    }
    
    # 2) updating shrinkage priors
    if (model=="mixBART" || model=="BVAR"){
      #Update HS for the VAR coefficients
      if(prior=="HS"){
        hs_draw <- get.hs(bdraw=as.numeric(PHI_draw - A_prior),lambda.hs=lambda.VAR,nu.hs=nu.VAR,tau.hs=tau.VAR,zeta.hs=zeta.VAR)
        lambda.VAR <- hs_draw$lambda
        nu.VAR <- hs_draw$nu
        tau.VAR <- hs_draw$tau
        zeta.VAR <- hs_draw$zeta
        theta_A <- matrix(hs_draw$psi,K,M)
      }else if(prior=="Minn"){
        shrink.1.prop <- exp(rnorm(1,0,scale1))*shrink.1
        theta.prop <- get_V(shrink.1.prop,shrink.2,ind1=ind,sigma_sq1=sigma_sq,p1=p,set.cons=cons,M1=M,K1=K)
        
        post1.prop <- sum(dnorm(as.vector(PHI_draw),as.vector(A_prior),sqrt(as.vector(theta.prop)),log=TRUE)) + dgamma(shrink.1.prop,0.01,0.01,log=TRUE)
        post1 <- sum(dnorm(as.vector(PHI_draw),as.vector(A_prior),sqrt(as.vector(theta_A)),log=TRUE)) + dgamma(shrink.1,.01,0.01,log=TRUE)
        if ((post1.prop-post1)>log(runif(1,0,1))){
          shrink.1 <- shrink.1.prop
          theta_A <- theta.prop
          accept <- accept+1
        }
        
        shrink.2.prop <- exp(rnorm(1,0,scale2))*shrink.2
        theta.prop <- get_V(shrink.1,shrink.2.prop,ind1=ind,sigma_sq1=sigma_sq,p1=p,set.cons=cons,M1=M,K1=K)
        
        post2.prop <- sum(dnorm(as.vector(PHI_draw),as.vector(A_prior),sqrt(as.vector(theta.prop)),log=TRUE))+dgamma(shrink.2.prop,0.01,0.01,log=TRUE)
        post2 <- sum(dnorm(as.vector(PHI_draw),as.vector(A_prior),sqrt(as.vector(theta_A)),log=TRUE))+dgamma(shrink.2,0.01,0.01,log=TRUE)
        if ((post2.prop-post2)>log(runif(1,0,1))){
          shrink.2 <- shrink.2.prop
          theta_A <- theta.prop
          accept2 <- accept2+1
        }
        
        if (irep<(0.5*nburn)){
          if ((accept/irep)>0.3) scale1 <- 1.01*scale1
          if ((accept/irep)<0.15) scale1 <- 0.99*scale1
          if ((accept2/irep)>0.3) scale2 <- 1.01*scale2
          if ((accept2/irep)<0.15) scale2 <- 0.99*scale2
        }
      }
      theta_A[theta_A < 1e-12] <- 1e-12
    }
    
    if(irep %in% thin.set){
      in.thin <- in.thin+1
      H_store[in.thin,,] <- H
      
      # include forecast loop here
      if(fhorz>0 & is.null(restr.var)){
        Yfc <- matrix(NA,fhorz,M)
        Hfc <- matrix(NA,fhorz,M)
        
        if (cons){
          X.hat <- c(Y[T,],X[T,1:(M*(p-1))],1)
        }else{
          X.hat <- c(Y[T,],X[T,1:(M*(p-1))])  
        }
        
        if (sv == "SV" || sv == "heteroBART"){
          HT <-  H[T, ] - log(as.numeric(sigma.mat)^2)
          OT <- Omega[T, ] - log(as.numeric(sigma.fac.mat)^2)
        }else{
          HT <- H[T,]
          OT <- Omega[T, ]
        }
        
        Sig_T <- Sig_t[T,,] # use final observation for Sigma
        tree.pred <- matrix(0, M)
        for (hh in seq_len(fhorz)){
          if (sv == "SV"){
            HT  <- log(as.numeric(sigma.mat)^2) + (sv_parms_mat[, 1] + sv_parms_mat[ , 2] * (HT - sv_parms_mat[,1]) + sv_parms_mat[ , 3]*rnorm(M))
            OT <- log(as.numeric(sigma.fac.mat)^2) + (svfac_parms_mat[,1] + svfac_parms_mat[,2] * (OT - svfac_parms_mat[,1]) + svfac_parms_mat[,3]*rnorm(Q))
          }else if (sv == "heteroBART"){
            vola.predict.tree <- matrix(0, M, 1)
            for (nn in seq_len(M)){
              vola.predict.tree[nn,] <- sampler.hetero.list[[nn]]$predict(c(T+hh, X.hat))#
            }
            HT <- as.numeric(log(as.numeric(sigma.mat)^2) + vola.predict.tree)
            
            vola.predict.fac <- matrix(0, Q, 1)
            for (qq in seq_len(Q)){
              vola.predict.fac[qq,] <- sampler.fac.list[[qq]]$predict(c(T+hh,X.hat))
            }
            OT <- as.numeric(log(as.numeric(sigma.fac.mat)^2) + vola.predict.fac)
          }
          Hfc[hh,] <- exp(HT)
          
          if (model == "mixBART" || model == "BART"){
            for (j in seq_len(M)) tree.pred[j] <- sampler.list[[j]]$predict(X.hat)
          }else{
            tree.pred <- rep(0, M)
          }
          Sig_T <- Lambda%*%tcrossprod(diag(exp(OT)),Lambda) + diag(exp(HT))
          
          if(fc.approx=="exact"){
            Y.tp1 <- try(as.numeric(X.hat%*%PHI_draw)+ as.numeric(tree.pred) + t(chol(Sig_T))%*%rnorm(M), silent=TRUE)
            if (is(Y.tp1, "try-error")) Y.tp1 <- as.numeric(mvtnorm::rmvnorm(1, tree.pred + as.numeric(X.hat%*%PHI_draw), Sig_T))
          }else if(fc.approx=="approx"){
            Y.tp1 <- try(as.numeric(X.hat%*%(A_approx+PHI_draw)) + t(chol(Sig_T))%*%rnorm(M), silent=TRUE)
            if (is(Y.tp1, "try-error")) Y.tp1 <- as.numeric(mvtnorm::rmvnorm(1, as.numeric(X.hat%*%(A_approx+PHI_draw)), Sig_T))
          }
          
          if (cons){
            X.hat <- c(Y.tp1, X.hat[1:(M*(p-1))],1)
          }else{
            X.hat <- c(Y.tp1, X.hat[1:(M*(p-1))])
          }
          Yfc[hh,] <- Y.tp1
        }
        fcst_store[in.thin,,] <- (Yfc*t(matrix(Ysd,M,fhorz)))+t(matrix(Ymu,M,fhorz))
        Hfcst_store[in.thin,,] <- Hfc*t(matrix(Ysd,M,fhorz))
      }else if(fhorz>0 & !is.null(restr.var)){
        Yfc <- array(NA,c(fhorz,M, R))
        Hfc <- array(NA,c(fhorz,M, R))
        sl.cond <- which(colnames(Y)==restr.var)
        
        for (r in seq_len(R)){
          cond.vals <- range.conditional[[r]]
          sl.cond.lags <- sl.cond#seq(sl.cond, K, by=M) /only first lag
          if (cons){
            X.hat <- c(Y[T,],X[T,1:(M*(p-1))],1)
          }else{
            X.hat <- c(Y[T,],X[T,1:(M*(p-1))])  
          }
          X.hat[sl.cond.lags] <- cond.vals
          
          if (sv == "SV" || sv == "heteroBART"){
            HT <-  H[T, ] - log(as.numeric(sigma.mat)^2)
            OT <- Omega[T, ] - log(as.numeric(sigma.fac.mat)^2)
          }else{
            HT <- H[T,]
            OT <- Omega[T, ]
          }
          
          Sig_T <- Sig_t[T,,] # use final observation for Sigma
          tree.pred <- matrix(0, M)
          for (hh in seq_len(fhorz)){
            if (sv == "SV"){
              HT  <- log(as.numeric(sigma.mat)^2) + (sv_parms_mat[, 1] + sv_parms_mat[ , 2] * (HT - sv_parms_mat[,1]) + sv_parms_mat[ , 3]*rnorm(M))
              OT <- log(as.numeric(sigma.fac.mat)^2) + (svfac_parms_mat[,1] + svfac_parms_mat[,2] * (OT - svfac_parms_mat[,1]) + svfac_parms_mat[,3]*rnorm(Q))
            }else if (sv == "heteroBART"){
              vola.predict.tree <- matrix(0, M, 1)
              for (nn in seq_len(M)){
                vola.predict.tree[nn,] <- sampler.hetero.list[[nn]]$predict(c(T+hh, X.hat))
              }
              HT <- as.numeric(log(as.numeric(sigma.mat)^2) + vola.predict.tree)
              
              vola.predict.fac <- matrix(0, Q, 1)
              for (qq in seq_len(Q)){
                vola.predict.fac[qq,] <- sampler.fac.list[[qq]]$predict(c(T+hh,X.hat))
              }
              OT <- as.numeric(log(as.numeric(sigma.fac.mat)^2) + vola.predict.fac)
            }
            Hfc[hh,,r] <- exp(HT)*Ysd
            
            if (model == "mixBART" || model == "BART"){
              for (j in seq_len(M)) tree.pred[j] <- sampler.list[[j]]$predict(X.hat)
              tree.pred[sl.cond] <- cond.vals
            }else{
              tree.pred <- rep(0, M)
            }
            Sig_T <- Lambda%*%tcrossprod(diag(exp(OT)),Lambda) + diag(exp(HT))
            Y.tp1 <- try(as.numeric(X.hat%*%PHI_draw)+ as.numeric(tree.pred) + t(chol(Sig_T))%*%rnorm(M), silent=TRUE)
            if (is(Y.tp1, "try-error")) Y.tp1 <- as.numeric(mvtnorm::rmvnorm(1, tree.pred + as.numeric(X.hat%*%PHI_draw), Sig_T))
             if (hh <= 6) Y.tp1[sl.cond] <- cond.vals #only for six months
            if (cons){
              X.hat <- c(Y.tp1, X.hat[1:(M*(p-1))],1)
            }else{
              X.hat <- c(Y.tp1, X.hat[1:(M*(p-1))])
            }
            Yfc[hh,,r] <- Y.tp1*Ysd + Ymu
          }
        }
        fcst_store[in.thin,,,] <- Yfc
        Hfcst_store[in.thin,,,] <- Hfc
      }
    }
    
    if(!quiet){
      setTxtProgressBar(pb, irep)
      if (irep %% iter.update==0){
        end <- Sys.time()
        message(paste0("\n Average time for single draw over last ",iter.update," draws ", round(as.numeric(end-start)/iter.update, digits=4), " seconds, currently at draw ", irep))
        start <- Sys.time()
        ts.plot(cbind(Y.fit.full, Y), col=c(rep(1, M), rep(2, M)))
      }
    }
  }
  dimnames(fcst_store) <- dimnames(Hfcst_store) <- list(paste0("mcmc",1:nthin),paste0("fhorz", 1:fhorz),colnames(Y))
  return_obj <- list("fcst"=fcst_store,"Hfcst"=Hfcst_store)
  return(return_obj)
}
