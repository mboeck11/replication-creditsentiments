###############################################
### Auxiliary functions                     ###
### The Impact of Credit Market Sentiments  ###
### Maximilian Boeck                        ###
### 30/04/2021                              ###
###############################################

ar_svf <- function(Yraw, p=1, nburn=500, nsave=500) {
  library(stochvol)
  Xraw <- mlag(Yraw, 1)
  y    <- Yraw[-1,,drop=F]
  X    <- Xraw[-1,,drop=F]
  ntot  <- nburn+nsave
  T <- nrow(y)
  K <- ncol(X)
  
  # initial draws
  svdraw <- list(para = c(mu = -10, phi = 0.9, sigma = 0.2, latent0=-10),
                 latent = rep(-10, T))
  beta.draw <- 0
  
  # prior values
  b0    <- 0
  B0inv <- diag(1)/10
  
  priormu <- c(0, 100)
  priorphi <- c(20, 1.5)
  Bsigma <- 1
  
  Sv_priors <- specify_priors(mu=sv_normal(mean=priormu[1], sd=priormu[2]), 
                              phi=sv_beta(priorphi[1],priorphi[2]), 
                              sigma2=sv_gamma(shape=0.5,rate=1/(2*Bsigma)))
  
  # container
  beta.store <- array(NA, c(nsave, K))
  sv.store   <- array(NA, c(nsave, T))
  
  for(irep in 1:ntot){
    ytilde <- y - X %*% beta.draw
    para   <- as.list(svdraw[["para"]])
    para$nu = Inf; para$rho=0; para$beta<-0
    temp <- svsample_fast_cpp(y=ytilde, draws=1, burnin=0, designmatrix=matrix(NA_real_),
                                priorspec=Sv_priors, thinpara=1, thinlatent=1, keeptime="all",
                                startpara=para, startlatent=svdraw[["latent"]],
                                keeptau=FALSE, print_settings=list(quiet=TRUE, n_chains=1, chain=1),
                                correct_model_misspecification=FALSE, interweave=TRUE, myoffset=0,
                                fast_sv=default_fast_sv)
    svdraw[["para"]] <- c(temp$para[1,"mu"],temp$para[1,"phi"],temp$para[1,"sigma"],temp$latent0[1,"h_0"])
    names(svdraw[["para"]])<-c("mu","phi","sigma","latent0")
    svdraw[["latent"]] <- temp$latent
    normalizer <- as.numeric(exp(-svdraw[["latent"]] / 2))
    Xnew <- X*normalizer
    ynew <- y*normalizer
    Sigma <- solve(crossprod(Xnew) + B0inv)
    mu <- Sigma %*% (crossprod(Xnew, ynew) + B0inv %*% b0)
    beta.draw <- as.numeric(mvtnorm::rmvnorm(1, mu, Sigma))
    if(irep>nburn){
      beta.store[irep-nburn,] <- beta.draw
      sv.store[irep-nburn,]   <- svdraw[["latent"]]
    }
    #if(irep%%100==0) print(irep)
  }
  ar1 <- apply(beta.store,2,mean)
  sv  <- exp(apply(sv.store,2,mean))
  
  return(list(ar=ar1,
              sv=sv))
}

#-----------------------------------------------------------------------------------------------

btvar_wishart_ng = function(Yraw, Zraw, plag, args) {
  #----------------------------------------INPUTS----------------------------------------------------#
  # prepare arguments
  draws = burnin = 5000; d=NULL; cons=TRUE; trend=FALSE; qtrend=FALSE; thin=1; Ex=NULL; eigen=FALSE; save.prior=FALSE
  if(!is.null(args)){
    for(aa in c("draws","burnin","d","cons","trend","qtrend","thin","Ex","eigen","save.prior")){
      if(aa%in%names(args)) assign(aa, args[[aa]])
    }
  }
  arglist=list(Yraw=Yraw, Zraw=Zraw, plag=plag, draws=draws, burnin=burnin, d=d, cons=cons, trend=trend,
               qtrend=qtrend, thin=thin, Ex=Ex, eigen=eigen, save.prior=save.prior)
  
  #---------------------------load Packages------------------------------------------#
  require(stochvol,quietly=TRUE)
  require(GIGrvg, quietly=TRUE)
  require(Rcpp, quietly = TRUE)
  require(LaplacesDemon, quietly=TRUE)
  sourceCpp("./scripts/do_rgig.cpp")
  
  #--------------------------------------------------------------------
  # regime stuff
  h    <- 2
  if(is.null(d)){
    d.min <- d.max <- 0
  }else{
    if(length(d) == 1){
      d.min <- d.max <- d
    }else if(length(d) == 2){
      d.min <- d[1]
      d.max <- d[2]
    }
  }
  
  Traw <- nrow(Yraw)
  n    <- ncol(Yraw)
  varNames <- colnames(Yraw)
  if(is.null(varNames)) varNames <- paste0("Y",seq(n))
  K     <- n*plag
  Ylag  <- mlag(Yraw,plag)
  varNameslags <- NULL
  for(pp in 1:plag) varNameslags <- c(varNameslags,paste(varNames,".lag",pp,sep=""))
  colnames(Ylag) <- varNameslags
  if(d.min == 0){
    if(d.max > 0) Zlag <- cbind(Zraw,mlag(Zraw,d.max)) else Zlag <- Zraw
  }else{
    Zlag <- mlag(Zraw,d.max)
  }
  colnames(Zlag)<-paste0("thrsh.lag",seq(d.min,d.max))
  
  X <- Ylag[(max(plag,d.max)+1):Traw,,drop=FALSE]
  Y <- Yraw[(max(plag,d.max)+1):Traw,,drop=FALSE]
  Z <- Zlag[(max(plag,d.max)+1):Traw,1,drop=FALSE]
  bigT <- nrow(Y)
  
  if(cons){
    X <- cbind(X,1)
    varNameslags <- c(varNameslags,"cons")
    colnames(X)[ncol(X)] <- "cons"
  }
  if(trend){
    X <- cbind(X,seq(1,bigT))
    varNameslags <- c(varNameslags,"trend")
    colnames(X)[ncol(X)] <- "trend"
  }
  if(qtrend){
    X <- cbind(X,seq(1,bigT)^2)
    varNameslags <- c(varNameslags,"qtrend")
    colnames(X)[ncol(X)] <- "qtrend"
  }
  
  k <- ncol(X)
  q <- k*n
  v <- (n*(n-1))/2
  
  #-----------------------------------------Initialize Gibbs sampler--------------#
  # OLS quantities and initial values
  A_OLS      <- array(NA, c(k, n, h))
  res_OLS    <- array(NA, c(bigT, n))
  SIGMA_OLS  <- array(NA, c(n, n, h))
  thrsh_mean <- mean(Z)
  S_mean <- matrix(1, bigT, 1)
  S_mean[Z>thrsh_mean,1] <- 2
  for(hh in 1:h){
    Y.s <- Y[S_mean==hh,,drop=FALSE]
    X.s <- X[S_mean==hh,,drop=FALSE]
    
    A_OLS[,,hh]          <- solve(crossprod(X.s)) %*% crossprod(X.s,Y.s)
    res_OLS[S_mean==hh,] <- Y.s - X.s %*% A_OLS[,,hh]
    SIGMA_OLS[,,hh]      <- crossprod(Y.s - X.s %*% A_OLS[,,hh])/(sum(S_mean==hh)-k)
  }
  #---------------------------------------------------------------------------------------------------------
  # Initial Values
  #---------------------------------------------------------------------------------------------------------
  A_draw     <- A_OLS
  SIGMA_draw <- array(SIGMA_OLS, c(n,n,h))
  res_draw   <- res_OLS
  gamma_draw <- thrsh_mean
  S_draw     <- S_mean
  
  #---------------------------------------------------------------------------------------------------------
  # PRIORS
  #---------------------------------------------------------------------------------------------------------
  # Priors on VAR coefs
  #-----------------------------
  prmean    = 1
  # prior mean
  a_prior <- matrix(0, k, n) # state independent
  diag(a_prior) <- prmean
  # prior variance
  A_prior <- array(.1, c(k, n, h)) # state dependent
  
  # NG stuff
  d_lambda  <- 0.01
  e_lambda  <- 0.01
  tau_start <- 0.7
  sample_tau <- TRUE
  
  lambda2_draw <- matrix(0.01, plag ,h, dimnames=list(paste0("lag.",seq(1,plag)),c("regime1","regime2")))
  tau_draw     <- matrix(tau_start, plag, h, dimnames=dimnames(lambda2_draw))
  tau_tuning   <- matrix(.43, plag, h, dimnames=dimnames(tau_draw))
  tau_accept   <- matrix(0, plag, h, dimnames=dimnames(tau_draw))
  
  # VCV matrix
  # compute AR residuals of AR(4)
  arsig2 <- get_resid_ar(Yraw, 4)
  
  S_prior  <- diag(c(arsig2))
  nu_prior <- n + 2
  
  #---------threshold----------------#
  min.share        <- .1 # minimum share of observations per regime
  d_draw           <- d.min
  s_d              <- 2.4^2
  scale            <- 0.7
  g_prior          <- mean(Zraw)
  G_prior          <- 10
  scale_g          <- scale
  accept_g         <- 0
  Ncrit            <- min.share * bigT
  Smat_draw        <- crSmat(S_draw, h)
  # adaptive MH
  gamma_store_burn <- array(NA, c(burnin, 1))
  epsilon          <- 0.00000001
  eigs             <- rep(NA,h)
  
  #---------------------------------------------------------------------------------------------------------
  # SAMPLER MISCELLANEOUS
  #---------------------------------------------------------------------------------------------------------
  ntot  <- draws+burnin
  
  # thinning
  count <- 0
  thindraws    <- draws/thin
  thin.draws   <- seq(burnin+1,ntot,by=thin)
  arglist      <- c(arglist, thindraws=thindraws)
  
  #---------------------------------------------------------------------------------------------------------
  # STORAGES
  #---------------------------------------------------------------------------------------------------------
  A_store       <- array(NA_real_, c(thindraws, k, n, h))
  SIGMA_store   <- array(NA_real_, c(thindraws, n, n, h))
  res_store     <- array(NA_real_, c(thindraws, bigT, n))
  # NG
  if(save.prior){
    Aprior_store  <- array(NA_real_, c(thindraws, k, n, h))
    lambda2_store <- array(NA_real_, c(thindraws, plag, h))
    tau_store     <- array(NA_real_, c(thindraws, plag, h))
  }else{
    Aprior_store <- lambda2_store <- tau_store <- NULL
  }
  # Threshold
  S_store     <- array(NA_real_, c(thindraws, bigT, 1))
  gamma_store <- array(NA_real_, c(thindraws, 1))
  d_store     <- array(NA_real_, c(thindraws, 1))
  
  #---------------------------------------------------------------------------------------------------------
  # Gibbs Sampler
  #---------------------------------------------------------------------------------------------------------
  for(irep in 1:ntot) {
    #----------------------------------------------------------------------------
    # Step 1: Sample regime-dependent coefficients
    for(hh in 1:h) {
      sl.state <- which(Smat_draw[,hh]==1)
      Y_hh     <- Y[sl.state,,drop=FALSE]
      X_hh     <- X[sl.state,,drop=FALSE]
      A_hh     <- A_draw[,,hh]
      SIGMA_hh <- SIGMA_draw[,,hh]
      Apr_hh   <- A_prior[,,hh]
      
      #----------------------------------------------------------------------------
      # Step 1/1: Sample coefficients
      
      A_OLS_hh <- solve(crossprod(X_hh)) %*% crossprod(X_hh,Y_hh)
      SIGMAinv_hh <- solve(SIGMA_hh)
      V_post <- try(chol2inv(chol(diag(as.vector(1/Apr_hh)) + SIGMAinv_hh %x% crossprod(X_hh))),silent=TRUE)
      if(is(V_post,"try-error")) V_post = try(solve(diag(as.vector(1/Apr_hh)) + SIGMAinv_hh %x% crossprod(X_hh)), silent=TRUE)
      if(is(V_post,"try-error")) V_post = ginv(diag(as.vector(1/Apr_hh)) + SIGMAinv_hh %x% crossprod(X_hh))
      a_post <- V_post %*% (diag(as.vector(1/Apr_hh)) %*% as.vector(a_prior) + SIGMAinv_hh %x% crossprod(X_hh) %*% as.vector(A_OLS_hh))
      
      A_hh <- matrix(a_post + t(chol(V_post)) %*% rnorm(q),k,n)
      dimnames(A_hh)<-list(varNameslags, varNames)
      
      #----------------------------------------------------------------------------
      # Step 1/2: Normal-Gamma for endogenous variables
      for (pp in 1:plag){
        slct.i  <- grep(paste0("\\.lag",pp), rownames(A_hh))
        A.lag   <- A_hh[slct.i,,drop=FALSE]
        a.prior <- a_prior[slct.i,,drop=FALSE]
        A.prior <- Apr_hh[slct.i,,drop=FALSE]
        
        n.end = nrow(A.lag)
        if(pp==1){
          lambda2_draw[pp,hh] <- rgamma(n = 1,
                                        shape = d_lambda + tau_draw[pp,hh]*n*n.end,
                                        rate = e_lambda + 0.5*tau_draw[pp,hh]*sum(A.prior))
        }else{
          lambda2_draw[pp,hh] <- rgamma(n = 1,
                                        shape = d_lambda + tau_draw[pp,hh]*n*n.end,
                                        rate = e_lambda + 0.5*tau_draw[pp,hh]*prod(lambda2_draw[1:(pp-1),hh])*sum(A.prior))
        }
        for(jj in 1:n.end){
          for(nn in 1:n){
            temp <- do_rgig1(lambda = tau_draw[pp,hh] - 0.5,
                             chi = (A.lag[jj,nn] - a.prior[jj,nn])^2,
                             psi = tau_draw[pp,hh]*prod(lambda2_draw[1:pp,hh]))
            A.prior[jj,nn] <- ifelse(temp<1e-8,1e-8,ifelse(temp>1e+8,1e+8,temp))
          }
        }
        Apr_hh[slct.i,] <- A.prior
        if (sample_tau){
          #Sample a_tau through a simple RWMH step (on-line tuning of the MH scaling within the first 50% of the burn-in phase)
          tau_prop      <- exp(rnorm(1,0,tau_tuning[pp,hh]))*tau_draw[pp,hh]
          post_tau_prop <- atau_post(atau=tau_prop, thetas=as.vector(A.prior), lambda2=prod(lambda2_draw[1:pp,hh]), k=length(A.prior))
          post_tau_old  <- atau_post(atau=tau_draw[pp,hh], thetas=as.vector(A.prior), lambda2=prod(lambda2_draw[1:pp,hh]), k=length(A.prior))
          post.diff     <- post_tau_prop-post_tau_old
          post.diff     <- ifelse(is.nan(post.diff),-Inf,post.diff)
          if (post.diff > log(runif(1,0,1))){
            tau_draw[pp,hh]   <- tau_prop
            tau_accept[pp,hh] <- tau_accept[pp,hh] + 1
          }
          if (irep<(0.5*burnin)){
            if ((tau_accept[pp,hh]/irep)>0.3)  tau_tuning[pp,hh] <- 1.01*tau_tuning[pp,hh]
            if ((tau_accept[pp,hh]/irep)<0.15) tau_tuning[pp,hh] <- 0.99*tau_tuning[pp,hh]
          }
        }
      }
      #----------------------------------------------------------------------------
      # Step 1/3: Sample variance-covariance-matrix
      res_hh <- Y_hh - X_hh %*% A_hh
      
      SIGMA_hh <- rinvwishart(nu = length(sl.state) + nu_prior, S = S_prior + crossprod(res_hh))
      
      #----------------------------------------------------------------------------
      # Step 1/4: Save everything
      A_draw[,,hh]        <- A_hh
      SIGMA_draw[,,hh]    <- SIGMA_hh
      res_draw[sl.state,] <- res_hh
      A_prior[,,hh]       <- Apr_hh
    }
    
    #--------------------------------------------------------------------
    # Step 2: Sample regime allocation
    gamma_prop <- gamma_draw + rnorm(1,0,sqrt(scale_g))
    S_prop     <- matrix(1, bigT, 1)
    S_prop[Z>gamma_prop,1] <- 2
    
    # current likelihood
    Lik        <- getvarpost(Y,X,Z,gamma=gamma_draw,A=A_draw,Sigma=SIGMA_draw,Ncrit)
    prior      <- dnorm(gamma_draw,g_prior,G_prior,log=TRUE)
    cond.post  <- Lik + prior
    
    # proposal likelihood
    Lik.prop    <- getvarpost(Y,X,Z,gamma=gamma_prop,A=A_draw,Sigma=SIGMA_draw,Ncrit)
    prior.prop  <- dnorm(gamma_prop,g_prior,G_prior,log=TRUE)
    cond.post.prop <- Lik.prop + prior.prop
    
    if((cond.post.prop - cond.post)>log(runif(1))){
      gamma_draw <- gamma_prop
      S_draw     <- S_prop
      Smat_draw  <- crSmat(S_draw, h)
      accept_g   <- accept_g+1
      cond.post  <- cond.post.prop
    }
    if(irep<burnin) gamma_store_burn[irep,] <- gamma_draw
    if (irep<burnin & irep>75) {
      scale_g <- s_d*var(gamma_store_burn,na.rm=TRUE) + s_d*epsilon
    }
    #--------------------------------------------------------------------
    # Step 5: sample delay paramter
    prob <- c()
    for(jj in d.min:d.max){
      if(jj==d_draw){
        prob <- c(prob,Lik+prior)
      }else{
        # if(jj>0) Z.temp <- mlag(Zraw,jj)[,jj,drop=F] else Z.temp <- Zraw
        # Z.temp <- Z.temp[(max(p,d.max)+1):(nrow(Z.temp)-fhorz),,drop=F]
        Z.temp <- Zlag[(max(plag,d.max)+1):Traw,paste0("thrsh.lag",jj),drop=FALSE]
        lik.temp <- getvarpost(Y,X,Z.temp,gamma_draw,A_draw,SIGMA_draw,Ncrit)
        prob <- c(prob,lik.temp+prior) #+log(dpois(jj,lambda=lambda_pois)) ?????
      }
    }
    prob<-exp(prob-max(prob)) # safe exponential
    prob<-prob/sum(prob)
    d_draw <- sample(d.min:d.max,size=1,prob=prob)
    Z <- Zlag[(max(plag,d.max)+1):Traw,paste0("thrsh.lag",d_draw),drop=FALSE]
    S_draw <- as.matrix(ifelse(Z>gamma_draw,2,1))
    Smat_draw <- crSmat(S_draw, h)
    #----------------------------------------------------------------------------
    # Step 4: store draws
    if(irep %in% thin.draws){
      count <- count+1
      A_store[count,,,]     <- A_draw
      SIGMA_store[count,,,] <- SIGMA_draw
      res_store[count,,]    <- res_draw
      # NG
      if(save.prior){
        Aprior_store[count,,,] <- A_prior
        lambda2_store[count,,] <- lambda2_draw
        tau_store[count,,]     <- tau_draw
      }
      # thresh stuff
      S_store[count,,] <- S_draw
      gamma_store[count,] <- gamma_draw
      d_store[count,] <- d_draw
    }
    # show some stuff
    if(irep%%50==0){
      print(paste0("Round: ",irep))
    }
  }
  dimnames(A_store)<-list(NULL,varNameslags,varNames,NULL)
  dimnames(SIGMA_store)<-list(NULL,varNames,varNames,NULL)
  dimnames(res_store)<-list(NULL,NULL,varNames)
  
  # define output
  out <- list(Y=Y, Z=Z, A=A_store, res=res_store, SIGMA=SIGMA_store,
              Aprior=Aprior_store, lambda2=lambda2_store, tau=tau_store,
              gamma=gamma_store, d=d_store, S=S_store,
              args=arglist)
  return(out)
}

#-----------------------------------------------------------------------------------------------

bvar_wishart_ng <- function(Yraw, plag = 1, args = NULL){
  #----------------------------------------INPUTS----------------------------------------------------#
  # prepare arguments
  draws = burnin = 5000; SV = FALSE; cons=TRUE; trend=FALSE; qtrend=FALSE; thin=1; Ex=NULL; save.prior=FALSE
  if(!is.null(args)){
    for(aa in c("draws","burnin","SV","cons","trend","qtrend","thin","Ex","save.prior")){
      if(aa%in%names(args)) assign(aa, args[[aa]])
    }
  }
  arglist=list(Yraw=Yraw, plag=plag, draws=draws, burnin=burnin, SV=SV, cons=cons, trend=trend, qtrend=qtrend,
               thin=thin, Ex=Ex, save.prior=save.prior)
  
  #----------------------------------------PACKAGES--------------------------------------------------#
  require(GIGrvg, quietly=TRUE)
  require(Rcpp, quietly=TRUE)
  require(MASS, quietly=TRUE)
  require(mvtnorm, quietly=TRUE)
  require(LaplacesDemon, quietly=TRUE)
  sourceCpp("./scripts/do_rgig.cpp")
  
  #-------------------------------------------START--------------------------------------------------#
  Traw  <- nrow(Yraw)
  n     <- ncol(Yraw)
  K     <- n*plag
  Ylag  <- mlag(Yraw,plag)
  varNames <- colnames(Yraw)
  if(is.null(varNames)) varNames <- paste0("Y",seq(n))
  varNameslags <- NULL
  for(pp in 1:plag) varNameslags <- c(varNameslags,paste(varNames,".lag",pp,sep=""))
  colnames(Ylag) <- varNameslags
  
  texo <- FALSE; nex <- 0; Exraw <- NULL
  if(!is.null(Ex)){
    Exraw <- Ex; nex <- ncol(Exraw)
    texo <- TRUE
    ExNames <- paste0("Ex.",colnames(Exraw))
    if(is.null(ExNames)) ExNames <- paste0("Ex.",seq(1,nex))
    varNameslags <- c(varNameslags, ExNames)
  }
  
  X <- cbind(Ylag,Exraw)
  X <- X[(plag+1):nrow(X),,drop=FALSE]
  Y <- Yraw[(plag+1):Traw,,drop=FALSE]
  bigT  <- nrow(X)
  
  if(cons){
    X <- cbind(X,1)
    varNameslags <- c(varNameslags,"cons")
    colnames(X)[ncol(X)] <- "cons"
  }
  if(trend){
    X <- cbind(X,seq(1,bigT))
    varNameslags <- c(varNameslags,"trend")
    colnames(X)[ncol(X)] <- "trend"
  }
  if(qtrend){
    X <- cbind(X,seq(1,bigT)^2)
    varNameslags <- c(varNameslags,"qtrend")
    colnames(X)[ncol(X)] <- "qtrend"
  }
  
  k = ncol(X)
  q = k*n
  v = (n*(n-1))/2
  
  #---------------------------------------------------------------------------------------------------------
  # HYPERPARAMETERS
  #---------------------------------------------------------------------------------------------------------
  # compute AR residuals of AR(4)
  arsig2 <- get_resid_ar(Yraw, 4)
  
  S_prior  <- diag(c(arsig2))
  nu_prior <- n + 2
  prmean    = 1
  # NG
  d_lambda   = 0.01
  e_lambda   = 0.01
  tau_start  = 0.7
  sample_tau = TRUE
  
  #---------------------------------------------------------------------------------------------------------
  # OLS Quantitites
  #---------------------------------------------------------------------------------------------------------
  XtXinv <- try(solve(crossprod(X)),silent=TRUE)
  if(is(XtXinv,"try-error")) XtXinv <- MASS::ginv(crossprod(X))
  A_OLS  <- XtXinv%*%(t(X)%*%Y)
  E_OLS  <- Y - X%*%A_OLS
  #a_OLS <- as.vector(A_OLS)
  #SSE  <-  t((Y - X%*%A_OLS))%*%(Y - X%*%A_OLS)
  SIGMA_OLS  <- crossprod(E_OLS)/(bigT-k)
  #IXY  <-   kronecker(diag(M),(t(X)%*%Y))
  
  #---------------------------------------------------------------------------------------------------------
  # Initial Values
  #---------------------------------------------------------------------------------------------------------
  A_draw     <- A_OLS
  SIGMA_draw <- SIGMA_OLS
  C_draw     <- t(chol(SIGMA_draw))
  
  #---------------------------------------------------------------------------------------------------------
  # PRIORS
  #---------------------------------------------------------------------------------------------------------
  # prior mean
  a_prior <- matrix(0,k,n)
  diag(a_prior) <- prmean
  
  # prior variance
  A_prior <- matrix(10,k,n)
  
  # NG stuff
  lambda2_draw <- matrix(0.01, plag, 1, dimnames=list(paste0("lag.",seq(plag)),NULL))
  tau_draw     <- matrix(tau_start, plag, 1, dimnames=dimnames(lambda2_draw))
  tau_tuning   <- matrix(.43, plag, 1, dimnames=dimnames(lambda2_draw))
  tau_accept   <- matrix(0, plag, 1, dimnames=dimnames(lambda2_draw))
  
  #---------------------------------------------------------------------------------------------------------
  # SAMPLER MISCELLANEOUS
  #---------------------------------------------------------------------------------------------------------
  ntot  <- draws+burnin
  
  # thinning
  count <- 0
  thindraws    <- draws/thin
  thin.draws   <- seq(burnin+1,ntot,by=thin)
  arglist      <- c(arglist, thindraws=thindraws)
  
  #---------------------------------------------------------------------------------------------------------
  # STORAGES
  #---------------------------------------------------------------------------------------------------------
  A_store       <- array(NA_real_, c(thindraws, k, n))
  SIGMA_store   <- array(NA_real_, c(thindraws, n, n))
  res_store     <- array(NA_real_, c(thindraws, bigT, n))
  # NG
  if(save.prior){
    Aprior_store  <- array(NA_real_, c(thindraws, k, n))
    lambda2_store <- array(NA_real_, c(thindraws, plag, 1))
    tau_store     <- array(NA_real_, c(thindraws, plag, 1))
  }else{
    Aprior_store <- lambda2_store <- tau_store <- NULL
  }
  
  #---------------------------------------------------------------------------------------------------------
  # MCMC LOOP
  #---------------------------------------------------------------------------------------------------------
  for (irep in 1:ntot){
    #----------------------------------------------------------------------------
    # Step 1: Sample coefficients
    
    SIGMAinv_draw <- solve(SIGMA_draw)
    V_post <- try(chol2inv(chol(diag(as.vector(1/A_prior)) + SIGMAinv_draw %x% crossprod(X))),silent=TRUE)
    if(is(V_post,"try-error")) V_post = try(solve(diag(as.vector(1/A_prior)) + SIGMAinv_draw %x% crossprod(X)), silent=TRUE)
    if(is(V_post,"try-error")) V_post = ginv(diag(as.vector(1/A_prior)) + SIGMAinv_draw %x% crossprod(X))
    a_post <- V_post %*% (diag(as.vector(1/A_prior)) %*% as.vector(a_prior) + SIGMAinv_draw %x% crossprod(X) %*% as.vector(A_OLS))
    
    A_draw <- matrix(a_post + t(chol(V_post)) %*% rnorm(q),k,n)
    dimnames(A_draw)<-list(varNameslags, varNames)
    
    #----------------------------------------------------------------------------
    # Step 2: different shrinkage prior setups
    #----------------------------------------------------------------------------
    # Normal-Gamma for endogenous variables
    for (pp in 1:plag){
      slct.i  <- grep(paste0("\\.lag",pp), rownames(A_draw))
      A.lag   <- A_draw[slct.i,,drop=FALSE]
      a.prior <- a_prior[slct.i,,drop=FALSE]
      A.prior <- A_prior[slct.i,,drop=FALSE]
      
      n.end = nrow(A.lag)
      if(pp==1){
        lambda2_draw[pp,1] <- rgamma(n = 1,
                                     shape = d_lambda + tau_draw[pp,1]*n*n.end,
                                     rate = e_lambda + 0.5*tau_draw[pp,1]*sum(A.prior))
      }else{
        lambda2_draw[pp,1] <- rgamma(n = 1,
                                     shape = d_lambda + tau_draw[pp,1]*n*n.end,
                                     rate = e_lambda + 0.5*tau_draw[pp,1]*prod(lambda2_draw[1:(pp-1),1])*sum(A.prior))
      }
      for(jj in 1:n.end){
        for(nn in 1:n){
          temp <- do_rgig1(lambda = tau_draw[pp,1] - 0.5,
                           chi = (A.lag[jj,nn] - a.prior[jj,nn])^2,
                           psi = tau_draw[pp,1]*prod(lambda2_draw[1:pp,1]))
          A.prior[jj,nn] <- ifelse(temp<1e-8,1e-8,ifelse(temp>1e+8,1e+8,temp))
        }
      }
      A_prior[slct.i,] <- A.prior
      if (sample_tau){
        #Sample a_tau through a simple RWMH step (on-line tuning of the MH scaling within the first 50% of the burn-in phase)
        tau_prop      <- exp(rnorm(1,0,tau_tuning[pp,1]))*tau_draw[pp,1]
        post_tau_prop <- atau_post(atau=tau_prop, thetas=as.vector(A.prior), lambda2=prod(lambda2_draw[1:pp,1]), k=length(A.prior))
        post_tau_old  <- atau_post(atau=tau_draw[pp,1], thetas=as.vector(A.prior), lambda2=prod(lambda2_draw[1:pp,1]), k=length(A.prior))
        post.diff     <- post_tau_prop-post_tau_old
        post.diff     <- ifelse(is.nan(post.diff),-Inf,post.diff)
        if (post.diff > log(runif(1,0,1))){
          tau_draw[pp,1] <- tau_prop
          tau_accept[pp,1] <- tau_accept[pp,1] + 1
        }
        if (irep<(0.5*burnin)){
          if ((tau_accept[pp,1]/irep)>0.3)  tau_tuning[pp,1] <- 1.01*tau_tuning[pp,1]
          if ((tau_accept[pp,1]/irep)<0.15) tau_tuning[pp,1] <- 0.99*tau_tuning[pp,1]
        }
      }
    }
    #----------------------------------------------------------------------------
    # Step 3: Sample variances
    #----------------------------------------------------------------------------
    res_draw <- Y - X %*% A_draw
    
    SIGMA_draw <- rinvwishart(nu = bigT + nu_prior, S = S_prior + crossprod(res_draw))
    #----------------------------------------------------------------------------
    # Step 4: store draws
    #----------------------------------------------------------------------------
    if(irep %in% thin.draws){
      count <- count+1
      
      A_store[count,,]       <- A_draw
      SIGMA_store[count,,]   <- SIGMA_draw
      res_store[count,,]     <- res_draw
      # NG
      if(save.prior){
        Aprior_store[count,,]  <- A_prior
        lambda2_store[count,,] <- lambda2_draw
        tau_store[count,,]     <- tau_draw
      }
    }
    if(irep%%50==0) print(paste0("Round: ",irep))
  }
  #---------------------------------------------------------------------------------------------------------
  # END ESTIMATION
  #---------------------------------------------------------------------------------------------------------
  dimnames(A_store)=list(NULL,varNameslags,varNames)
  dimnames(SIGMA_store)=list(NULL,varNames,varNames)
  ret <- list(Y=Y, X=X, A=A_store, SIGMA=SIGMA_store, res=res_store, 
              Aprior=Aprior_store, lambda2=lambda2_store, tau=tau_store,
              args=arglist)
  return(ret)
}

#-----------------------------------------------------------------------------------------------

getvarpost <- function(Y,X,Z,gamma,A,Sigma,Ncrit){
  bigT <- nrow(Y)
  S <- matrix(1,bigT,1)
  S[Z>gamma,1] <- 2
  if(any(tabulate(S)<Ncrit)){
    lik <- -Inf
  }else{
    wrapper <- function(tt,Y,X,Z,S,A,Mu,Sigma){
      require(mvnfast)
      state.tt <- S[tt,1]
      mean.tt  <- X[tt,] %*% A[,,state.tt]
      sigm.tt  <- Sigma[,,state.tt]
      return(mvtnorm::dmvnorm(Y[tt,], mean=mean.tt, sigma=sigm.tt, log=TRUE))
    }
    lik <- sum(sapply(1:bigT, wrapper,Y,X,Z,S,A,Mu,Sigma))
  }
  return(lik)
}

#-----------------------------------------------------------------------------------------------

crSmat <- function(S, H) {
  Smat <- matrix(0, nrow = length(S), ncol = H)
  for(h in 1:H) Smat[S==h,h] <- 1
  return(Smat)
}

#-----------------------------------------------------------------------------------------------

get_resid_ar <- function(Yraw,plag){
  # get dimensions
  n    = ncol(Yraw)
  bigT = nrow(Yraw)-plag
  
  # sample variance of AR(plag) process
  sigma_sq  <- matrix(0,n,1)
  for(nn in 1:n){
    Ylag_nn        = mlag(Yraw[,nn],plag)
    Ylag_nn        = Ylag_nn[(plag+1):nrow(Ylag_nn),,drop=FALSE]
    Y_nn           = Yraw[(plag+1):nrow(Yraw),nn,drop=FALSE]
    Ylag_nn        = cbind(Ylag_nn,1)
    alpha_nn       = solve(crossprod(Ylag_nn))%*%crossprod(Ylag_nn,Y_nn)
    sigma_sq[nn,1] = (1/bigT)*t(Y_nn-Ylag_nn%*%alpha_nn)%*%(Y_nn-Ylag_nn%*%alpha_nn)
  }
  
  return(sigma_sq)
}

#-----------------------------------------------------------------------------------------------

atau_post <- function(atau,lambda2,thetas,k,rat=1){
  logpost <- sum(dgamma(thetas,atau,(atau*lambda2/2),log=TRUE))+dexp(atau,rate=rat,log=TRUE)
  return(logpost)
}

#-----------------------------------------------------------------------------------------------

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

#-----------------------------------------------------------------------------------------------

gen_compMat <- function(A, M, p){
  Jm          <- matrix(0, M*p, M)
  Jm[1:M,1:M] <- diag(M)
  
  Cm  <- matrix(0, M*p, M*p)
  if(p==1) Cm <- t(A) else {
    for(j in 1:(p-1)){
      Cm[(j*M+1):(M*(j+1)),(M*(j-1)+1):(j*M)] <- diag(M)
    }
  }
  bbtemp <- A[1:(M*p),]
  splace <- 0
  for(ii in 1:p){
    for(iii in 1:M) {
      Cm[iii,((ii-1)*M+1):(ii*M)] <- t(bbtemp[(splace+1):(splace+M),iii])
    }
    splace <- splace+M
  }
  return(list(Cm=Cm,
              Jm=Jm))
}

#-----------------------------------------------------------------------------------------------

mlag <- function(X,lag){
  p <- lag
  X <- as.matrix(X)
  Traw <- nrow(X)
  N <- ncol(X)
  Xlag <- matrix(0,Traw,p*N)
  for (ii in 1:p){
    Xlag[(p+1):Traw,(N*(ii-1)+1):(N*ii)] <- X[(p+1-ii):(Traw-ii),(1:N)]
  }
  colnames(Xlag) <- paste0(colnames(X),".lag",rep(seq(p),each=N))
  return(Xlag)  
}

#-----------------------------------------------------------------------------------------------

pct <- function(x,p,f) {log(x/dplyr::lag(x,n=p))*100*f}

#-----------------------------------------------------------------------------------------------

.atau_post <- function(atau,lambda2,thetas,k,rat=1){
  logpost <- sum(dgamma(thetas,atau,(atau*lambda2/2),log=TRUE))+dexp(atau,rate=rat,log=TRUE)
  return(logpost)
}

#-----------------------------------------------------------------------------------------------

.construct.arglist = function (funobj, envir = NULL){
  namedlist = formals(funobj)
  argnames = names(namedlist)
  if (!is.environment(envir))
    envir = sys.frame(-1)
  for (argn in 1:length(namedlist)) {
    testval = as.logical(try(exists(argnames[argn], envir = envir),
                             silent = TRUE))
    if (is.na(testval))
      testval = FALSE
    if (testval) {
      testout = try(get(argnames[argn], envir = envir),silent = TRUE)
      if (is.null(testout)) {
        namedlist[[argn]] = "list(NULL)blabla"
      } else {
        namedlist[[argn]] = testout
      }
    }
  }
  namedlist = lapply(namedlist,function(x) if (any(x=="list(NULL)blabla")) NULL else x)
  lapply(namedlist, function(l) if(any(l=="list(NULL)blabla")){NULL}else{l})
  return(namedlist)
}

#-----------------------------------------------------------------------------------------------

transx <- function(x,tcode,lag){
  if(!is.matrix(x)) x<-as.matrix(x)
  # transform x series
  # -- tcodes:
  #       1 Level
  #       2 First Difference
  #       3 Second Difference
  #       4 Log-Level
  #       5 Log-First-Difference
  #       6 Log-Second-Difference
  #       7 Detrend Log Using 1-sided HP detrending for Monthly data
  #       8 Detrend Log Using 1-sided HP detrending for Quarterly data
  #       9 (1-L)/(1-L^12)
  #
  #       Translated from the Gauss procs of Stock&Watson(2005),'Implications of
  #       dynamic factor models for VAR analysis'
  #       Dimitris Korobilis, June 2007
  #       
  #       Translated from Matlab Code by Dimitris Korobilis to R
  #       Maximilian Böck, November 2019
  
  small   = 1e-06
  relvarm = .00000075
  relvarq = .000625
  
  n <- nrow(x)
  y <- matrix(NA,n,1)
  
  if(tcode==1){
    y<-x
  }else if(tcode==2){
    y[(lag+1):n,] <- x[(lag+1):n,] - x[1:(n-lag),]
  }else if(tcode==3){
    y[(lag*2+1):n,] <- x[(lag*2+1):n,] - 2*x[(lag+1):(n-lag),] + x[1:(n-lag*2),]
  }else if(tcode==4||tcode==5||tcode==6||tcode==7||tcode==8||tcode==9){
    if(min(x) < small){
      y[1:n,]<-NA
    }else{
      if(tcode==4) y[1:n,] <- log(x)
      if(tcode==5) y[(lag+1):n,] <- log(x)[(lag+1):n,]-log(x)[1:(n-lag),]
      if(tcode==6) y[(2*lag+1):n,] <- diff(log(x), differences=2)
      if(tcode==7) y[1:n,] <- detrend(log(x),relvarm)[["xc"]]
      if(tcode==8) y[1:n,] <- detrend(log(x),relvarq)[["xc"]]
      if(tcode==9) y[14:n,] <- log(x)[14:n,]-log(x)[13:(n-1),]-log(x)[2:(n-12)]+log(x)[1:(n-13),]
    }
  }
  return(y)
}

#-----------------------------------------------------------------------------------------------

extract <- function(data,k){
  t <- nrow(data);n <- ncol(data)
  xx <- crossprod(data)
  eigs <- eigen(xx)
  evec <- eigs$vectors; eval <- eigs$values
  
  lam <- sqrt(n)*evec[,1:k]
  fac <- data%*%lam/n
  fac <- apply(fac, 2, scale)
  
  return(list(fac,lam,eval))
}

#-----------------------------------------------------------------------------------------------

# btvar <- function(Yraw, plag, d.min, d.max, Zraw, nsave=5000, nburn=5000, thin=1,
#                   cons=FALSE, trend=FALSE, sv=FALSE, eigen=TRUE) {
#   args <- .construct.arglist(btvar)
#   #---------------------------load Packages------------------------------------------#
#   require(stochvol,quietly=TRUE)
#   require(GIGrvg, quietly=TRUE)
#   require(Rcpp, quietly = TRUE)
#   sourceCpp("./scripts/do_rgig.cpp")
#   #--------------------------------------------------------------------
#   Traw <- nrow(Yraw)
#   M    <- ncol(Yraw)
#   varNames <- colnames(Yraw)
#   if(is.null(varNames)) varNames <- paste0("Y",seq(M))
#   K     <- M*plag
#   Ylag  <- mlag(Yraw,plag)
#   nameslags <- NULL
#   for(pp in 1:plag) nameslags <- c(nameslags,paste(varNames,".lag",pp,sep=""))
#   colnames(Ylag) <- nameslags
#   Zlag <- mlag(Zraw,d.max)
#   colnames(Zlag)<-paste0("thrsh.lag",seq(d.max))
#   
#   X <- Ylag[(max(plag,d.max)+1):Traw,,drop=FALSE]
#   Y <- Yraw[(max(plag,d.max)+1):Traw,,drop=FALSE]
#   Z <- Zlag[(max(plag,d.max)+1):Traw,1,drop=FALSE]
#   bigT  <- nrow(Y)
#   
#   if(cons){
#     X <- cbind(X,1)
#     colnames(X)[ncol(X)] <- "cons"
#   }
#   if(trend){
#     X <- cbind(X,seq(1,bigT))
#     colnames(X)[ncol(X)] <- "trend"
#   }
#   
#   k     <- ncol(X)
#   n <- k*M
#   v <- (M*(M-1))/2
#   #---------------------------------------------------------------------------------------------------------
#   # HYPERPARAMETERS
#   #---------------------------------------------------------------------------------------------------------
#   prmean    <- 1
#   a_1       <- 0.01
#   b_1       <- 0.01
#   # SV
#   Bsigma    <- 1
#   a0        <- 1.5
#   b0        <- 25
#   bmu       <- 0
#   Bmu       <- 25
#   # prior == 3: NG
#   d_lambda  <- 0.01
#   e_lambda  <- 0.01
#   tau_start <- 0.7
#   sample_tau <- TRUE
#   # threshold stuff
#   h         <- 2 # number of regimes
#   #-----------------------------------------Initialize Gibbs sampler--------------#
#   # OLS quantities and initial values
#   A_OLS <- array(NA, c(k, M, h))
#   E_OLS <- array(NA, c(bigT, M))
#   SIGMA_OLS <- array(NA, c(M, M, h))
#   thrsh_mean <- mean(Z)
#   S_mean <- matrix(1, bigT, 1)
#   S_mean[Z>thrsh_mean,1] <- 2
#   for(hh in 1:h){
#     Y.s <- Y[S_mean==hh,,drop=F]
#     X.s <- X[S_mean==hh,,drop=F]
#     
#     A_OLS[,,hh]        <- solve(crossprod(X.s)) %*% crossprod(X.s,Y.s)
#     E_OLS[S_mean==hh,] <- Y.s - X.s %*% A_OLS[,,hh]
#     SIGMA_OLS[,,hh]    <- crossprod(Y.s - X.s %*% A_OLS[,,hh])/(sum(S_mean==hh)-k)
#   }
#   #---------------------------------------------------------------------------------------------------------
#   # Initial Values
#   #---------------------------------------------------------------------------------------------------------
#   A_draw  <- A_OLS
#   SIGMA   <- array(SIGMA_OLS, c(M,M,h))
#   Em_draw <- Em_str <- E_OLS
#   L_draw  <- L_drawinv <- array(0, c(M, M, h))
#   for(hh in 1:h) L_draw[,,hh] <- diag(M)
#   for(hh in 1:h) L_drawinv[,,hh] <- solve(L_draw[,,hh])
#   gamma_draw <- thrsh_mean
#   S_draw <- S_mean
#   #---------------------------------------------------------------------------------------------------------
#   # PRIORS
#   #---------------------------------------------------------------------------------------------------------
#   # Priors on VAR coefs
#   #-----------------------------
#   # prior mean
#   a_prior <- matrix(0, k, M) # state independent
#   diag(a_prior) <- prmean
#   # prior variance
#   A_prior <- array(.1, c(k, M, h)) # state dependent
#   
#   # NG stuff
#   lambda2_draw <- matrix(0.01, plag ,h, dimnames=list(paste0("lag.",seq(1,plag)),c("regime1","regime2")))
#   tau_draw     <- matrix(tau_start, plag, h, dimnames=dimnames(lambda2_draw))
#   tau_tuning   <- matrix(.43, plag, h, dimnames=dimnames(tau_draw))
#   tau_accept   <- matrix(0, plag, h, dimnames=dimnames(tau_draw))
#   #------------------------------------
#   # Priors on coefs in H matrix of VCV
#   #------------------------------------
#   # prior mean
#   l_prior <- matrix(0,M,M) # state independent
#   
#   # prior variance
#   L_prior <- diag(M); diag(L_prior) <- 0; L_prior[lower.tri(L_prior)] <- .1
#   L_prior <- array(L_prior, c(M, M, h)) # state dependent
#   
#   # NG
#   lambda2_draw <- rbind(lambda2_draw,matrix(0.01,1,h))
#   tau_draw     <- rbind(tau_draw,matrix(tau_start,1,h))
#   tau_tuning   <- rbind(tau_tuning, matrix(.43, 1, h))
#   tau_accept   <- rbind(tau_accept, matrix(0, 1, h))
#   rownames(tau_draw)[nrow(tau_draw)] <- rownames(lambda2_draw)[nrow(lambda2_draw)] <- 
#     rownames(tau_tuning)[nrow(tau_tuning)] <- rownames(tau_accept)[nrow(tau_accept)] <- "cov"
#   #------------------------------------
#   # SV quantities
#   #------------------------------------
#   Diags_draw <- matrix(-3, M, h)
#   Sv_draw <- matrix(-3, bigT, M)
#   svdraw <- list(para=c(mu=-10,phi=.9,sigma=.2),latent=rep(-3,bigT))
#   svl <- list()
#   for (mm in 1:M) svl[[varNames[mm]]] <- svdraw
#   pars_var <- matrix(c(-3,.9,.2,-3),4,M,dimnames=list(c("mu","phi","sigma","latent0"),varNames))
#   
#   hv <- svdraw$latent
#   para <- list(mu=-3,phi=.9,sigma=.2)
#   Sv_priors <- specify_priors(mu=sv_normal(mean=bmu, sd=Bmu), phi=sv_beta(a0,b0), sigma2=sv_gamma(shape=0.5,rate=1/(2*Bsigma)))
#   #---------threshold----------------#
#   d_draw <- 1
#   s_d=2.4^2
#   scale=0.7
#   g_prior                 <- mean(Zraw)
#   G_prior                 <- 10
#   scale_g                 <- scale
#   accept_g                <- 0
#   Ncrit                   <- 1/20*bigT
#   Smat_draw               <- crSmat(S_draw, h)
#   # adaptive MH
#   gamma_store_burn        <- array(NA, c(nburn, 1))
#   epsilon                 <- 0.00000001
#   eigs                    <- rep(NA,h)
#   #---------------------------------------------------------------------------------------------------------
#   # SAMPLER MISCELLANEOUS
#   #---------------------------------------------------------------------------------------------------------
#   ntot  <- nsave+nburn
#   
#   # thinning
#   count <- 0
#   thindraws    <- nsave/thin
#   thin.draws   <- seq(nburn+1,ntot,by=thin)
#   #---------------------------------------------------------------------------------------------------------
#   # STORAGES
#   #---------------------------------------------------------------------------------------------------------
#   A_store       <- array(NA_real_, c(thindraws, k, M, h))
#   L_store       <- array(NA_real_, c(thindraws, M, M, h))
#   res_store     <- array(NA_real_, c(thindraws, bigT, M))
#   # variances
#   Sv_store    <- array(NA_real_, c(thindraws, bigT, M))
#   pars_store  <- array(NA_real_, c(thindraws, 4, M))
#   Diags_store <- array(NA_real_, c(thindraws, M, h))
#   # NG
#   Aprior_store  <- array(NA_real_, c(thindraws, k, M, h))
#   Lprior_store  <- array(NA_real_, c(thindraws, M, M, h))
#   lambda2_store <- array(NA_real_, c(thindraws, plag+1, h))
#   tau_store     <- array(NA_real_, c(thindraws, plag+1, h))
#   # Threshold
#   S_store    <- array(NA_real_, c(thindraws, bigT, 1))
#   gamma_store   <- array(NA_real_, c(thindraws, 1))
#   d_store       <- array(NA_real_, c(thindraws, 1))
#   for(irep in 1:ntot) {
#     #----------------------------------------------------------------------------
#     # Step 1: Sample coefficients
#     for(hh in 1:h) {
#       sl.state <- which(Smat_draw[,hh]==1)
#       Y.s      <- Y[sl.state,,drop=FALSE]
#       X.s      <- X[sl.state,,drop=FALSE]
#       Sv.s     <- Sv_draw[sl.state,,drop=F]
#       A.s      <- A_draw[,,hh]
#       L.s      <- L_draw[,,hh]
#       L.s.inv  <- L_drawinv[,,hh]
#       # Step 1: Sample coefficients - carrier/clark/marcellino 2019 JoE 
#       # for (mm in 1:M){
#       #   if (mm==1){
#       #     Y.i <- Y.s[,mm]*exp(-0.5*Sv.s[,mm])
#       #     X.i <- X.s*exp(-0.5*Sv.s[,mm])
#       #     
#       #     V_post <- try(chol2inv(chol(crossprod(X.i)+diag(1/A_prior[,mm,hh]))),silent=TRUE)
#       #     if (is(V_post,"try-error")) V_post <- MASS::ginv(crossprod(X.i)+diag(1/A_prior[,mm,hh]))
#       #     A_post <- V_post%*%(crossprod(X.i,Y.i)+diag(1/A_prior[,mm,hh])%*%a_prior[,mm])
#       #     
#       #     A.draw.i <- try(A_post+t(chol(V_post))%*%rnorm(ncol(X.i)),silent=TRUE)
#       #     if (is(A.draw.i,"try-error")) A.draw.i <- mvtnorm::mvrnorm(1,A_post,V_post)
#       #     A_draw[,mm,hh] <- A.draw.i
#       #     Em[sl.state,mm] <- Em_str[sl.state,mm] <- Y.s[,mm] - X.s%*%A.draw.i
#       #   }else{
#       #     Y.i <- Y.s[,mm]*exp(-0.5*Sv.s[,mm])
#       #     X.i <- cbind(X.s,Em[sl.state,1:(mm-1)])*exp(-0.5*Sv.s[,mm])
#       #     
#       #     V_post <- try(chol2inv(chol((crossprod(X.i)+diag(1/c(A_prior[,mm,hh],L_prior[mm,1:(mm-1),hh]))))),silent=TRUE)
#       #     if (is(V_post,"try-error")) V_post <- MASS::ginv((crossprod(X.i)+diag(1/c(A_prior[,mm,hh],L_prior[mm,1:(mm-1),hh]))))
#       #     A_post <- V_post%*%(crossprod(X.i,Y.i)+diag(1/c(A_prior[,mm,hh],L_prior[mm,1:(mm-1),hh]))%*%c(a_prior[,mm],l_prior[mm,1:(mm-1)]))
#       #     
#       #     A.draw.i <- try(A_post+t(chol(V_post))%*%rnorm(ncol(X.i)),silent=TRUE)
#       #     if (is(A.draw.i,"try-error")) A.draw.i <- mvtnorm::mvrnorm(1,A_post,V_post)
#       #     
#       #     A_draw[,mm,hh] <- A.draw.i[1:ncol(X)]
#       #     Em[sl.state,mm] <- Y.s[,mm] - X.s%*%A.draw.i[1:ncol(X)]
#       #     Em_str[sl.state,mm] <- Y.s[,mm] - X.s%*%A.draw.i[1:ncol(X)] - Em[sl.state,1:(mm-1),drop=FALSE]%*%A.draw.i[(ncol(X)+1):ncol(X.i),drop=FALSE]
#       #     L_draw[mm,1:(mm-1),hh] <- A.draw.i[(ncol(X)+1):ncol(X.i)]
#       #   }
#       # }
#       #----------------------------------------------------------------------------
#       # Chan, Carriero, Clark and Marcellino 2021 Corrigendum to Carriero/Clark/Marcellino 2019 JoE
#       # Step 1a: Sample coefficients of A matrix
#       for(mm in 1:M){
#         A0_draw = A.s
#         A0_draw[,mm] <- 0
#         
#         ztilde <- as.vector((Y.s - X.s%*%A0_draw)%*%t(L.s.inv[mm:M,,drop=FALSE])) * exp(-0.5*as.vector(Sv.s[,mm:M,drop=FALSE]))
#         xtilde <- (L.s.inv[mm:M,mm,drop=FALSE] %x% X.s) * exp(-0.5*as.vector(Sv.s[,mm:M,drop=FALSE]))
#         
#         V_post <- try(chol2inv(chol(crossprod(xtilde)+diag(1/A_prior[,mm,hh]))),silent=TRUE)
#         if(is(V_post,"try-error")) V_post <- solve(crossprod(xtilde)+diag(1/A_prior[,mm,hh]))
#         if(is(V_post,"try-error")) V_post <- ginv(crossprod(xtilde)+diag(1/A_prior[,mm,hh]))
#         A_post <- V_post%*%(crossprod(xtilde,ztilde)+diag(1/A_prior[,mm,hh])%*%a_prior[,mm])
#         
#         A.draw.i <- try(A_post+t(chol(V_post))%*%rnorm(ncol(X)),silent=TRUE)
#         if(is(A.draw.i,"try-error")) A.draw.i <- mvrnorm(1,A_post,V_post)
#         A_draw[,mm,hh] <- A.draw.i
#         Em_draw[sl.state,mm] <- Y.s[,mm]-X.s%*%A.draw.i
#       }
#       # Step 1b: Sample coefficients in L matrix
#       if(M > 1){
#         for(mm in 2:M){
#           eps.m <- Em_draw[sl.state,mm]*exp(-0.5*Sv.s[,mm,drop=TRUE])
#           eps.x <- Em_draw[sl.state,1:(mm-1),drop=FALSE]*exp(-0.5*Sv.s[,mm,drop=TRUE])
#           
#           L_post <- try(chol2inv(chol(crossprod(eps.x)+diag(1/L_prior[mm,1:(mm-1),hh],mm-1,mm-1))),silent=TRUE)
#           if(is(L_post,"try-error")) L_post <- solve(crossprod(eps.x)+diag(1/L_prior[mm,1:(mm-1),hh],mm-1,mm-1))
#           if(is(L_post,"try-error")) L_post <- ginv(crossprod(eps.x)+diag(1/L_prior[mm,1:(mm-1),hh],mm-1,mm-1))
#           l_post <- L_post%*%(crossprod(eps.x,eps.m)+diag(1/L_prior[mm,1:(mm-1),hh],mm-1,mm-1)%*%l_prior[mm,1:(mm-1)])
#           
#           L.draw.i <- try(l_post+t(chol(L_post))%*%rnorm(length(1:(mm-1))),silent=TRUE)
#           if(is(L.draw.i,"try-error")) L.draw.i <- mvrnorm(1,l_post,L_post)
#           L_draw[mm,1:(mm-1),hh] <- L.draw.i
#         }
#       }
#       # Step 1c: Compute Em_str
#       L_drawinv[,,hh] = solve(L_draw[,,hh])
#       Em_str[sl.state,] = Y.s%*%t(L_drawinv[,,hh]) - X.s%*%A_draw[,,hh]%*%t(L_drawinv[,,hh])
#     }
#     rownames(A_draw) <- colnames(X)
#     #----------------------------------------------------------------------------------------
#     # Step 2: Normal-Gamma prior 
#     for(hh in 1:h){
#       L.prior.h <- L_prior[,,hh]
#       # Normal-Gamma for Covariances
#       lambda2_draw["cov",hh] <- rgamma(n = 1,
#                                        shape = d_lambda + tau_draw["cov",hh]*v,
#                                        rate = e_lambda + 0.5*tau_draw["cov",hh]*sum(L.prior.h[lower.tri(L.prior.h)]))
#       #Step VI: Sample the prior scaling factors for covariances from GIG
#       for(mm in 2:M){
#         for(ii in 1:(mm-1)){
#           temp <- do_rgig1(lambda = tau_draw["cov",hh] - 0.5, 
#                            chi = (L_draw[mm,ii,hh] - l_prior[mm,ii])^2, 
#                            psi = tau_draw["cov",hh]*lambda2_draw["cov",hh])
#         L_prior[mm,ii,hh] = ifelse(temp<1e-8,1e-8,ifelse(temp>1e+8,1e+8,temp))
#         }
#       }
#       if(sample_tau){
#         #Sample L_tau through a simple RWMH step
#         tau_prop      <- exp(rnorm(1,0,tau_tuning["cov",hh]))*tau_draw["cov",hh]
#         post_tau_prop <- .atau_post(atau=tau_prop, thetas=L.prior.h[lower.tri(L.prior.h)], k=v, lambda2=lambda2_draw["cov",hh])
#         post_tau_old  <- .atau_post(atau=tau_draw["cov",1], thetas=L.prior.h[lower.tri(L.prior.h)], k=v, lambda2=lambda2_draw["cov",hh])
#         post.diff     <- post_tau_prop - post_tau_old
#         post.diff     <- ifelse(is.nan(post.diff),-Inf,post.diff)
#         if (post.diff > log(runif(1,0,1))){
#           tau_draw["cov",hh]   <- tau_prop
#           tau_accept["cov",hh] <- tau_accept["cov",hh] + 1
#         }
#         if (irep<(0.5*nburn)){
#           if ((tau_accept["cov",hh]/irep)>0.3)  tau_tuning["cov",hh] <- 1.01*tau_tuning["cov",hh]
#           if ((tau_accept["cov",hh]/irep)<0.15) tau_tuning["cov",hh] <- 0.99*tau_tuning["cov",hh]
#         }
#       }
#       # Normal-Gamma for endogenous variables
#       for (pp in 1:plag){
#         slct.i  <- grep(paste0("\\.lag",pp,"$"), rownames(A_draw))
#         A.lag   <- A_draw[slct.i,,hh]
#         a.prior <- a_prior[slct.i,,drop=FALSE]
#         A.prior <- A_prior[slct.i,,hh]
#         
#         if (pp==1){
#           lambda2_draw[pp,hh] <- rgamma(n = 1,
#                                         shape = d_lambda + tau_draw[pp,hh]*M^2,
#                                         rate = e_lambda + 0.5*tau_draw[pp,hh]*sum(A.prior))
#         }else{
#           lambda2_draw[pp,hh] <- rgamma(n = 1,
#                                         shape = d_lambda + tau_draw[pp,hh]*M^2,
#                                         rate = e_lambda + 0.5*tau_draw[pp,hh]*prod(lambda2_draw[1:(pp-1),hh])*sum(A.prior))
#         }
#         for (jj in 1:M){
#           for (ii in 1:M){
#           temp <- do_rgig1(lambda = tau_draw[pp,hh] - 0.5,
#                            chi = (A.lag[jj,ii] - a.prior[jj,ii])^2,
#                            psi = tau_draw[pp,hh]*prod(lambda2_draw[1:pp,hh]))
#           A.prior[jj,ii] = ifelse(temp<1e-8,1e-8,ifelse(temp>1e+8,1e+8,temp))
#           }
#         }
#         A_prior[slct.i,,hh] <- A.prior
#         if (sample_tau){
#           #Sample a_tau through a simple RWMH step (on-line tuning of the MH scaling within the first 50% of the burn-in phase)
#           tau_prop      <- exp(rnorm(1,0,tau_tuning[pp,hh]))*tau_draw[pp,hh]
#           post_tau_prop <- .atau_post(atau=tau_prop, thetas=as.vector(A.prior), lambda2=prod(lambda2_draw[1:pp,hh]), k=length(A.prior))
#           post_tau_old  <- .atau_post(atau=tau_draw[pp,hh], thetas=as.vector(A.prior), lambda2=prod(lambda2_draw[1:pp,hh]), k=length(A.prior))
#           post.diff     <- post_tau_prop-post_tau_old
#           post.diff     <- ifelse(is.nan(post.diff),-Inf,post.diff)
#           if (post.diff > log(runif(1,0,1))){
#             tau_draw[pp,hh] <- tau_prop
#             tau_accept[pp,hh] <- tau_accept[pp,hh] + 1
#           }
#           if (irep<(0.5*nburn)){
#             if ((tau_accept[pp,hh]/irep)>0.3)  tau_tuning[pp,hh] <- 1.01*tau_tuning[pp,hh]
#             if ((tau_accept[pp,hh]/irep)<0.15) tau_tuning[pp,hh] <- 0.99*tau_tuning[pp,hh]
#           }
#         }
#       }
#     }
#     #----------------------------------------------------------------------------------------  
#     # Step 3: Sample variances
#     if(sv) {
#       for (mm in 1:M){
#         para   <- as.list(pars_var[,mm])
#         para$nu = Inf; para$rho=0; para$beta<-0
#         svdraw <- svsample_fast_cpp(y=Em_str[,mm], draws=1, burnin=0, designmatrix=matrix(NA_real_),
#                                     priorspec=Sv_priors, thinpara=1, thinlatent=1, keeptime="all",
#                                     startpara=para, startlatent=Sv_draw[,mm],
#                                     keeptau=FALSE, print_settings=list(quiet=TRUE, n_chains=1, chain=1),
#                                     correct_model_misspecification=FALSE, interweave=TRUE, myoffset=0,
#                                     fast_sv=default_fast_sv)
#         svl[[mm]]     <- svdraw
#         h_            <- exp(svdraw$latent[1,])
#         para$mu       <- svdraw$para[1,"mu"]
#         para$phi      <- svdraw$para[1,"phi"]
#         para$sigma    <- svdraw$para[1,"sigma"]
#         para$latent0  <- svdraw$latent0[1,"h_0"]
#         pars_var[,mm] <- unlist(para[c("mu","phi","sigma","latent0")])
#         Sv_draw[,mm]  <- log(h_)
#       }
#       # compute SIGMA
#       for(hh in 1:H) {
#         SIGMA[,,hh] <- L_draw[,,hh] %*%
#           diag(exp(Sv_draw[bigT,])) %*% t(L_draw[,,hh])
#       }
#     } else {
#       for(hh in 1:h){
#         sl.state <- which(Smat_draw[,hh]==1)
#         for(mm in 1:M){
#           S_1 <- a_1 + 0.5*length(sl.state)
#           S_2 <- b_1 + 0.5*crossprod(Em_str[sl.state,mm])
#           
#           sig_eta <- 1/rgamma(1,S_1,S_2)
#           Diags_draw[mm,hh] <- sig_eta
#           Sv_draw[sl.state,mm] <- log(sig_eta)
#         }
#         SIGMA[,,hh] <- L_draw[,,hh] %*%
#           diag(exp(Sv_draw[sl.state[1],])) %*% t(L_draw[,,hh])
#       }
#     }
#     #--------------------------------------------------------------------
#     # Step 4: Sample regime allocation
#     gamma_prop <- gamma_draw + rnorm(1,0,sqrt(scale_g))
#     Lik        <- getvarpost(Y,X,Z,gamma=gamma_draw,A=A_draw,Sigma=SIGMA,Ncrit)
#     Lik.new    <- getvarpost(Y,X,Z,gamma=gamma_prop,A=A_draw,Sigma=SIGMA,Ncrit)
#     prior      <- dnorm(gamma_draw,g_prior,G_prior,log=TRUE)
#     prior.new  <- dnorm(gamma_prop,g_prior,G_prior,log=TRUE)
#     cond.post.old <- Lik + prior
#     cond.post.new <- Lik.new + prior.new
#     
#     if((cond.post.new-cond.post.old)>log(runif(1))){
#       gamma_draw <- gamma_prop
#       S_draw     <- matrix(1, bigT, 1)
#       S_draw[Z>gamma_draw,1] <- 2
#       Smat_draw <- crSmat(S_draw, h)
#       accept_g <- accept_g+1
#       Lik <- Lik.new
#       prior <- prior.new
#     }
#     if(irep<nburn) gamma_store_burn[irep,] <- gamma_draw
#     if (irep<nburn & irep>75) {
#       scale_g <- s_d*var(gamma_store_burn,na.rm=TRUE) + s_d*epsilon
#     }
#     #--------------------------------------------------------------------
#     # Step 5: sample delay paramter
#     prob <- c()
#     for(jj in d.min:d.max){
#       if(jj==d_draw){
#         prob <- c(prob,Lik+prior)
#       }else{
#         # if(jj>0) Z.temp <- mlag(Zraw,jj)[,jj,drop=F] else Z.temp <- Zraw
#         # Z.temp <- Z.temp[(max(p,d.max)+1):(nrow(Z.temp)-fhorz),,drop=F]
#         Z.temp <- Zlag[(max(plag,d.max)+1):Traw,paste0("thrsh.lag",jj),drop=FALSE]
#         lik.temp <- getvarpost(Y,X,Z.temp,gamma_draw,A_draw,SIGMA,Ncrit)
#         prob <- c(prob,lik.temp+prior) #+log(dpois(jj,lambda=lambda_pois)) ?????
#       }
#     }
#     prob<-exp(prob-max(prob)) # safe exponential
#     prob<-prob/sum(prob)
#     d_draw <- sample(d.min:d.max,size=1,prob=prob)
#     Z <- Zlag[(max(plag,d.max)+1):Traw,paste0("thrsh.lag",d_draw),drop=FALSE]
#     S_draw <- as.matrix(ifelse(Z>gamma_draw,2,1))
#     Smat_draw <- crSmat(S_draw, h)
#     #----------------------------------------------------------------------------
#     # Step 4: store draws
#     if(irep %in% thin.draws){
#       count <- count+1
#       A_store[count,,,] <- A_draw
#       L_store[count,,,] <- L_draw
#       res_store[count,,] <- Em_draw
#       # variances
#       if(sv){
#         Sv_store[count,,] <- Sv_draw
#         pars_store[count,,] <- pars_var
#       }else{
#         Diags_store[count,,] <- Diags_draw
#       }
#       # NG
#       Aprior_store[count,,,]  <- A_prior
#       Lprior_store[count,,,]  <- L_prior
#       lambda2_store[count,,] <- lambda2_draw
#       tau_store[count,,]     <- tau_draw
#       # thresh stuff
#       S_store[count,,] <- S_draw
#       gamma_store[count,] <- gamma_draw
#       d_store[count,] <- d_draw
#     }
#     if(irep%%50==0) print(paste0("Round: ",irep))
#   }
#   #------------------------EX POST STUFF-------------------------------------#
#   A.eigen.1 <- sapply(1:thindraws,function(irep){
#     Cm <- gen_compMat(A_store[irep,,,1], M, plag)$Cm
#     return(max(abs(Re(eigen(Cm)$values))))
#   })
#   A.eigen.2 <- sapply(1:thindraws,function(irep){
#     Cm <- gen_compMat(A_store[irep,,,2], M, plag)$Cm
#     return(max(abs(Re(eigen(Cm)$values))))
#   })
#   trim_eigen <- which(A.eigen.1<1 & A.eigen.1<1)
#   print(paste0("Model yields ", length(trim_eigen), " (", round(length(trim_eigen)/thindraws*100,2),"%) stable draws out of ", thindraws, " total draws."))
#   # kick all other draws
#   A_store <- A_store[trim_eigen,,,]
#   L_store <- L_store[trim_eigen,,,]
#   res_store <- res_store[trim_eigen,,]
#   Sv_store <- Sv_store[trim_eigen,,]
#   Diags_store <- Diags_store[trim_eigen,,]
#   pars_store <- pars_store[trim_eigen,,]
#   Aprior_store <- Aprior_store[trim_eigen,,,]
#   Lprior_store <- Lprior_store[trim_eigen,,,]
#   lambda2_store <- lambda2_store[trim_eigen,,]
#   tau_store <- tau_store[trim_eigen,,]
#   S_store <- S_store[trim_eigen,,]
#   gamma_store <- gamma_store[trim_eigen,]
#   d_store <- d_store[trim_eigen,]
#   
#   args$thindraws <- length(trim_eigen)
#   
#   # define output
#   out <- list(store=list(A_store=A_store,L_store=L_store,res_store=res_store,Sv_store=Sv_store,pars_store=pars_store,
#                          Aprior_store=Aprior_store,Lprior_store=Lprior_store,lambda2_store=lambda2_store,tau_store=tau_store,
#                          gamma_store=gamma_store,d_store=d_store,Diags_store=Diags_store,S_store=S_store),
#               args=args)
#   return(out)
# }

# -----------------------------------------------------------------------------------------------
# 
# bvar <- function(Yraw, plag, nsave=5000, nburn=5000, thin=1, cons=FALSE, trend=FALSE, sv=TRUE, eigen=TRUE){
#   args <- .construct.arglist(bvar)
#   #---------------------------load Packages------------------------------------------#
#   require(stochvol,quietly=TRUE)
#   require(GIGrvg, quietly=TRUE)
#   require(Rcpp, quietly = TRUE)
#   sourceCpp("./scripts/do_rgig.cpp")
#   #------------------------------Setup----------------------------------------------#
#   Traw  <- nrow(Yraw)
#   M     <- ncol(Yraw)
#   varNames <- colnames(Yraw)
#   if(is.null(varNames)) varNames <- paste0("Y",seq(M))
#   K     <- M*plag
#   Ylag  <- mlag(Yraw,plag)
#   nameslags <- NULL
#   for(pp in 1:plag) nameslags <- c(nameslags,paste(varNames,".lag",pp,sep=""))
#   colnames(Ylag) <- nameslags
#   
#   X <- Ylag[(plag+1):Traw,,drop=FALSE]
#   Y <- Yraw[(plag+1):Traw,,drop=FALSE]
#   bigT  <- nrow(Y)
#   
#   if(cons){
#     X <- cbind(X,1)
#     colnames(X)[ncol(X)] <- "cons"
#   }
#   if(trend){
#     X <- cbind(X,seq(1,bigT))
#     colnames(X)[ncol(X)] <- "trend"
#   }
#   
#   k     <- ncol(X)
#   n <- k*M
#   v <- (M*(M-1))/2
#   #---------------------------------------------------------------------------------------------------------
#   # HYPERPARAMETERS
#   #---------------------------------------------------------------------------------------------------------
#   prmean    <- 1
#   a_1       <- 0.01
#   b_1       <- 0.01
#   # SV
#   Bsigma    <- 1
#   a0        <- 1.5
#   b0        <- 25
#   bmu       <- 0
#   Bmu       <- 25
#   # prior == 3: NG
#   d_lambda  <- 0.01
#   e_lambda  <- 0.01
#   tau_start <- 0.7
#   sample_tau<- TRUE
#   #---------------------------------------------------------------------------------------------------------
#   # OLS Quantitites
#   #---------------------------------------------------------------------------------------------------------
#   XtXinv <- try(solve(crossprod(X)),silent=TRUE)
#   if(is(XtXinv,"try-error")) XtXinv <- ginv(crossprod(X))
#   A_OLS  <- XtXinv%*%(t(X)%*%Y)
#   E_OLS  <- Y - X%*%A_OLS
#   #a_OLS <- as.vector(A_OLS)
#   #SSE  <-  t((Y - X%*%A_OLS))%*%(Y - X%*%A_OLS)
#   SIGMA_OLS  <- crossprod(E_OLS)/(bigT-k)
#   #IXY  <-   kronecker(diag(M),(t(X)%*%Y))
#   #---------------------------------------------------------------------------------------------------------
#   # Initial Values
#   #---------------------------------------------------------------------------------------------------------
#   A_draw  <- A_OLS
#   SIGMA   <- array(SIGMA_OLS, c(M,M,bigT))
#   Em_draw <- Em_str <- E_OLS
#   L_draw  <- diag(M)
#   L_drawinv <- solve(L_draw)
#   #---------------------------------------------------------------------------------------------------------
#   # PRIORS
#   #---------------------------------------------------------------------------------------------------------
#   # Priors on VAR coefs
#   #-----------------------------
#   # prior mean
#   a_prior <- matrix(0,k,M)
#   diag(a_prior) <- prmean
#   # prior variance
#   A_prior <- matrix(10,k,M)
#   
#   # NG stuff
#   lambda2_draw <- matrix(0.01, plag, 1, dimnames=list(paste0("lag.",seq(plag)),NULL))
#   tau_draw     <- matrix(tau_start, plag, 1, dimnames=dimnames(lambda2_draw))
#   tau_tuning   <- matrix(.43, plag, 1)
#   tau_accept   <- matrix(0, plag, 1)
#   #------------------------------------
#   # Priors on coefs in H matrix of VCV
#   #------------------------------------
#   # prior mean
#   l_prior <- matrix(0,M,M)
#   
#   # prior variance
#   L_prior <- matrix(10,M,M)
#   L_prior[upper.tri(L_prior)] <- 0; diag(L_prior) <- 0
#   
#   # NG
#   lambda2_draw <- rbind(lambda2_draw, matrix(0.01, 1, 1))
#   tau_draw     <- rbind(tau_draw, matrix(tau_start, 1, 1))
#   tau_accept   <- rbind(tau_accept, matrix(0, 1, 1))
#   tau_tuning   <- rbind(tau_tuning, matrix(.43, 1, 1))
#   rownames(tau_draw)[nrow(tau_draw)] <- rownames(lambda2_draw)[nrow(lambda2_draw)] <- 
#     rownames(tau_tuning)[nrow(tau_tuning)] <- rownames(tau_accept)[nrow(tau_accept)] <- "cov"
#   #------------------------------------
#   # SV quantities
#   #------------------------------------
#   Sv_draw <- matrix(-3,bigT,M)
#   svdraw <- list(para=c(mu=-10,phi=.9,sigma=.2),latent=rep(-3,bigT))
#   svl <- list()
#   for (jj in 1:M) svl[[jj]] <- svdraw
#   pars_var <- matrix(c(-3,.9,.2,-3),4,M,dimnames=list(c("mu","phi","sigma","latent0"),NULL))
#   
#   hv <- svdraw$latent
#   para <- list(mu=-3,phi=.9,sigma=.2)
#   Sv_priors <- specify_priors(mu=sv_normal(mean=bmu, sd=Bmu), phi=sv_beta(a0,b0), sigma2=sv_gamma(shape=0.5,rate=1/(2*Bsigma)))
#   eta <- list()
#   #---------------------------------------------------------------------------------------------------------
#   # SAMPLER MISCELLANEOUS
#   #---------------------------------------------------------------------------------------------------------
#   ntot  <- nsave+nburn
#   
#   # thinning
#   count <- 0
#   thindraws    <- nsave/thin
#   thin.draws   <- seq(nburn+1,ntot,by=thin)
#   #---------------------------------------------------------------------------------------------------------
#   # STORAGES
#   #---------------------------------------------------------------------------------------------------------
#   A_store       <- array(NA_real_, c(thindraws, k, M))
#   L_store       <- array(NA_real_, c(thindraws, M, M))
#   res_store     <- array(NA_real_, c(thindraws, bigT, M))
#   # SV
#   Sv_store      <- array(NA_real_, c(thindraws, bigT, M))
#   pars_store    <- array(NA_real_, c(thindraws, 4, M))
#   # NG
#   Aprior_store  <- array(NA_real_, c(thindraws, k, M))
#   Lprior_store  <- array(NA_real_, c(thindraws, M, M))
#   lambda2_store <- array(NA_real_, c(thindraws, plag+1, 1))
#   tau_store     <- array(NA_real_, c(thindraws, plag+1, 1))
#   for(irep in 1:ntot) {
#     #----------------------------------------------------------------------------
#     # Step 1: Sample coefficients - carrier/clark/marcellino 2019 JoE 
#     # for (mm in 1:M){
#     #   if (mm==1){
#     #     Y.i <- Y[,mm]*exp(-0.5*Sv_draw[,mm])
#     #     X.i <- X*exp(-0.5*Sv_draw[,mm])
#     #     
#     #     V_post <- try(chol2inv(chol(crossprod(X.i)+diag(1/A_prior[,mm]))),silent=TRUE)
#     #     if (is(V_post,"try-error")) V_post <- MASS::ginv(crossprod(X.i)+diag(1/A_prior[,mm]))
#     #     A_post <- V_post%*%(crossprod(X.i,Y.i)+diag(1/A_prior[,mm])%*%a_prior[,mm])
#     #     
#     #     A.draw.i <- try(A_post+t(chol(V_post))%*%rnorm(ncol(X.i)),silent=TRUE)
#     #     if (is(A.draw.i,"try-error")) A.draw.i <- mvtnorm::mvrnorm(1,A_post,V_post)
#     #     A_draw[,mm] <- A.draw.i
#     #     Em[,mm] <-  Em_str[,mm] <- Y[,mm]-X%*%A.draw.i
#     #   }else{
#     #     Y.i <- Y[,mm]*exp(-0.5*Sv_draw[,mm])
#     #     X.i <- cbind(X,Em[,1:(mm-1)])*exp(-0.5*Sv_draw[,mm])
#     #     
#     #     V_post <- try(chol2inv(chol((crossprod(X.i)+diag(1/c(A_prior[,mm],L_prior[mm,1:(mm-1)]))))),silent=TRUE)
#     #     if (is(V_post,"try-error")) V_post <- MASS::ginv((crossprod(X.i)+diag(1/c(A_prior[,mm],L_prior[mm,1:(mm-1)]))))
#     #     A_post <- V_post%*%(crossprod(X.i,Y.i)+diag(1/c(A_prior[,mm],L_prior[mm,1:(mm-1)]))%*%c(a_prior[,mm],l_prior[mm,1:(mm-1)]))
#     #     
#     #     A.draw.i <- try(A_post+t(chol(V_post))%*%rnorm(ncol(X.i)),silent=TRUE)
#     #     if (is(A.draw.i,"try-error")) A.draw.i <- mvtnorm::mvrnorm(1,A_post,V_post)
#     #     
#     #     A_draw[,mm] <- A.draw.i[1:ncol(X)]
#     #     Em[,mm] <- Y[,mm]-X%*%A.draw.i[1:ncol(X)]
#     #     Em_str[,mm] <- Y[,mm]-X%*%A.draw.i[1:ncol(X)]-Em[,1:(mm-1),drop=FALSE]%*%A.draw.i[(ncol(X)+1):ncol(X.i),drop=FALSE]
#     #     L_draw[mm,1:(mm-1)] <- A.draw.i[(ncol(X)+1):ncol(X.i)]
#     #   }
#     # }
#     # rownames(A_draw) <- colnames(X)
#     #----------------------------------------------------------------------------
#     # Chan, Carriero, Clark and Marcellino 2021 Corrigendum to Carriero/Clark/Marcellino 2019 JoE
#     # Step 1a: Sample coefficients of A matrix
#     for(mm in 1:M){
#       A0_draw = A_draw
#       A0_draw[,mm] <- 0
#       
#       ztilde <- as.vector((Y - X%*%A0_draw)%*%t(L_drawinv[mm:M,,drop=FALSE])) * exp(-0.5*as.vector(Sv_draw[,mm:M,drop=FALSE]))
#       xtilde <- (L_drawinv[mm:M,mm,drop=FALSE] %x% X) * exp(-0.5*as.vector(Sv_draw[,mm:M,drop=FALSE]))
#       
#       V_post <- try(chol2inv(chol(crossprod(xtilde)+diag(1/A_prior[,mm]))),silent=TRUE)
#       if(is(V_post,"try-error")) V_post <- solve(crossprod(xtilde)+diag(1/A_prior[,mm]))
#       if(is(V_post,"try-error")) V_post <- ginv(crossprod(xtilde)+diag(1/A_prior[,mm]))
#       A_post <- V_post%*%(crossprod(xtilde,ztilde)+diag(1/A_prior[,mm])%*%a_prior[,mm])
#       
#       A.draw.i <- try(A_post+t(chol(V_post))%*%rnorm(ncol(X)),silent=TRUE)
#       if(is(A.draw.i,"try-error")) A.draw.i <- mvrnorm(1,A_post,V_post)
#       A_draw[,mm] <- A.draw.i
#       Em_draw[,mm] <- Y[,mm]-X%*%A.draw.i
#     }
#     rownames(A_draw) <- colnames(X)
#     # Step 1b: Sample coefficients in L matrix
#     if(M > 1){
#       for(mm in 2:M){
#         eps.m <- Em_draw[,mm]*exp(-0.5*Sv_draw[,mm,drop=TRUE])
#         eps.x <- Em_draw[,1:(mm-1),drop=FALSE]*exp(-0.5*Sv_draw[,mm,drop=TRUE])
#         
#         L_post <- try(chol2inv(chol(crossprod(eps.x)+diag(1/L_prior[mm,1:(mm-1)],mm-1,mm-1))),silent=TRUE)
#         if(is(L_post,"try-error")) L_post <- solve(crossprod(eps.x)+diag(1/L_prior[mm,1:(mm-1)],mm-1,mm-1))
#         if(is(L_post,"try-error")) L_post <- ginv(crossprod(eps.x)+diag(1/L_prior[mm,1:(mm-1)],mm-1,mm-1))
#         l_post <- L_post%*%(crossprod(eps.x,eps.m)+diag(1/L_prior[mm,1:(mm-1)],mm-1,mm-1)%*%l_prior[mm,1:(mm-1)])
#         
#         L.draw.i <- try(l_post+t(chol(L_post))%*%rnorm(length(1:(mm-1))),silent=TRUE)
#         if(is(L.draw.i,"try-error")) L.draw.i <- mvrnorm(1,l_post,L_post)
#         L_draw[mm,1:(mm-1)] <- L.draw.i
#       }
#     }
#     # Step 1c: Compute Em_str
#     L_drawinv <- solve(L_draw)
#     Em_str <- Y%*%t(L_drawinv) - X%*%A_draw%*%t(L_drawinv)
#     #----------------------------------------------------------------------------------------
#     # Step 2: Normal-Gamma prior 
#     # Normal-Gamma for Covariances
#     lambda2_draw["cov",1] <- rgamma(n = 1,
#                                     shape = d_lambda + tau_draw["cov",1]*v,
#                                     rate = e_lambda + 0.5*tau_draw["cov",1]*sum(L_prior[lower.tri(L_prior)]))
#     #Step VI: Sample the prior scaling factors for covariances from GIG
#     for(mm in 2:M){
#       for(ii in 1:(mm-1)){
#         temp <- do_rgig1(lambda = tau_draw["cov",1] - 0.5, 
#                         chi = (L_draw[mm,ii] - l_prior[mm,ii])^2, 
#                         psi = tau_draw["cov",1]*lambda2_draw["cov",1])
#         L_prior[mm,ii] = ifelse(temp<1e-8,1e-8,ifelse(temp>1e+8,1e+8,temp))
#       }
#     }
#     if(sample_tau){
#       #Sample L_tau through a simple RWMH step
#       tau_prop      <- exp(rnorm(1,0,tau_tuning["cov",1]))*tau_draw["cov",1]
#       post_tau_prop <- .atau_post(atau=tau_prop, thetas=L_prior[lower.tri(L_prior)], k=v, lambda2=lambda2_draw["cov",1])
#       post_tau_old  <- .atau_post(atau=tau_draw["cov",1], thetas=L_prior[lower.tri(L_prior)], k=v, lambda2=lambda2_draw["cov",1])
#       post.diff     <- post_tau_prop - post_tau_old
#       post.diff     <- ifelse(is.nan(post.diff),-Inf,post.diff)
#       if (post.diff > log(runif(1,0,1))){
#         tau_draw["cov",1]   <- tau_prop
#         tau_accept["cov",1] <- tau_accept["cov",1] + 1
#       }
#       if (irep<(0.5*nburn)){
#         if ((tau_accept["cov",1]/irep)>0.3)  tau_tuning["cov",1] <- 1.01*tau_tuning["cov",1]
#         if ((tau_accept["cov",1]/irep)<0.15) tau_tuning["cov",1] <- 0.99*tau_tuning["cov",1]
#       }
#     }
#     # Normal-Gamma for endogenous variables
#     for (pp in 1:plag){
#       slct.i  <- grep(paste0("\\.lag",pp), rownames(A_draw))
#       A.lag   <- A_draw[slct.i,,drop=FALSE]
#       a.prior <- a_prior[slct.i,,drop=FALSE]
#       A.prior <- A_prior[slct.i,,drop=FALSE]
#       
#       if (pp==1){
#         lambda2_draw[pp,1] <- rgamma(n = 1,
#                                      shape = d_lambda + tau_draw[pp,1]*M^2,
#                                      rate = e_lambda + 0.5*tau_draw[pp,1]*sum(A.prior))
#       }else{
#         lambda2_draw[pp,1] <- rgamma(n = 1,
#                                      shape = d_lambda + tau_draw[pp,1]*M^2,
#                                      rate = e_lambda + 0.5*tau_draw[pp,1]*prod(lambda2_draw[1:(pp-1),1])*sum(A.prior))
#       }
#       for (jj in 1:M){
#         for (ii in 1:M){
#           temp <- do_rgig1(lambda = tau_draw[pp,1] - 0.5,
#                            chi = (A.lag[jj,ii] - a.prior[jj,ii])^2,
#                            psi = tau_draw[pp,1]*prod(lambda2_draw[1:pp,1]))
#           A.prior[jj,ii] <- ifelse(temp<1e-8,1e-8,ifelse(temp>1e+8,1e+8,temp))
#         }
#       }
#       A_prior[slct.i,] <- A.prior
#       if (sample_tau){
#         #Sample a_tau through a simple RWMH step (on-line tuning of the MH scaling within the first 50% of the burn-in phase)
#         tau_prop      <- exp(rnorm(1,0,tau_tuning[pp,1]))*tau_draw[pp,1]
#         post_tau_prop <- .atau_post(atau=tau_prop, thetas=as.vector(A.prior), lambda2=prod(lambda2_draw[1:pp,1]), k=length(A.prior))
#         post_tau_old  <- .atau_post(atau=tau_draw[pp,1], thetas=as.vector(A.prior), lambda2=prod(lambda2_draw[1:pp,1]), k=length(A.prior))
#         post.diff     <- post_tau_prop-post_tau_old
#         post.diff     <- ifelse(is.nan(post.diff),-Inf,post.diff)
#         if (post.diff > log(runif(1,0,1))){
#           tau_draw[pp,1] <- tau_prop
#           tau_accept[pp,1] <- tau_accept[pp,1] + 1
#         }
#         if (irep<(0.5*nburn)){
#           if ((tau_accept[pp,1]/irep)>0.3)  tau_tuning[pp,1] <- 1.01*tau_tuning[pp,1]
#           if ((tau_accept[pp,1]/irep)<0.15) tau_tuning[pp,1] <- 0.99*tau_tuning[pp,1]
#         }
#       }
#     }
#     #----------------------------------------------------------------------------
#     # Step 3: Sample variances
#     if (sv){
#       for (mm in 1:M){
#         para   <- as.list(pars_var[,mm])
#         para$nu = Inf; para$rho=0; para$beta<-0
#         svdraw <- svsample_fast_cpp(y=Em_draw[,mm], draws=1, burnin=0, designmatrix=matrix(NA_real_),
#                                     priorspec=Sv_priors, thinpara=1, thinlatent=1, keeptime="all",
#                                     startpara=para, startlatent=Sv_draw[,mm],
#                                     keeptau=FALSE, print_settings=list(quiet=TRUE, n_chains=1, chain=1),
#                                     correct_model_misspecification=FALSE, interweave=TRUE, myoffset=0,
#                                     fast_sv=default_fast_sv)
#         svl[[mm]]     <- svdraw
#         h_            <- exp(svdraw$latent[1,])
#         para$mu       <- svdraw$para[1,"mu"]
#         para$phi      <- svdraw$para[1,"phi"]
#         para$sigma    <- svdraw$para[1,"sigma"]
#         para$latent0  <- svdraw$latent0[1,"h_0"]
#         pars_var[,mm] <- unlist(para[c("mu","phi","sigma","latent0")])
#         Sv_draw[,mm]  <- log(h_)
#       }
#     }else{
#       for (mm in 1:M){
#         S_1 <- a_1+bigT/2
#         S_2 <- b_1+crossprod(Em_draw[,mm])/2
#         
#         sig_eta <- 1/rgamma(1,S_1,S_2)
#         Sv_draw[,mm] <- log(sig_eta)
#       }
#     }
#     #----------------------------------------------------------------------------
#     # Step 4: store draws
#     if(irep %in% thin.draws){
#       count <- count+1
#       A_store[count,,]       <- A_draw
#       L_store[count,,]       <- L_draw
#       res_store[count,,]     <- Em_draw
#       # SV
#       Sv_store[count,,]      <- Sv_draw
#       pars_store[count,,]    <- pars_var
#       # NG
#       Aprior_store[count,,]  <- A_prior
#       Lprior_store[count,,]  <- L_prior
#       lambda2_store[count,,] <- lambda2_draw
#       tau_store[count,,]     <- tau_draw
#     }
#     if(irep%%50==0) print(paste0("Round: ",irep))
#   }
#   #------------------------EX POST STUFF-------------------------------------#
#   A.eigen <- sapply(1:thindraws,function(irep){
#     Cm <- gen_compMat(A_store[irep,,], M, plag)$Cm
#     return(max(abs(Re(eigen(Cm)$values))))
#   })
#   trim_eigen <- which(A.eigen<1)
#   print(paste0("Model yields ", length(trim_eigen), " (", round(length(trim_eigen)/thindraws*100,2),"%) stable draws out of ", thindraws, " total draws."))
#   # kick all other draws
#   A_store <- A_store[trim_eigen,,]
#   L_store <- L_store[trim_eigen,,]
#   res_store <- res_store[trim_eigen,,]
#   Sv_store <- Sv_store[trim_eigen,,]
#   pars_store <- pars_store[trim_eigen,,]
#   Aprior_store <- Aprior_store[trim_eigen,,]
#   Lprior_store <- Lprior_store[trim_eigen,,]
#   lambda2_store <- lambda2_store[trim_eigen,,]
#   tau_store <- tau_store[trim_eigen,,]
#   
#   args$thindraws <- length(trim_eigen)
#   
#   # define output
#   out <- list(store=list(A_store=A_store,L_store=L_store,res_store=res_store,Sv_store=Sv_store,pars_store=pars_store,
#                          Aprior_store=Aprior_store,Lprior_store=Lprior_store,lambda2_store=lambda2_store,tau_store=tau_store),
#               args=args)
#   return(out)
# }