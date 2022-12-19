#-------------------------------------------------------------------------------#
# Threshold Model - Extended Information                                        #
#                                                                               #
# The Impact of Credit Market Sentiments                                        #
#                                                                               #
# Maximilian Boeck, Vienna School of International Studies                      #
# Created: 02/05/2021                                                           #
# Last Edit: 28/11/2022                                                         #
#-------------------------------------------------------------------------------#

Yraw1 <- as.matrix(dataset_est[,vars])
Qraw1 <- as.matrix(dataset_est[,proxyvar])
Zraw1 <- as.matrix(dataset_est[,thrshvar])
Xraw1 <- as.matrix(dataset_est[,2:129])
M     <- ncol(Yraw1)

# transformations
for(mm in 1:M){
  Yraw1[,mm] <- transx(Yraw1[,mm], tcode=tcode[mm], lag=diff)*tperc[mm]
}
for(kk in 1:ncol(Xraw1)){
  if(any(is.na(Xraw1[,kk]))){
    Xraw1[,kk] <- NA_real_
  }else{
    Xraw1[,kk] <- transx(Xraw1[,kk], tcode=dataset_tcode[mm], lag=diff)
  }
}
Xraw1 <- Xraw1[-c(1:diff),]
Xraw1 <- Xraw1[,!apply(Xraw1,2,function(x)any(is.na(x)))]
Yraw1 <- Yraw1[-c(1:diff),]
Qraw1 <- Qraw1[-c(1:diff),,drop=FALSE]
Zraw1 <- Zraw1[-c(1:diff),,drop=FALSE]
if(do_scale) Yraw1 <- apply(Yraw1, 2, scale)
rownames(Yraw1)<-rownames(Qraw1)<-rownames(Xraw1)<-as.character(time_sample)

# extract factors
temp <- extract(Xraw1,k=10)
factors <- temp[[1]]
rownames(factors) <- time_sample
colnames(factors) <- paste0("fact.",seq(1,10))

# transformations
# Yraw1[,"INDPRO"]   <- pct(Yraw1[,"INDPRO"],p=diff,f=12)
# Yraw1[,"BUSLOANS"] <- pct(Yraw1[,"BUSLOANS"],p=diff,f=12)
# Yraw1[,"CPIAUCSL"] <- pct(Yraw1[,"CPIAUCSL"],p=diff,f=12)
# Yraw1[-c(1:diff),"FFRWXSR"]  <- diff(Yraw1[,"FFRWXSR"], lag=diff)
# Yraw1[-c(1:diff),"FEDFUNDS"] <- diff(Yraw1[,"FEDFUNDS"], lag=diff)


# original estimation
# load("../02 data/US/DiagnosticExpectationsBAAT10_forecast.rda")
# Qraw1 <- as.matrix(dataset_est$BAAT10[-1] - DE[-c(565:576),"BAAT10.sv.m"])

args = list(draws = draws, burnin = burnin, d=c(1,4), thin = thin, cons = TRUE, save.prior=TRUE)

run_tvarext <- btvar_wishart_ng(Yraw = cbind(Yraw1, factors[,1:q]), Zraw = Zraw1, plag = plag, args)

# round(apply(run_tvar$store$A_store,c(2,3,4),median),4)
# round(apply(run_tvar$store$L_store,c(2,3,4),median),4)
# round(apply(run_tvar$store$Aprior_store,c(2,3,4),median),4)
# round(apply(run_tvar$store$Lprior_store,c(2,3,4),median),4)
# round(median(run_tvar$store$gamma_store),4)
# summary(run_tvar$store$d_store)
# round(apply(run_tvar$store$Diags_store,c(2,3),median),4)
# round(apply(run_tvar$store$Sv_store,c(2,3),median),4)

#------ Identification
thindraws <- run_tvarext$args$thindraws
K         <- M+q
varnames  <- c(colnames(Yraw1),colnames(factors[,1:q]))

irftvarext_chol_store    <- array(NA_real_, c(thindraws, K, K, nhor, h),
                               dimnames=list(NULL, varnames, varnames, seq(nhor), c("regime 1", "regime 2")))
irftvarext_ext_store     <- array(NA_real_, c(thindraws, K, K, nhor, h),
                               dimnames=list(NULL, varnames, varnames, seq(nhor), c("regime 1", "regime 2")))
for(irep in 1:thindraws){
  for(hh in 1:h){
    temp    <- gen_compMat(A=run_tvarext$A[irep,,,hh], M=K, p=plag)
    compMat <- temp$Cm
    Jm      <- temp$Jm
    SIGMA   <- run_tvarext$SIGMA[irep,,,hh]
    
    if(max(abs(Re(eigen(compMat)$values)))>1) next
    
    # Cholesky Identification
    shock <- t(chol(SIGMA))
    shock <- solve(diag(diag(shock)))%*%shock
    
    impresp1 <- array(NA_real_, c(K, K, nhor))
    impresp1[,,1] <- shock
    compMati <- compMat
    for(ihor in 2:nhor){
      impresp1[,,ihor] <- t(Jm) %*% compMati %*% Jm %*% shock
      compMati <- compMati %*% compMat
    }
    irftvarext_chol_store[irep,,,,hh] <- impresp1
    
    # External Instruments
    res <- run_tvarext$res[irep,,]
    Q <- Qraw1[(plag+1):nrow(Qraw1),,drop=FALSE]
    # sl.state <- which(run_tvar$store$Smat_store[irep,,hh]==1)
    # res.s <- res[sl.state,]
    # Q.s <- Q[sl.state,]
    reg0 <- lm(res[,1] ~ Q - 1)
    fit.res <- fitted(reg0)
    b21ib11 <- t(lm(res[,-1] ~ fit.res - 1)$coef)
    Sig11   <- matrix(SIGMA[1, 1], 1, 1)
    Sig21   <- matrix(SIGMA[2:K, 1], K-1, 1)
    Sig12   <- matrix(SIGMA[1, 2:K], 1, K-1)
    Sig22   <- matrix(SIGMA[2:K, 2:K], K-1, K-1)
    ZZp     <- b21ib11%*%Sig11%*%t(b21ib11) - Sig21%*%t(b21ib11) + b21ib11%*%t(Sig21) + Sig22
    b12b12p <- t(Sig21 - b21ib11%*%Sig11) %*% solve(ZZp) %*% (Sig21 - b21ib11%*%Sig11)
    b11b11p <- Sig11 - b12b12p
    b11     <- sqrt(b11b11p)
    impact  <- c(b11, b21ib11*c(b11))
    impact  <- impact/impact[1] # normalization
    
    # create shock
    shock <- diag(K)
    shock[,1] <- impact
    
    impresp2 <- array(NA_real_, c(K, K, nhor))
    impresp2[,,1] <- shock
    compMati <- compMat
    for(ihor in 2:nhor){
      impresp2[,,ihor] <- t(Jm) %*% compMati %*% Jm %*% shock
      compMati <- compMati %*% compMat
    }
    irftvarext_ext_store[irep,,,,hh] <- impresp2
  }
  # show
  if(irep %% 50 == 0)
    cat(paste0("Round: ", irep, ".\n"))
}
idx1                  <- which(!is.na(irftvarext_chol_store[,1,1,1,1]))
idx2                  <- which(!is.na(irftvarext_chol_store[,1,1,1,2]))
idx                   <- base::intersect(idx1,idx2)
thindraws             <- length(idx)
irftvarext_chol_store <- irftvarext_chol_store[idx,,,,]
irftvarext_ext_store  <- irftvarext_ext_store[idx,,,,]

irftvarext_chol <- apply(irftvarext_chol_store[,,"BAAT10",,], c(2,3,4), quantile, c(.05,.10,.16,.50,.84,.90,0.95))
irftvarext_ext  <- apply(irftvarext_ext_store[,,"BAAT10",,], c(2,3,4), quantile, c(.05,.10,.16,.50,.84,.90,.95))

#------ Convergence Diagnostics

A       <- run_tvarext$A
SIGMA   <- run_tvarext$SIGMA
Aprior  <- run_tvarext$Aprior
lambda2 <- run_tvarext$lambda2
tau     <- run_tvarext$tau

Ineff_A <- array(NA_real_, c(K*plag+1, K, h))
raftd_A <- array(NA_real_, c(K*plag+1, K, h))
gewek_A <- array(NA_real_, c(K*plag+1, K, h))
for(hh in 1:h){
  for(kk in 1:(K*plag+1)){
    for(ii in 1:K){
      temp <- as.mcmc(A[,kk,ii,hh])
      Ineff_A[kk,ii,hh] <- thindraws/ess(temp)
      raftd_A[kk,ii,hh] <- raftery.diag(temp,r=0.015)$resmatrix[,"I"]
      gewek_A[kk,ii,hh] <- geweke.diag(temp, frac1=0.1, frac2=0.5)$z
    }
  }
}

Ineff_SIGMA <- array(NA_real_, c(M, M, h))
raftd_SIGMA <- array(NA_real_, c(M, M, h))
gewek_SIGMA <- array(NA_real_, c(M, M, h))
for(hh in 1:h){
  for(kk in 1:M){
    for(mm in 1:M){
      temp <- as.mcmc(SIGMA[,kk,mm,hh])
      Ineff_SIGMA[kk,mm,hh] <- thindraws/ess(temp)
      raftd_SIGMA[kk,mm,hh] <- raftery.diag(temp,r=0.015)$resmatrix[,"I"]
      gewek_SIGMA[kk,mm,hh] <- geweke.diag(temp, frac1=0.1, frac2=0.5)$z
    }
  }
}

Ineff_Aprior <- array(NA_real_, c(K*plag, K, h))
raftd_Aprior <- array(NA_real_, c(K*plag, K, h))
gewek_Aprior <- array(NA_real_, c(K*plag, K, h))
for(hh in 1:h){
  for(kk in 1:(K*plag)){
    for(ii in 1:K){
      temp <- as.mcmc(Aprior[,kk,ii,hh])
      Ineff_Aprior[kk,ii,hh] <- thindraws/ess(temp)
      raftd_Aprior[kk,ii,hh] <- raftery.diag(temp,r=0.015)$resmatrix[,"I"]
      gewek_Aprior[kk,ii,hh] <- geweke.diag(temp, frac1=0.1, frac2=0.5)$z
    }
  }
}

Ineff_lambda2 <- array(NA_real_, c(dim(lambda2)[2],h))
raftd_lambda2 <- array(NA_real_, c(dim(lambda2)[2],h))
gewek_lambda2 <- array(NA_real_, c(dim(lambda2)[2],h))
for(hh in 1:h){
  for(mm in 1:dim(lambda2)[2]){
    temp <- as.mcmc(lambda2[,mm,hh])
    Ineff_lambda2[mm,hh] <- thindraws/ess(temp)
    raftd_lambda2[mm,hh] <- raftery.diag(temp,r=0.015)$resmatrix[,"I"]
    gewek_lambda2[mm,hh] <- geweke.diag(temp, frac1=0.1, frac2=0.5)$z
  }
}

Ineff_tau <- array(NA_real_, c(dim(tau)[2],h))
raftd_tau <- array(NA_real_, c(dim(tau)[2],h))
gewek_tau <- array(NA_real_, c(dim(tau)[2],h))
for(hh in 1:h){
  for(mm in 1:dim(tau)[2]){
    temp <- as.mcmc(tau[,mm,hh])
    Ineff_tau[mm,hh] <- thindraws/ess(temp)
    raftd_tau[mm,hh] <- raftery.diag(temp,r=0.015)$resmatrix[,"I"]
    gewek_tau[mm,hh] <- geweke.diag(temp, frac1=0.1, frac2=0.5)$z
  }
}

tvarext_conv_reg1 = list(Ineff=mean(c(Ineff_A[,,1],Ineff_SIGMA[,,1],Ineff_Aprior[,,1],Ineff_lambda2[,1],Ineff_tau[,1]),na.rm=TRUE),
                         raftd=mean(c(raftd_A[,,1],raftd_SIGMA[,,1],raftd_Aprior[,,1],raftd_lambda2[,1],raftd_tau[,1]),na.rm=TRUE),
                         gewek=mean(c(abs(gewek_A[,,1])>1.96,abs(gewek_SIGMA[,,1])>1.96,abs(gewek_Aprior[,,1])>1.96,abs(gewek_lambda2[,1])>1.96,abs(gewek_tau[,1])>1.96),na.rm=TRUE),
                         percd=thindraws/draws)
tvarext_conv_reg2 = list(Ineff=mean(c(Ineff_A[,,2],Ineff_SIGMA[,,2],Ineff_Aprior[,,2],Ineff_lambda2[,2],Ineff_tau[,2]),na.rm=TRUE),
                         raftd=mean(c(raftd_A[,,2],raftd_SIGMA[,,2],raftd_Aprior[,,2],raftd_lambda2[,2],raftd_tau[,2]),na.rm=TRUE),
                         gewek=mean(c(abs(gewek_A[,,2])>1.96,abs(gewek_SIGMA[,,2])>1.96,abs(gewek_Aprior[,,2])>1.96,abs(gewek_lambda2[,2])>1.96,abs(gewek_tau[,2])>1.96),na.rm=TRUE),
                         percd=thindraws/draws)

rm(Yraw1, Qraw1, fit.res, ihor, impact, impresp1, impresp2, irep, irftvarext_chol_store, irftvarext_ext_store, 
   Q, thindraws, b11, b11b11p, b12b12p, b21ib11, compMat, compMati, Jm, reg0, res, shock, Sig11,
   Sig12, Sig21, Sig22, SIGMA, temp, ZZp, mm, hh, Ineff_A, Ineff_SIGMA, Ineff_Aprior, Ineff_lambda2, Ineff_tau, raftd_A, raftd_SIGMA, 
   raftd_Aprior, raftd_lambda2, raftd_tau, gewek_A, gewek_SIGMA, gewek_Aprior, gewek_lambda2, gewek_tau)
