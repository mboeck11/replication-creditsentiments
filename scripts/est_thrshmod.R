#-------------------------------------------------------------------------------#
# Threshold Model Estimation                                                    #
#                                                                               #
# The Impact of Credit Market Sentiments                                        #
#                                                                               #
# Maximilian Boeck, Vienna School of International Studies                      #
# Created: 02/05/2021                                                           #
# Last Edit: 27/11/2022                                                         #
#-------------------------------------------------------------------------------#

Yraw1 <- as.matrix(dataset_est[,vars])
Qraw1 <- as.matrix(dataset_est[,proxyvar])
Zraw1 <- as.matrix(dataset_est[,thrshvar])
M     <- ncol(Yraw1)

#------ Transformations
for(mm in 1:M){
  Yraw1[,mm] <- transx(Yraw1[,mm], tcode=tcode[mm], lag=diff)*tperc[mm]
}
Yraw1 <- Yraw1[-c(1:diff),]
Qraw1 <- Qraw1[-c(1:diff),,drop=FALSE]
Zraw1 <- Zraw1[-c(1:diff),,drop=FALSE]
if(do_scale) Yraw1 <- apply(Yraw1, 2, scale)
rownames(Yraw1)<-rownames(Qraw1)<-as.character(time_sample)

#------ Estimation

args = list(draws=draws, burnin=burnin, d=c(1,4), thin=thin, cons=TRUE, save.prior=TRUE)

run_tvar <- btvar_wishart_ng(Yraw1, Zraw1, plag, args)

#------ Identification
thindraws <- run_tvar$args$thindraws

irftvar_chol_store    <- array(NA_real_, c(thindraws, M, M, nhor, h),
                               dimnames=list(NULL, colnames(Yraw1), colnames(Yraw1), seq(nhor), c("regime 1", "regime 2")))
irftvar_ext_store     <- array(NA_real_, c(thindraws, M, M, nhor, h),
                               dimnames=list(NULL, colnames(Yraw1), colnames(Yraw1), seq(nhor), c("regime 1", "regime 2")))
for(irep in 1:thindraws){
  for(hh in 1:h){
    temp    <- gen_compMat(A=run_tvar$A[irep,,,hh], M=M, p=plag)
    compMat <- temp$Cm
    Jm      <- temp$Jm
    SIGMA   <- run_tvar$SIGMA[irep,,,hh]
    
    # check stability
    if(max(Re(eigen(compMat)$values)) > 1)
      next
    
    # Cholesky Identification
    shock <- t(chol(SIGMA))
    shock <- solve(diag(diag(shock)))%*%shock
    
    impresp1 <- array(NA_real_, c(M, M, nhor))
    impresp1[,,1] <- shock
    compMati <- compMat
    for(ihor in 2:nhor){
      impresp1[,,ihor] <- t(Jm) %*% compMati %*% Jm %*% shock
      compMati <- compMati %*% compMat
    }
    irftvar_chol_store[irep,,,,hh] <- impresp1
    
    # External Instruments
    res <- run_tvar$res[irep,,]
    Q <- Qraw1[(plag+1):nrow(Qraw1),,drop=FALSE]
    # sl.state <- which(run_tvar$store$Smat_store[irep,,hh]==1)
    # res.s <- res[sl.state,]
    # Q.s <- Q[sl.state,]
    reg0 <- lm(res[,1] ~ Q - 1)
    fit.res <- fitted(reg0)
    b21ib11 <- t(lm(res[,-1] ~ fit.res - 1)$coef)
    Sig11   <- matrix(SIGMA[1, 1], 1, 1)
    Sig21   <- matrix(SIGMA[2:M, 1], M-1, 1)
    Sig12   <- matrix(SIGMA[1, 2:M], 1, M-1)
    Sig22   <- matrix(SIGMA[2:M, 2:M], M-1, M-1)
    ZZp     <- b21ib11%*%Sig11%*%t(b21ib11) - Sig21%*%t(b21ib11) + b21ib11%*%t(Sig21) + Sig22
    b12b12p <- t(Sig21 - b21ib11%*%Sig11) %*% solve(ZZp) %*% (Sig21 - b21ib11%*%Sig11)
    b11b11p <- Sig11 - b12b12p
    b11     <- sqrt(b11b11p)
    impact  <- c(b11, b21ib11*c(b11))
    impact  <- impact/impact[1] # normalization
    
    # create shock
    shock <- diag(M)
    shock[,1] <- impact
    
    impresp2 <- array(NA_real_, c(M, M, nhor))
    impresp2[,,1] <- shock
    compMati <- compMat
    for(ihor in 2:nhor){
      impresp2[,,ihor] <- t(Jm) %*% compMati %*% Jm %*% shock
      compMati <- compMati %*% compMat
    }
    irftvar_ext_store[irep,,,,hh] <- impresp2
  }
}
idx1               <- which(!is.na(irftvar_chol_store[,1,1,1,1]))
idx2               <- which(!is.na(irftvar_chol_store[,1,1,1,2]))
idx                <- base::intersect(idx1,idx2)
thindraws          <- length(idx)
irftvar_chol_store <- irftvar_chol_store[idx,,,,]
irftvar_ext_store  <- irftvar_ext_store[idx,,,,]

irftvar_chol <- apply(irftvar_chol_store[,,"BAAT10",,], c(2,3,4), quantile, c(.05,.10,.16,.50,.84,.90,0.95))
irftvar_ext  <- apply(irftvar_ext_store[,,"BAAT10",,], c(2,3,4), quantile, c(.05,.10,.16,.50,.84,.90,.95))

#------ Robustness Identification

irftvar_robust <- array(NA_real_, c(thindraws, M, h, r))

Qrawl <- as.matrix(dataset_est[,proxyrob])
Qrawl <- Qrawl[-c(1:diff),,drop=FALSE]
rownames(Qrawl) <- time_sample

for(irep in 1:thindraws){
  for(rr in 1:r){
    for(hh in 1:h){
      SIGMA   <- run_tvar$SIGMA[irep,,,hh]
      
      # External Instruments
      res <- run_tvar$res[irep,,]
      Q <- Qrawl[(plag+1):nrow(Qraw1),rr,drop=FALSE]
      # sl.state <- which(run_tvar$store$Smat_store[irep,,hh]==1)
      # res.s <- res[sl.state,]
      # Q.s <- Q[sl.state,]
      reg0 <- lm(res[,1] ~ Q - 1)
      fit.res <- fitted(reg0)
      b21ib11 <- t(lm(res[,-1] ~ fit.res - 1)$coef)
      Sig11   <- matrix(SIGMA[1, 1], 1, 1)
      Sig21   <- matrix(SIGMA[2:M, 1], M-1, 1)
      Sig12   <- matrix(SIGMA[1, 2:M], 1, M-1)
      Sig22   <- matrix(SIGMA[2:M, 2:M], M-1, M-1)
      ZZp     <- b21ib11%*%Sig11%*%t(b21ib11) - Sig21%*%t(b21ib11) + b21ib11%*%t(Sig21) + Sig22
      b12b12p <- t(Sig21 - b21ib11%*%Sig11) %*% solve(ZZp) %*% (Sig21 - b21ib11%*%Sig11)
      b11b11p <- Sig11 - b12b12p
      b11     <- sqrt(b11b11p)
      impact  <- c(b11, b21ib11*c(b11))
      impact  <- impact/impact[1] # normalization
      
      irftvar_robust[irep,,hh,rr] <- impact
    }
  }
}

#------ Convergence Diagnostics

A       <- run_tvar$A
SIGMA   <- run_tvar$SIGMA
Aprior  <- run_tvar$Aprior
lambda2 <- run_tvar$lambda2
tau     <- run_tvar$tau

Ineff_A <- array(NA_real_, c(M*plag+1, M, h))
raftd_A <- array(NA_real_, c(M*plag+1, M, h))
gewek_A <- array(NA_real_, c(M*plag+1, M, h))
for(hh in 1:h){
  for(kk in 1:(M*plag+1)){
    for(mm in 1:M){
      temp <- as.mcmc(A[,kk,mm,hh])
      Ineff_A[kk,mm,hh] <- draws/ess(temp)
      raftd_A[kk,mm,hh] <- raftery.diag(temp,r=0.015)$resmatrix[,"I"]
      gewek_A[kk,mm,hh] <- geweke.diag(temp, frac1=0.1, frac2=0.5)$z
    }
  }
}

Ineff_SIGMA <- array(NA_real_, c(M, M, h))
raftd_SIGMA <- array(NA_real_, c(M, M, h))
gewek_SIGMA <- array(NA_real_, c(M, M, h))
for(hh in 1:h){
  for(mm in 1:M){
    for(jj in 1:M){
      temp <- as.mcmc(SIGMA[,jj,mm,hh])
      Ineff_SIGMA[jj,mm,hh] <- draws/ess(temp)
      raftd_SIGMA[jj,mm,hh] <- raftery.diag(temp,r=0.015)$resmatrix[,"I"]
      gewek_SIGMA[jj,mm,hh] <- geweke.diag(temp, frac1=0.1, frac2=0.5)$z
    }
  }
}

Ineff_Aprior <- array(NA_real_, c(M*plag, M, h))
raftd_Aprior <- array(NA_real_, c(M*plag, M, h))
gewek_Aprior <- array(NA_real_, c(M*plag, M, h))
for(hh in 1:h){
  for(kk in 1:(M*plag)){
    for(mm in 1:M){
      temp <- as.mcmc(Aprior[,kk,mm,hh])
      Ineff_Aprior[kk,mm,hh] <- thindraws/ess(temp)
      raftd_Aprior[kk,mm,hh] <- raftery.diag(temp,r=0.015)$resmatrix[,"I"]
      gewek_Aprior[kk,mm,hh] <- geweke.diag(temp, frac1=0.1, frac2=0.5)$z
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

tvar_conv_reg1 = list(Ineff=mean(c(Ineff_A[,,1],Ineff_SIGMA[,,1],Ineff_Aprior[,,1],Ineff_lambda2[,1],Ineff_tau[,1]),na.rm=TRUE),
                      raftd=mean(c(raftd_A[,,1],raftd_SIGMA[,,1],raftd_Aprior[,,1],raftd_lambda2[,1],raftd_tau[,1]),na.rm=TRUE),
                      gewek=mean(c(abs(gewek_A[,,1])>1.96,abs(gewek_SIGMA[,,1])>1.96,abs(gewek_Aprior[,,1])>1.96,abs(gewek_lambda2[,1])>1.96,abs(gewek_tau[,1])>1.96),na.rm=TRUE),
                      percd=thindraws/draws)
tvar_conv_reg2 = list(Ineff=mean(c(Ineff_A[,,2],Ineff_SIGMA[,,2],Ineff_Aprior[,,2],Ineff_lambda2[,2],Ineff_tau[,2]),na.rm=TRUE),
                      raftd=mean(c(raftd_A[,,2],raftd_SIGMA[,,2],raftd_Aprior[,,2],raftd_lambda2[,2],raftd_tau[,2]),na.rm=TRUE),
                      gewek=mean(c(abs(gewek_A[,,2])>1.96,abs(gewek_SIGMA[,,2])>1.96,abs(gewek_Aprior[,,2])>1.96,abs(gewek_lambda2[,2])>1.96,abs(gewek_tau[,2])>1.96),na.rm=TRUE),
                      percd=thindraws/draws)

rm(Yraw1, Qraw1, Qrawl, fit.res, ihor, impact, impresp1, impresp2, irep, irftvar_chol_store, irftvar_ext_store, 
   Q, thindraws, b11, b11b11p, b12b12p, b21ib11, compMat, compMati, Jm, reg0, res, shock, Sig11,
   Sig12, Sig21, Sig22, SIGMA, temp, ZZp, rr, mm, hh, Ineff_A, Ineff_SIGMA, Ineff_Aprior, Ineff_lambda2, 
   Ineff_tau, raftd_A, raftd_SIGMA, raftd_Aprior, raftd_lambda2, raftd_tau, gewek_A, gewek_SIGMA, gewek_Aprior, gewek_lambda2, gewek_tau)
