#-------------------------------------------------------------------------------#
# Linear Model Estimation                                                       #
#                                                                               #
# The Impact of Credit Market Sentiments                                        #
#                                                                               #
# Maximilian Boeck, Vienna School of International Studies                      #
# Created: 02/05/2021                                                           #
# Last Edit: 27/11/2022                                                         #
#-------------------------------------------------------------------------------#

Yraw1 <- as.matrix(dataset_est[,vars])
Qraw1 <- as.matrix(dataset_est[,proxyvar])
M     <- ncol(Yraw1)

#------ Transformations
for(mm in 1:M){
  Yraw1[,mm] <- transx(Yraw1[,mm], tcode=tcode[mm], lag=diff)*tperc[mm]
}
Yraw1 <- Yraw1[-c(1:diff),]
Qraw1 <- Qraw1[-c(1:diff),,drop=FALSE]
if(do_scale) Yraw1 <- apply(Yraw1, 2, scale)
rownames(Yraw1)<-rownames(Qraw1)<-as.character(time_sample)

#------ Estimation
# run_var <- bvar(Yraw = Yraw1, plag = plag, nsave = draws, nburn = burnin, thin = thin, 
#                 cons = TRUE, trend = FALSE, sv = FALSE, eigen = TRUE)

args = list(draws = draws, burnin = burnin, thin = thin, cons = TRUE, trend = FALSE,
            save.prior = TRUE)

run_var <- bvar_wishart_ng(Yraw1, plag, args)

#------ Identification
thindraws <- run_var$args$thindraws
M         <- ncol(Yraw1)

irfvar_chol_store <- array(NA_real_, c(thindraws, M, M, nhor),
                           dimnames=list(NULL, colnames(Yraw1), colnames(Yraw1), seq(nhor)))
irfvar_extInstr_store <- array(NA_real_, c(thindraws, M, M, nhor),
                               dimnames=list(NULL, colnames(Yraw1), colnames(Yraw1), seq(nhor)))
for(irep in 1:thindraws){
  temp    <- gen_compMat(A=run_var$A[irep,,], M=M, p=plag)
  compMat <- temp$Cm
  Jm      <- temp$Jm
  SIGMA   <- run_var$SIGMA[irep,,]
  res     <- run_var$res[irep,,]
  
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
  irfvar_chol_store[irep,,,] <- impresp1
  
  # External Instruments
  Q <- Qraw1[(plag+1):nrow(Qraw1),,drop=TRUE]
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
  irfvar_extInstr_store[irep,,,] <- impresp2
}
idx                    <- which(!is.na(irfvar_chol_store[,1,1,1]))
thindraws              <- length(idx)
irfvar_chol_store      <- irfvar_chol_store[idx,,,]
irfvar_extInstr_store  <- irfvar_extInstr_store[idx,,,]

irfvar_chol <- apply(irfvar_chol_store[,,"BAAT10",], c(2,3), quantile, c(.05,.10,.16,.50,.84,.90,0.95))
irfvar_ext  <- apply(irfvar_extInstr_store[,,"BAAT10",], c(2,3), quantile, c(.05,.10,.16,.50,.84,.90,.95))

#------ Convergence Diagnostics

A <- run_var$A
SIGMA <- run_var$SIGMA
Aprior <- run_var$Aprior
lambda2 <- run_var$lambda2
tau <- run_var$tau

Ineff_A <- array(NA_real_, c(M*plag+1, M))
raftd_A <- array(NA_real_, c(M*plag+1, M))
gewek_A <- array(NA_real_, c(M*plag+1, M))
for(kk in 1:(M*plag+1)){
  for(mm in 1:M){
    temp <- as.mcmc(A[,kk,mm])
    #Ineff_A[kk,mm] <- draws/effectiveSize(temp)
    Ineff_A[kk,mm] <- draws/ess(temp)
    raftd_A[kk,mm] <- raftery.diag(temp,r=0.015)$resmatrix[,"I"]
    gewek_A[kk,mm] <- geweke.diag(temp, frac1=0.1, frac2=0.5)$z
  }
}

Ineff_SIGMA <- array(NA_real_, c(M, M))
raftd_SIGMA <- array(NA_real_, c(M, M))
gewek_SIGMA <- array(NA_real_, c(M, M))
for(jj in 1:M){
  for(mm in 1:M){
    temp <- as.mcmc(SIGMA[,mm,jj])
    Ineff_SIGMA[mm,jj] <- draws/ess(temp)
    raftd_SIGMA[mm,jj] <- raftery.diag(temp,r=0.015)$resmatrix[,"I"]
    gewek_SIGMA[mm,jj] <- geweke.diag(temp, frac1=0.1, frac2=0.5)$z
  }
}

Ineff_Aprior <- array(NA_real_, c(M*plag, M))
raftd_Aprior <- array(NA_real_, c(M*plag, M))
gewek_Aprior <- array(NA_real_, c(M*plag, M))
for(kk in 1:(M*plag)){
  for(mm in 1:M){
    temp <- as.mcmc(Aprior[,kk,mm])
    Ineff_Aprior[kk,mm] <- draws/ess(temp)
    raftd_Aprior[kk,mm] <- raftery.diag(temp,r=0.015)$resmatrix[,"I"]
    gewek_Aprior[kk,mm] <- geweke.diag(temp, frac1=0.1, frac2=0.5)$z
  }
}

Ineff_lambda2 <- array(NA_real_, c(dim(lambda2)[2]))
raftd_lambda2 <- array(NA_real_, c(dim(lambda2)[2]))
gewek_lambda2 <- array(NA_real_, c(dim(lambda2)[2]))
for(mm in 1:dim(lambda2)[2]){
  temp <- as.mcmc(lambda2[,mm,1])
  Ineff_lambda2[mm] <- draws/ess(temp)
  raftd_lambda2[mm] <- raftery.diag(temp,r=0.015)$resmatrix[,"I"]
  gewek_lambda2[mm] <- geweke.diag(temp, frac1=0.1, frac2=0.5)$z
}

Ineff_tau <- array(NA_real_, c(dim(tau)[2]))
raftd_tau <- array(NA_real_, c(dim(tau)[2]))
gewek_tau <- array(NA_real_, c(dim(tau)[2]))
for(mm in 1:dim(tau)[2]){
  temp <- as.mcmc(tau[,mm,1])
  Ineff_tau[mm] <- draws/effectiveSize(temp)
  raftd_tau[mm] <- raftery.diag(temp,r=0.015)$resmatrix[,"I"]
  gewek_tau[mm] <- geweke.diag(temp, frac1=0.1, frac2=0.5)$z
}

var_conv = list(Ineff=mean(c(Ineff_A,Ineff_SIGMA,Ineff_Aprior,Ineff_lambda2,Ineff_tau),na.rm=TRUE),
                raftd=mean(c(raftd_A,raftd_SIGMA,raftd_Aprior,raftd_lambda2,raftd_tau),na.rm=TRUE),
                gewek=mean(c(abs(gewek_A)>1.96,abs(gewek_SIGMA)>1.96,abs(gewek_Aprior)>1.96,abs(gewek_lambda2)>1.96,abs(gewek_tau)>1.96),na.rm=TRUE),
                percd=thindraws/draws)

rm(Yraw1, Qraw1, fit.res, ihor, impact, impresp1, impresp2, irep, irfvar_chol_store, irfvar_extInstr_store, 
   Q, thindraws, b11, b11b11p, b12b12p, b21ib11, compMat, compMati, Jm, reg0, res, shock, Sig11,
   Sig12, Sig21, Sig22, SIGMA, temp, ZZp, Ineff_A, Ineff_SIGMA, Ineff_Aprior, Ineff_lambda2, Ineff_tau,
   raftd_A, raftd_SIGMA, raftd_Aprior, raftd_lambda2, raftd_tau, gewek_A, gewek_SIGMA, gewek_Aprior, gewek_lambda2, gewek_tau)

