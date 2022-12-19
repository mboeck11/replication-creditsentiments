#-------------------------------------------------------------------------------#
# Linear Model - Extended Information                                           #
#                                                                               #
# The Impact of Credit Market Sentiments                                        #
#                                                                               #
# Maximilian Boeck, Vienna School of International Studies                      #
# Created: 08/07/2021                                                           #
# Last Edit: 28/11/2022                                                         #
#-------------------------------------------------------------------------------#

Yraw1 <- as.matrix(dataset_est[,vars])
Qraw1 <- as.matrix(dataset_est[,proxyvar])
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
if(do_scale) Yraw1 <- apply(Yraw1, 2, scale)
rownames(Yraw1)<-rownames(Qraw1)<-rownames(Xraw1)<-as.character(time_sample)

# extract factors
temp <- extract(Xraw1,k=10)
factors <- temp[[1]]
rownames(factors) <- time_sample
colnames(factors) <- paste0("fact.",seq(1,10))

# transformations
#Yraw1[,"INDPRO"]   <- pct(Yraw1[,"INDPRO"],p=diff,f=12)
#Yraw1[,"BUSLOANS"] <- pct(Yraw1[,"BUSLOANS"],p=diff,f=12)
#Yraw1[,"CPIAUCSL"] <- pct(Yraw1[,"CPIAUCSL"],p=diff,f=12)
# Yraw1[-c(1:diff),"FFRWXSR"]  <- diff(Yraw1[,"FFRWXSR"], lag=diff)

# original estimation
# load("../02 data/US/DiagnosticExpectationsBAAT10_forecast.rda")
# Qraw1 <- as.matrix(dataset_est$BAAT10[-1] - DE[-c(565:576),"BAAT10.sv.m"])

args = list(draws = draws, burnin = burnin, thin = thin, cons = TRUE, save.prior = TRUE)

run_varext <- bvar_wishart_ng(Yraw = cbind(Yraw1, factors[,1:q]), plag = plag, args)

#------ Identification
thindraws <- run_varext$args$thindraws
K         <- M+q
varnames  <- c(colnames(Yraw1),colnames(factors[,1:q]))

irfvarext_chol_store <- array(NA_real_, c(thindraws, K, K, nhor),
                             dimnames=list(NULL, varnames, varnames, seq(nhor)))
irfvarext_extInstr_store <- array(NA_real_, c(thindraws, K, K, nhor),
                                  dimnames=list(NULL, varnames, varnames, seq(nhor)))
for(irep in 1:thindraws){
  temp    <- gen_compMat(A=run_varext$A[irep,,], M=K, p=plag)
  compMat <- temp$Cm
  Jm      <- temp$Jm
  SIGMA   <- run_varext$SIGMA[irep,,]
  res     <- run_varext$res[irep,,]
  
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
  irfvarext_chol_store[irep,,,] <- impresp1
  
  # External Instruments
  Q <- Qraw1[(plag+1):nrow(Qraw1),,drop=TRUE]
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
  irfvarext_extInstr_store[irep,,,] <- impresp2
  
  # show
  if(irep %% 50 == 0)
    cat(paste0("Round: ", irep, ".\n"))
}
idx                      <- which(!is.na(irfvarext_chol_store[,1,1,1]))
thindraws                <- length(idx)
irfvarext_chol_store     <- irfvarext_chol_store[idx,,,]
irfvarext_extInstr_store <- irfvarext_extInstr_store[idx,,,]

irfvarext_chol <- apply(irfvarext_chol_store[,,"BAAT10",], c(2,3), quantile, c(.05,.10,.16,.50,.84,.90,0.95))
irfvarext_ext  <- apply(irfvarext_extInstr_store[,,"BAAT10",], c(2,3), quantile, c(.05,.10,.16,.50,.84,.90,.95))

#------ Convergence Diagnostics

A       <- run_varext$A
SIGMA   <- run_varext$SIGMA
Aprior  <- run_varext$Aprior
lambda2 <- run_varext$lambda2
tau     <- run_varext$tau

Ineff_A <- array(NA_real_, c(K*plag+1, K))
raftd_A <- array(NA_real_, c(K*plag+1, K))
gewek_A <- array(NA_real_, c(K*plag+1, K))
for(kk in 1:(K*plag+1)){
  for(ii in 1:K){
    temp <- as.mcmc(A[,kk,ii])
    Ineff_A[kk,ii] <- thindraws/ess(temp)
    raftd_A[kk,ii] <- raftery.diag(temp,r=0.015)$resmatrix[,"I"]
    gewek_A[kk,ii] <- geweke.diag(temp, frac1=0.1, frac2=0.5)$z
  }
}

Ineff_SIGMA <- array(NA_real_, c(M, M))
raftd_SIGMA <- array(NA_real_, c(M, M))
gewek_SIGMA <- array(NA_real_, c(M, M))
for(mm in 1:M){
  for(jj in 1:M){
    temp <- as.mcmc(SIGMA[,jj,mm])
    Ineff_SIGMA[jj,mm] <- thindraws/ess(temp)
    raftd_SIGMA[jj,mm] <- raftery.diag(temp,r=0.015)$resmatrix[,"I"]
    gewek_SIGMA[jj,mm] <- geweke.diag(temp, frac1=0.1, frac2=0.5)$z
  }
}

Ineff_Aprior <- array(NA_real_, c(K*plag, K))
raftd_Aprior <- array(NA_real_, c(K*plag, K))
gewek_Aprior <- array(NA_real_, c(K*plag, K))
for(kk in 1:(K*plag)){
  for(ii in 1:K){
    temp <- as.mcmc(Aprior[,kk,ii])
    Ineff_Aprior[kk,ii] <- thindraws/ess(temp)
    raftd_Aprior[kk,ii] <- raftery.diag(temp,r=0.015)$resmatrix[,"I"]
    gewek_Aprior[kk,ii] <- geweke.diag(temp, frac1=0.1, frac2=0.5)$z
  }
}

Ineff_lambda2 <- array(NA_real_, c(dim(lambda2)[2]))
raftd_lambda2 <- array(NA_real_, c(dim(lambda2)[2]))
gewek_lambda2 <- array(NA_real_, c(dim(lambda2)[2]))
for(kk in 1:dim(lambda2)[2]){
  temp <- as.mcmc(lambda2[,kk,1])
  Ineff_lambda2[kk] <- thindraws/ess(temp)
  raftd_lambda2[kk] <- raftery.diag(temp,r=0.015)$resmatrix[,"I"]
  gewek_lambda2[kk] <- geweke.diag(temp, frac1=0.1, frac2=0.5)$z
}

Ineff_tau <- array(NA_real_, c(dim(tau)[2]))
raftd_tau <- array(NA_real_, c(dim(tau)[2]))
gewek_tau <- array(NA_real_, c(dim(tau)[2]))
for(kk in 1:dim(tau)[2]){
  temp <- as.mcmc(tau[,kk,1])
  Ineff_tau[kk] <- thindraws/ess(temp)
  raftd_tau[kk] <- raftery.diag(temp,r=0.015)$resmatrix[,"I"]
  gewek_tau[kk] <- geweke.diag(temp, frac1=0.1, frac2=0.5)$z
}

varext_conv = list(Ineff=mean(c(Ineff_A,Ineff_SIGMA,Ineff_Aprior,Ineff_lambda2,Ineff_tau),na.rm=TRUE),
                   raftd=mean(c(raftd_A,raftd_SIGMA,raftd_Aprior,raftd_lambda2,raftd_tau),na.rm=TRUE),
                   gewek=mean(c(abs(gewek_A)>1.96,abs(gewek_SIGMA)>1.96,abs(gewek_Aprior)>1.96,abs(gewek_lambda2)>1.96,abs(gewek_tau)>1.96),na.rm=TRUE),
                   percd=thindraws/draws)

rm(Yraw1, Qraw1, fit.res, ihor, impact, impresp1, impresp2, irep, irfvarext_chol_store, irfvarext_extInstr_store, 
   Q, thindraws, b11, b11b11p, b12b12p, b21ib11, compMat, compMati, Jm, reg0, res, shock, Sig11,
   Sig12, Sig21, Sig22, SIGMA, temp, ZZp, Ineff_A, Ineff_SIGMA, Ineff_Aprior, Ineff_lambda2, Ineff_tau,
   raftd_A, raftd_SIGMA, raftd_Aprior, raftd_lambda2, raftd_tau, gewek_A, gewek_SIGMA, gewek_Aprior, gewek_lambda2, gewek_tau)

