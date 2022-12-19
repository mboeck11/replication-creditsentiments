#-------------------------------------------------------------------------------#
# Linear Model - Internal Instruments Approach                                  #
#                                                                               #
# The Impact of Credit Market Sentiments                                        #
#                                                                               #
# Maximilian Boeck, Vienna School of International Studies                      #
# Created: 11/07/2021                                                           #
# Last Edit: 28/11/2022                                                         #
#-------------------------------------------------------------------------------#

Yraw1 <- as.matrix(dataset_est[,vars])
Qraw1 <- as.matrix(dataset_est[,proxyvar])
M     <- ncol(Yraw1)

# transformations
for(mm in 1:M){
  Yraw1[,mm] <- transx(Yraw1[,mm], tcode=tcode[mm], lag=diff)*tperc[mm]
}
Yraw1 <- Yraw1[-c(1:diff),]
Qraw1 <- Qraw1[-c(1:diff),,drop=FALSE]
if(do_scale) Yraw1 <- apply(Yraw1, 2, scale)
rownames(Yraw1)<-rownames(Qraw1)<-as.character(time_sample)

args <- list(draws = draws, burnin = burnin, thin = thin, cons = TRUE)

run_varint <- bvar_wishart_ng(Yraw = cbind(Qraw1,Yraw1), plag = plag, args = args)

#------ Identification
thindraws <- run_varint$args$thindraws
K         <- ncol(Yraw1)+ncol(Qraw1)
varNames  <- c("Instrument",colnames(Yraw1))

irfvar_int_store <- array(NA_real_, c(thindraws, K, K, nhor),
                          dimnames=list(NULL, varNames, varNames, seq(nhor)))
for(irep in 1:thindraws){
  temp    <- gen_compMat(A=run_varint$A[irep,,], M=K, p=plag)
  compMat <- temp$Cm
  Jm      <- temp$Jm
  SIGMA   <- run_varint$SIGMA[irep,,]
  
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
  irfvar_int_store[irep,,,] <- impresp1
}
idx <- which(!is.na(irfvar_int_store[,1,1,1]))

irfvar_int  <- apply(irfvar_int_store[idx,,"Instrument",], c(2,3), quantile, c(.05,.10,.16,.50,.84,.90,.95))

rm(Yraw1, Qraw1, ihor, impresp1, irep, irfvar_int_store, 
   thindraws, compMat, compMati, Jm, varNames)

