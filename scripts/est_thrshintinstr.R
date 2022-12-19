#-------------------------------------------------------------------------------#
# Threshold Model - Internal Instruments Approach                               #
#                                                                               #
# The Impact of Credit Market Sentiments                                        #
#                                                                               #
# Maximilian Boeck, Vienna School of International Studies                      #
# Created: 11/07/2021                                                           #
# Last Edit: 28/11/2022                                                         #
#-------------------------------------------------------------------------------#

Yraw1 <- as.matrix(dataset_est[,vars])
Qraw1 <- as.matrix(dataset_est[,proxyvar])
Zraw1 <- as.matrix(dataset_est[,thrshvar])
M     <- ncol(Yraw1)

#------ Transformation
for(mm in 1:M){
  Yraw1[,mm] <- transx(Yraw1[,mm], tcode=tcode[mm], lag=diff)*tperc[mm]
}
Yraw1 <- Yraw1[-c(1:diff),]
Qraw1 <- Qraw1[-c(1:diff),,drop=FALSE]
Zraw1 <- Zraw1[-c(1:diff),,drop=FALSE]
if(do_scale) Yraw1 <- apply(Yraw1, 2, scale)
rownames(Yraw1)<-rownames(Qraw1)<-as.character(time_sample)

#------ Estimation

args <- list(draws = draws, burnin = burnin, d = c(1,4), thin = thin, cons = TRUE)

run_tvarint <- btvar_wishart_ng(Yraw = cbind(Qraw1,Yraw1), Zraw = Zraw1, plag = plag, args = args)

#------ Identification
thindraws <- run_tvarint$args$thindraws
varNames  <- c("Instrument",colnames(Yraw1))
K         <- ncol(Yraw1)+ncol(Qraw1)

irftvar_int_store     <- array(NA_real_, c(thindraws, K, K, nhor, h),
                               dimnames=list(NULL, varNames, varNames, seq(nhor), c("regime 1", "regime 2")))
for(irep in 1:thindraws){
  for(hh in 1:h){
    temp    <- gen_compMat(A=run_tvarint$A[irep,,,hh], M=K, p=plag)
    compMat <- temp$Cm
    Jm      <- temp$Jm
    SIGMA   <- run_tvarint$SIGMA[irep,,,hh]
    
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
    irftvar_int_store[irep,,,,hh] <- impresp1
  }
}
idx1 <- which(!is.na(irftvar_int_store[,1,1,1,1]))
idx2 <- which(!is.na(irftvar_int_store[,1,1,1,2]))
idx  <- base::intersect(idx1,idx2)

irftvar_int <- apply(irftvar_int_store[idx,,"Instrument",,], c(2,3,4), quantile, c(.05,.10,.16,.50,.84,.90,0.95))


rm(Yraw1, Qraw1, ihor, impresp1, irep, irftvar_int_store, 
   thindraws, compMat, compMati, Jm, shock)
