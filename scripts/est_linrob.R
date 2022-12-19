#-------------------------------------------------------------------------------#
# Linear Model - Robustness                                                     #
#                                                                               #
# The Impact of Credit Market Sentiments                                        #
#                                                                               #
# Maximilian Boeck, Vienna School of International Studies                      #
# Created: 11/07/2021                                                           #
# Last Edit: 28/11/2022                                                         #
#-------------------------------------------------------------------------------#

Yraw1 <- as.matrix(dataset_est[,vars])
M     <- ncol(Yraw1)

# transformations
for(mm in 1:M){
  Yraw1[,mm] <- transx(Yraw1[,mm], tcode=tcode[mm], lag=diff)*tperc[mm]
}
Yraw1 <- Yraw1[-c(1:diff),]
if(do_scale) Yraw1 <- apply(Yraw1, 2, scale)
rownames(Yraw1)<-as.character(time_sample)

irfvarrob_chol <- lapply(order, function(oo){
  #------ Model Estimation
  args <- list(draws = draws, burnin = burnin, thin = thin, cons = TRUE, save.prior=TRUE)
  
  run_varrob <- bvar_wishart_ng(Yraw = Yraw1[,oo], plag = plag, args)
  
  #------ Identification
  thindraws <- run_varrob$args$thindraws
  K         <- M
  varnames  <- colnames(Yraw1[,oo])
  
  irfvarrob_chol_store <- array(NA_real_, c(thindraws, K, K, nhor),
                                  dimnames=list(NULL, varnames, varnames, seq(nhor)))
  for(irep in 1:thindraws){
    temp    <- gen_compMat(A=run_varrob$A[irep,,], M=K, p=plag)
    compMat <- temp$Cm
    Jm      <- temp$Jm
    SIGMA   <- run_varrob$SIGMA[irep,,]
    res     <- run_varrob$res[irep,,]
    
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
    irfvarrob_chol_store[irep,,,] <- impresp1
  }
  idx <- which(!is.na(irfvarrob_chol_store[,1,1,1]))
  
  return(apply(irfvarrob_chol_store[idx,,"BAAT10",], c(2,3), quantile, c(.05,.10,.16,.50,.84,.90,0.95)))
})

