#-------------------------------------------------------------------------------#
# Threshold Model - Robustness                                                  #
#                                                                               #
# The Impact of Credit Market Sentiments                                        #
#                                                                               #
# Maximilian Boeck, Vienna School of International Studies                      #
# Created: 02/05/2021                                                           #
# Last Edit: 28/11/2022                                                         #
#-------------------------------------------------------------------------------#

Yraw1 <- as.matrix(dataset_est[,vars])
Zraw1 <- as.matrix(dataset_est[,thrshvar])
M     <- ncol(Yraw1)

# transformations
for(mm in 1:M){
  Yraw1[,mm] <- transx(Yraw1[,mm], tcode=tcode[mm], lag=diff)*tperc[mm]
}
Yraw1 <- Yraw1[-c(1:diff),]
Zraw1 <- Zraw1[-c(1:diff),,drop=FALSE]
if(do_scale) Yraw1 <- apply(Yraw1, 2, scale)
rownames(Yraw1)<-as.character(time_sample)

irftvarrob_chol <- lapply(order, function(oo){
  #------ Identification
  args <- list(draws = draws, burnin = burnin, d = c(1,4), thin=thin, cons=TRUE)
  
  run_tvarrob <- btvar_wishart_ng(Yraw = Yraw1[,oo], Zraw=Zraw1, plag = plag, args=args)
  
  #------ Identification
  thindraws <- run_tvarrob$args$thindraws
  K         <- M
  varnames  <- colnames(Yraw1[,oo])
  
  irftvarrob_chol_store    <- array(NA_real_, c(thindraws, K, K, nhor, h),
                                    dimnames=list(NULL, varnames, varnames, seq(nhor), c("regime 1", "regime 2")))
  for(irep in 1:thindraws){
    for(hh in 1:h){
      temp    <- gen_compMat(A=run_tvarrob$A[irep,,,hh], M=K, p=plag)
      compMat <- temp$Cm
      Jm      <- temp$Jm
      SIGMA   <- run_tvarrob$SIGMA[irep,,,hh]
      
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
      irftvarrob_chol_store[irep,,,,hh] <- impresp1
    }
  }
  idx1 <- which(!is.na(irftvarrob_chol_store[,1,1,1,1]))
  idx2 <- which(!is.na(irftvarrob_chol_store[,1,1,1,2]))
  idx  <- base::intersect(idx1,idx2)
    
    return(apply(irftvarrob_chol_store[idx,,"BAAT10",,], c(2,3,4), quantile, c(.05,.10,.16,.50,.84,.90,0.95)))
})

