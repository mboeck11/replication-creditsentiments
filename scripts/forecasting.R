###############################################
### Forecast Heuristics                     ###
### The Impact of Credit Market Sentiments  ###
### Maximilian Boeck                        ###
### 01/05/2021                              ###
###############################################

# define start / end of training sample
train_start_m <- as.Date("1959-01-01", format = "%Y-%m-%d")
train_end_m   <- as.Date("1967-12-01", format = "%Y-%m-%d")
train_period  <- seq.Date(train_start_m,train_end_m, by = "1 month")

forecasts = matrix(NA_real_, Traw, 7,
                  dimnames=list(as.character(time_sample), c("RE","DE","ADA","WTR","STR","LAA","AA")))
coef      = matrix(NA_real_, Traw, 3, dimnames=list(as.character(time_sample), c("ar1","var","ind")))

for(tt in 1:length(time_sample)){
  trainT <- length(train_period)
  ytrain <- dataset_train$BAAT10[dataset_train$sasdate %in% train_period]
  
  ar_sv <- ar_svf(as.matrix(ytrain), p=1, nburn=burnin, nsave=draws)
  ar_sv_m <- ar_sv$ar
  ar_sv_v <- c(NA,ar_sv$sv)
  coef[tt, "ar1"] <- ar_sv_m
  coef[tt, "var"] <- ar_sv_v[trainT]
  coef[tt, "ind"] <- ifelse(ar_sv_m>=1,1,0)
  
  ####### w SV
  
  mu_0  <- ar_sv_m * ytrain[trainT]
  mu_1  <- ar_sv_m^2 * ytrain[trainT-1]
  sig_0 <- ar_sv_v[trainT] 
  sig_1 <- ar_sv_v[trainT-1]
  
  mu_theta  <- mu_0 + (theta*sig_0)/(sig_1 + theta*(sig_1-sig_0))*(mu_0-mu_1)
  sig_theta <- (sig_0*sig_1)/(sig_1 + theta*(sig_1-sig_0))
  
  forecasts[tt, "RE"] = mu_0
  forecasts[tt, "DE"] = mu_theta
  
  ######## Heuristics
  if(tt == 1) {
    temp <- 0
    for(ttt in 2:trainT) temp <- 0.65 * ytrain[ttt] + 0.35 * temp
    forecasts[tt, "ADA"] <- temp
  } else {
    forecasts[tt, "ADA"] <- 0.65 * ytrain[trainT] + 0.35 * forecasts[tt-1, "ADA"]
  }
  forecasts[tt, "WTR"]   <- ytrain[trainT] + 0.4 * (ytrain[trainT] - ytrain[trainT - 1])
  forecasts[tt, "STR"]   <- ytrain[trainT] + 1.3 * (ytrain[trainT] - ytrain[trainT - 1])
  forecasts[tt, "LAA"]   <- 0.5 *(mean(ytrain[(trainT-11):trainT]) + ytrain[trainT]) + (ytrain[trainT] - ytrain[trainT - 1])
  forecasts[tt, "AA"]    <- 0.5 *(mean(ytrain) + ytrain[trainT]) + (ytrain[trainT] - ytrain[trainT - 1])
  
  train_period  <- seq.Date(train_start_m,time_sample[tt], by = "1 month")
  
  print(tt)
}

# add forecast errors to dataset
dataset_est$FE.BAAT10.RE  <- c(rep(NA_real_,diff), dataset_est$BAAT10[(diff+1):nrow(dataset_est)] - forecasts[,"RE"])
dataset_est$FE.BAAT10.DE  <- c(rep(NA_real_,diff), dataset_est$BAAT10[(diff+1):nrow(dataset_est)] - forecasts[,"DE"])
dataset_est$FE.BAAT10.ADA <- c(rep(NA_real_,diff), dataset_est$BAAT10[(diff+1):nrow(dataset_est)] - forecasts[,"ADA"])
dataset_est$FE.BAAT10.WTR <- c(rep(NA_real_,diff), dataset_est$BAAT10[(diff+1):nrow(dataset_est)] - forecasts[,"WTR"])
dataset_est$FE.BAAT10.STR <- c(rep(NA_real_,diff), dataset_est$BAAT10[(diff+1):nrow(dataset_est)] - forecasts[,"STR"])
dataset_est$FE.BAAT10.LAA <- c(rep(NA_real_,diff), dataset_est$BAAT10[(diff+1):nrow(dataset_est)] - forecasts[,"LAA"])
dataset_est$FE.BAAT10.AA  <- c(rep(NA_real_,diff), dataset_est$BAAT10[(diff+1):nrow(dataset_est)] - forecasts[,"AA"])

rm(list=c("ar_sv","ar_sv_m","ar_sv_v","mu_0","mu_1","mu_theta","sig_0","sig_1","sig_theta",
          "temp","train_end_m","train_period","train_start_m","trainT","tt","ttt","ytrain"))
