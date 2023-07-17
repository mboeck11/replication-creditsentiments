#-------------------------------------------------------------------------------#
# Main Replication Files: Creates Figures and Tables for the Paper              #
#                                                                               #
# The Impact of Credit Market Sentiments                                        #
#                                                                               #
# Maximilian Boeck, Vienna School of International Studies                      #
# Created:   02/05/2021                                                         #
# Last Edit: 17/07/2023                                                         #
#-------------------------------------------------------------------------------#
rm(list=ls())
set.seed(571)

# set working directory - delete this later
setwd("/users/mboeck/dropbox/projects/!credit-sentiments/replication-creditsentiments")

# 1) Load Stuff / Settings
source("./scripts/aux.R")
library(coda)
library(mcmcse)

# draws and burn-in's used for all involved samplers
draws  = 10000 # full: 10.000
burnin = 15000 # full: 15.000

# define sample
begin_sample <- as.Date("1968-01-01",format = "%Y-%m-%d")
end_sample   <- as.Date("2014-12-01",format = "%Y-%m-%d")
time_sample  <- seq.Date(begin_sample, end_sample, by = "1 month")
Traw         <- length(time_sample)

# variables in VAR
vars     <- c("BAAT10","INDPRO","BUSLOANS","CPIAUCSL","FFRWXSR")
tcode    <- c(1,5,5,5,1)
tperc    <- c(1,100,100,100,1)
proxyvar <- "FE.BAAT10.DE"
proxyrob <- c("FE.BAAT10.DE","FE.BAAT10.RE","FE.BAAT10.ADA","FE.BAAT10.WTR","FE.BAAT10.STR","FE.BAAT10.LAA")
thrshvar <- "BAAT10"
M <- length(vars) # number of variables in the VAR
r <- length(proxyrob)

# parameters
theta <- 0.91   # diagnosticity parameter
diff  <- 12     # difference parameter: 1 = month-on-month growth rate, 12 = year-on-year growth rate
plag  <- 13      # number of lags in the VAR
thin  <- 1      # thinning factor
nhor  <- 61     # impulse response horizon
h     <- 2      # number of regimes
q_est <- c(3,7) # number of factors
order <- list(rob1=c(2,3,4,1,5),rob2=c(2,3,4,5,1))
do_scale <- FALSE

# graphical settings
varnames_plot <- c("Credit Spread", "Economic Activity", "Credit Volume", "Prices", "Interest Rate")
width <- 3000
height <- 1400
height_higher <- 2000

# 2) Build Dataset
source("./scripts/data_prep.R")

# 3) Forecasting Credit Spreads
if(file.exists(paste0("./saves/forecasting_diff=",diff,"_plag=",plag,"_draws=",draws+burnin,".rda"))){
  load(paste0("./saves/forecasting_diff=",diff,"_plag=",plag,"_draws=",draws+burnin,".rda"))
}else{
  source("./scripts/forecasting.R")
  save(dataset_est, coef, forecasts, file=paste0("./saves/forecasting_diff=",diff,"_plag=",plag,"_draws=",draws+burnin,".rda"))
}

# 4) Baseline Model

# Linear
if(file.exists(paste0("./saves/est_linmod_diff=",diff,"_plag=",plag,"_scale=",as.numeric(do_scale),"_draws=",draws+burnin,".rda"))){
  load(paste0("./saves/est_linmod_diff=",diff,"_plag=",plag,"_scale=",as.numeric(do_scale),"_draws=",draws+burnin,".rda"))
}else{
  source("./scripts/est_linmod.R")
  save(run_var, irfvar_chol, irfvar_ext, var_conv,
       file=paste0("./saves/est_linmod_diff=",diff,"_plag=",plag,"_scale=",as.numeric(do_scale),"_draws=",draws+burnin,".rda"))
}

# Threshold
if(file.exists(paste0("./saves/est_thrshmod_diff=",diff,"_plag=",plag,"_scale=",as.numeric(do_scale),"_draws=",draws+burnin,".rda"))){
  load(paste0("./saves/est_thrshmod_diff=",diff,"_plag=",plag,"_scale=",as.numeric(do_scale),"_draws=",draws+burnin,".rda"))
}else{
  source("./scripts/est_thrshmod.R")
  save(run_tvar, irftvar_chol, irftvar_ext, irftvar_robust, tvar_conv_reg1, tvar_conv_reg2,
       file=paste0("./saves/est_thrshmod_diff=",diff,"_plag=",plag,"_scale=",as.numeric(do_scale),"_draws=",draws+burnin,".rda"))
}

# 5) Extending the Information Set

for(q in q_est){
  # Linear Model
  if(file.exists(paste0("./saves/est_linextinfo_diff=",diff,"_plag=",plag,"_scale=",as.numeric(do_scale),"_q=",q,"_draws=",draws+burnin,".rda"))){
    load(paste0("./saves/est_linextinfo_diff=",diff,"_plag=",plag,"_scale=",as.numeric(do_scale),"_q=",q,"_draws=",draws+burnin,".rda"))
  }else{
    source("./scripts/est_linextinfo.R")
    save(run_varext, irfvarext_chol, irfvarext_ext, varext_conv,
         file=paste0("./saves/est_linextinfo_diff=",diff,"_plag=",plag,"_scale=",as.numeric(do_scale),"_q=",q,"_draws=",draws+burnin,".rda"))
  }
  assign(paste0("run_varext_",q), run_varext)
  assign(paste0("irfvarext_chol_",q), irfvarext_chol)
  assign(paste0("irfvarext_ext_",q), irfvarext_ext)
  assign(paste0("varext_conv_",q), varext_conv)
  
  # Threshold Model
  if(file.exists(paste0("./saves/est_thrshextinfo_diff=",diff,"_plag=",plag,"_scale=",as.numeric(do_scale),"_q=",q,"_draws=",draws+burnin,".rda"))){
    load(paste0("./saves/est_thrshextinfo_diff=",diff,"_plag=",plag,"_scale=",as.numeric(do_scale),"_q=",q,"_draws=",draws+burnin,".rda"))
  }else{
    source("./scripts/est_thrshextinfo.R")
    save(run_tvarext, irftvarext_chol, irftvarext_ext, tvarext_conv_reg1, tvarext_conv_reg2,
         file=paste0("./saves/est_thrshextinfo_diff=",diff,"_plag=",plag,"_scale=",as.numeric(do_scale),"_q=",q,"_draws=",draws+burnin,".rda"))
  }
  assign(paste0("run_tvarext_",q), run_tvarext)
  assign(paste0("irftvarext_chol_",q), irftvarext_chol)
  assign(paste0("irftvarext_ext_",q), irftvarext_ext)
  assign(paste0("tvarext_conv_reg1_",q), tvarext_conv_reg1)
  assign(paste0("tvarext_conv_reg2_",q), tvarext_conv_reg2)
}
rm(run_varext, irfvarext_chol, irfvarext_ext, varext_conv, run_tvarext, irftvarext_chol, irftvarext_ext, tvarext_conv_reg1, tvarext_conv_reg2)

# 6a) Robustness

# Linear Model
if(file.exists(paste0("./saves/est_linrob_diff=",diff,"_plag=",plag,"_scale=",as.numeric(do_scale),"_draws=",draws+burnin,".rda"))){
  load(paste0("./saves/est_linrob_diff=",diff,"_plag=",plag,"_scale=",as.numeric(do_scale),"_draws=",draws+burnin,".rda"))
}else{
  source("./scripts/est_linrob.R")
  save(irfvarrob_chol,
       file=paste0("./saves/est_linrob_diff=",diff,"_plag=",plag,"_scale=",as.numeric(do_scale),"_draws=",draws+burnin,".rda"))
}

# Threshold Model
if(file.exists(paste0("./saves/est_thrshrob_diff=",diff,"_plag=",plag,"_scale=",as.numeric(do_scale),"_draws=",draws+burnin,".rda"))){
  load(paste0("./saves/est_thrshrob_diff=",diff,"_plag=",plag,"_scale=",as.numeric(do_scale),"_draws=",draws+burnin,".rda"))
}else{
  source("./scripts/est_thrshrob.R")
  save(irftvarrob_chol,
       file=paste0("./saves/est_thrshrob_diff=",diff,"_plag=",plag,"_scale=",as.numeric(do_scale),"_draws=",draws+burnin,".rda"))
}

# 6b) Robustness - Internal Instruments Approach

# Linear
if(file.exists(paste0("./saves/est_linintinstr_diff=",diff,"_plag=",plag,"_scale=",as.numeric(do_scale),"_draws=",draws+burnin,".rda"))){
  load(paste0("./saves/est_linintinstr_diff=",diff,"_plag=",plag,"_scale=",as.numeric(do_scale),"_draws=",draws+burnin,".rda"))
}else{
  source("./scripts/est_linintinstr.R")
  save(run_var, irfvar_int,
       file=paste0("./saves/est_linintinstr_diff=",diff,"_plag=",plag,"_scale=",as.numeric(do_scale),"_draws=",draws+burnin,".rda"))
}

# Threshold
if(file.exists(paste0("./saves/est_thrshintinstr_diff=",diff,"_plag=",plag,"_scale=",as.numeric(do_scale),"_draws=",draws+burnin,".rda"))){
  load(paste0("./saves/est_thrshintinstr_diff=",diff,"_plag=",plag,"_scale=",as.numeric(do_scale),"_draws=",draws+burnin,".rda"))
}else{
  source("./scripts/est_thrshintinstr.R")
  save(run_tvar, irftvar_int,
       file=paste0("./saves/est_thrshintinstr_diff=",diff,"_plag=",plag,"_scale=",as.numeric(do_scale),"_draws=",draws+burnin,".rda"))
}

####### ONLY ADDITIONAL NOT IN THE PAPER
####### ATTENTION: VERY TIME-CONSUMING

# Robustness -  Ordering in the VAR / TVAR

# if(file.exists(paste0("./saves/rob_ordering_diff=",diff,"_plag=",plag,"_draws=",draws+burnin,".rda"))){
#   load(paste0("./saves/rob_ordering_diff=",diff,"_plag=",plag,"_draws=",draws+burnin,".rda"))
# }else{
#   source("./scripts/rob_ordering.R")
#   save(corr_var, corr_tvar,
#        file=paste0("./saves/rob_ordering_diff=",diff,"_plag=",plag,"_draws=",draws+burnin,".rda"))
# }

# 7) Figures

####################################
## Figure 1                       ##
####################################

png("./MS19-470Fig1.png", type = "cairo", width = width, height = height_higher, res = 300)
par(fig = c(0,1,0,1), mfrow=c(1,1), mar=c(3,2,1,1))
baat10_plot = dataset_est$BAAT10[(diff+1):nrow(dataset_est)]
plot.ts(baat10_plot, xaxt="n", yaxt="n",
        panel.first = rect(nbermat_common[,1], nbermat_common[,2], nbermat_common[,3], nbermat_common[,4], col='grey80', border=NA), 
        xlab = "", ylab = "",lwd=3, ylim = c(0.5,8.4))
lines(forecasts[,"DE"],col="grey30",lty=2,lwd=2)
rect(xleft = which(begin_zoom_m == time_sample), ybottom = 1.3, 
     xright = which(end_zoom_m == time_sample), ytop = 4.15, col=NA,
     border = TRUE, lty = 2)
axis(2, lwd = 2, cex = 1.5, las=2)
axis(1, at  = seq(1,Traw,by=20), 
     labels = format(time_sample[seq(1,Traw,by=20)],"%Y"), lwd = 2, cex = 1.5)
box(lwd=2)
par(fig = c(0.24,0.78,0.45,0.96), new = T) 
plot.ts(baat10_plot[time_sample%in%time_zoom], axes=FALSE,
        panel.first = rect(nbermat_zoom[,1], nbermat_zoom[,2], nbermat_zoom[,3],
                           nbermat_zoom[,4], col='grey80', border=NA), 
        xlab = "", ylab = "",lwd=3,ylim=c(1.4,4.2),
        main = "Close-up", cex.main = 1)
lines(forecasts[time_sample%in%time_zoom,"DE"],col="grey30",lty=2,lwd=2)
axis(1,lty=2, at = c(-1,100), labels = c("",""), lwd.ticks = 0)
axis(3,lty=2, at = c(-1,100), labels = c("",""), lwd.ticks = 0)
axis(2,lty=2, at = seq(from = 1.3, to = 4.4, by = 0.5), las=2)
axis(4,lty=2, at = c(1.3,4.4), labels = c("",""), lwd.ticks = 0)
dev.off()

####################################
## Figure 2                       ##
####################################

png("./MS19-470Fig2a.png", type = "cairo", width = width, height = height/2, res = 300)
par(mfrow=c(1,M), mar=c(2.2,3,2,1))
for(mm in 1:M){
  plot.ts(irfvar_ext[4,mm,], col = "black", ylim = range(irfvar_ext[,mm,]),
          xaxt = "n", yaxt = "n", axes=FALSE,
          main = varnames_plot[mm],
          lty = 5, lwd = 3)
  polygon(c(1:nhor,rev(1:nhor)), c(irfvar_ext[1,mm,],rev(irfvar_ext[7,mm,])),
          col = "grey90", border=NA)
  polygon(c(1:nhor,rev(1:nhor)), c(irfvar_ext[2,mm,],rev(irfvar_ext[6,mm,])),
          col = "grey70", border=NA)
  polygon(c(1:nhor,rev(1:nhor)), c(irfvar_ext[3,mm,],rev(irfvar_ext[5,mm,])), 
          col = "grey30", border=NA)
  lines(irfvar_ext[4,mm,], col="black", lwd=3)
  axis(1, at = seq(1,nhor+1,by=12), labels = seq(0,nhor,by=12), lwd=2)
  axis(2, lwd=2, las=2)
  abline(h=0, col = "black", lty=2, lwd=2)
  box(lwd=2, bty="l")
}
dev.off()

png("./MS19-470Fig2b.png", type = "cairo", width = width, height = height/2, res = 300)
par(mfrow=c(1,M), mar=c(2.2,3,2,1))
for(mm in 1:M){
  plot.ts(irfvar_chol[4,mm,], col = "black", ylim = range(irfvar_chol[,mm,]),
          axes=FALSE, xlab = "", ylab = "",
          main = varnames_plot[mm],
          lty = 5, lwd = 3)
  polygon(c(1:nhor,rev(1:nhor)), c(irfvar_chol[1,mm,],rev(irfvar_chol[7,mm,])),
          col = "grey90", border=NA)
  polygon(c(1:nhor,rev(1:nhor)), c(irfvar_chol[2,mm,],rev(irfvar_chol[6,mm,])),
          col = "grey70", border=NA)
  polygon(c(1:nhor,rev(1:nhor)), c(irfvar_chol[3,mm,],rev(irfvar_chol[5,mm,])), 
          col = "grey30", border=NA)
  lines(irfvar_chol[4,mm,], col="black", lty=1, lwd=3)
  axis(1, at = seq(1,nhor+1,by=12), labels = seq(0,nhor,by=12), lwd=2)
  axis(2, lwd=2, las=2)
  abline(h=0, col = "black", lty=2,lwd=2)
  box(lwd=2, bty="l")
}
dev.off()

####################################
## Figure 3                       ##
####################################

png("./MS19-470Fig3.png", type = "cairo", width = width, height = height_higher, res = 300)
len <- nrow(run_tvar$args$Zraw[-c(1:plag),,drop=FALSE])
regmat <- matrix(NA, len, 4)
colnames(regmat) <- c("xstart","ystart","xend","yend")
regmat[,"xstart"] <- (0:(len-1))/len
regmat[,"ystart"] <- rep(0,len)
regmat[,"xend"]   <- (1:len)/len
regmat[,"yend"]   <- apply(run_tvar$S-1,2,mean)

par(mfrow=c(1,1), mar=c(2,2,1,2))
plot.ts(run_tvar$args$Zraw[-c(1:plag),,drop=FALSE], axes=FALSE,
        panel.first = rect(xleft = regmat[,1], ybottom = regmat[,2], 
                           xright = regmat[,3], ytop = regmat[,4],
                           col="grey70",border=NA), lwd=3, ylim=range(run_tvar$args$Zraw[-c(1:plag),,drop=FALSE]))
segments(1,quantile(run_tvar$gamma, .50), len, col="black", lty=1, lwd=2)
segments(1,quantile(run_tvar$gamma, .95), len, col="black", lty=2, lwd=1.5)
segments(1,quantile(run_tvar$gamma, .05), len, col="black", lty=2, lwd=1.5)
axis(2, lwd = 2, cex = 1.5, las=2)
axis(1, at=seq(1,len,by=20), labels=format(time_sample[seq(plag+1,Traw,by=20)],"%Y"),
     lwd = 2, cex = 1.5)
axis(4, at=seq(1,6,length.out=5), labels = c(0,0.25,0.50,0.75,1),
     lwd = 2, cex = 1.5)
box(lwd=2)
dev.off()

####################################
## Figure 4                       ##
####################################

png("./MS19-470Fig4.png", type = "cairo", width = width, height = height, res = 300)
par(mfrow=c(h,M), mar = c(2.2,3,2,1))
for(hh in 1:h){
  for(mm in 1:M) {
    if(mm==1) par(mar=c(2.2,4.5,2,1)) else par(mar=c(2,2,2,1))
    ylim1 <- range(irftvar_ext[,mm,,])
    plot.ts(irftvar_ext[4,mm,,hh], ylim = ylim1, axes=FALSE,
            xlab = "", ylab = "", lty = 5, lwd=3, main=varnames_plot[mm])
    polygon(c(1:nhor,rev(1:nhor)), c(irftvar_ext[1,mm,,hh], rev(irftvar_ext[7,mm,,hh])),
            col = "grey85", border=NA)
    polygon(c(1:nhor,rev(1:nhor)), c(irftvar_ext[2,mm,,hh], rev(irftvar_ext[6,mm,,hh])),
            col = "grey65", border=NA)
    polygon(c(1:nhor,rev(1:nhor)), c(irftvar_ext[3,mm,,hh], rev(irftvar_ext[5,mm,,hh])), 
            col = "grey30", border=NA)
    lines(irftvar_ext[4,mm,,hh], lty=1, lwd=2.5)
    abline(h=0, col = "black", lty=2, lwd=2)
    axis(1, at = seq(1,nhor,by=12), labels = seq(0,nhor,by=12), lwd=2)
    axis(2, lwd=2, las=2)
    if(mm==1){
      if(hh==1)  {mtext("Optimistic Credit Regime", side=2, padj=-2.2)}
      if(hh==2) {mtext("Pessimistic Credit Regime", side=2, padj=-2.2)}
    }
    box(lwd=2, bty="l")
  }
}
dev.off()

####################################
## Figure 5                       ##
####################################

png("./MS19-470Fig5.png", type = "cairo", width = width, height = height, res = 300)
par(mfrow=c(h,M), mar = c(2.2,3,2,1))
for(hh in 1:h){
  for(mm in 1:M) {
    if(mm==1) par(mar=c(2.2,4.5,2,1)) else par(mar=c(2,2,2,1))
    ylim1 <- range(irftvar_chol[,mm,,])
    plot.ts(irftvar_chol[4,mm,,hh], ylim = ylim1, axes=FALSE,
            xlab = "", ylab = "", lty = 5, lwd=3, main=varnames_plot[mm])
    polygon(c(1:nhor,rev(1:nhor)), c(irftvar_chol[1,mm,,hh],rev(irftvar_chol[7,mm,,hh])),
            col = "grey85", border=NA)
    polygon(c(1:nhor,rev(1:nhor)), c(irftvar_chol[2,mm,,hh],rev(irftvar_chol[6,mm,,hh])),
            col = "grey65", border=NA)
    polygon(c(1:nhor,rev(1:nhor)), c(irftvar_chol[3,mm,,hh],rev(irftvar_chol[5,mm,,hh])), 
            col = "grey30", border=NA)
    lines(irftvar_chol[4,mm,,hh], lty=1, lwd=3)
    abline(h=0, col = "black", lty=2, lwd=2)
    axis(1, at = seq(1,nhor,by=12), labels = seq(0,nhor,by=12), lwd=2)
    axis(2, lwd=2, las=2)
    if(mm==1){
      if(hh==1)  {mtext("Optimistic Credit Regime", side=2, padj=-2.2)}
      if(hh==2) {mtext("Pessimistic Credit Regime", side=2, padj=-2.2)}
    }
    box(lwd=2, bty="l")
  }
}
dev.off()

####################################
## Figure 6                       ##
####################################

# regime 1

png("./MS19-470Fig6.png", type = "cairo", width = width, height = height_higher, res = 300)
par(mfrow=c(r-1,M-1), mar=c(2,2,2,1))
heur  <- c("RE","ADA", "WTR", "STR", "LAA")
for(rr in 2:r){
  for(mm in 2:M){
    temp1 <- hist(irftvar_robust[,mm,1,rr], plot=FALSE, breaks=20)
    temp2 <- hist(irftvar_robust[,mm,1,1], plot=FALSE, breaks=20)
    if(mm==2) par(mar=c(2,4,2,1)) else par(mar=c(2,2,2,1))
    hist(irftvar_robust[,mm,1,rr], xlab = "", ylab = "", main = "", freq=FALSE, breaks=20,
         col=rgb(0.8,0.8,0.8,0.5), xlim=range(temp1$breaks,temp2$breaks), ylim=extendrange(c(temp1$density,temp2$density), f=c(0.01,.05)))
    if(rr==2) title(main = varnames_plot[mm])
    hist(irftvar_robust[,mm,1,1],col=rgb(0.05,0.05,0.05,0.5),add=T,freq=FALSE,breaks=20)
    abline(v=0, col="black", lty=2, lwd=2)
    if(mm==2) mtext(heur[rr-1], side=2, padj=-2.2)
    box(lwd=2)
  }
}
dev.off()

# regime 2 - no differences as pointed out in the paper

# png("./MS19-470Fig6_regime2.png", type = "cairo", width = width, height = height_higher, res = 300)
# par(mfrow=c(r-1,M-1), mar=c(2,2,2,1))
# heur  <- c("RE","ADA", "WTR", "STR", "LAA")
# for(rr in 2:r){
#   for(mm in 2:M){
#     temp1 <- hist(irftvar_robust[,mm,2,rr], plot=FALSE, breaks=20)
#     temp2 <- hist(irftvar_robust[,mm,2,1], plot=FALSE, breaks=20)
#     if(mm==2) par(mar=c(2,4,2,1)) else par(mar=c(2,2,2,1))
#     hist(irftvar_robust[,mm,2,rr], xlab = "", ylab = "", main = "", freq=FALSE, breaks=20,
#          col=rgb(0.8,0.8,0.8,0.5), xlim=range(temp1$breaks,temp2$breaks), ylim=extendrange(c(temp1$density,temp2$density), f=c(.01,.05)))
#     if(rr==2) title(main = varnames_plot[mm-1])
#     hist(irftvar_robust[,mm,2,1], col=rgb(0.05,0.05,0.05,0.5), add=T, freq=FALSE, breaks=20)
#     abline(v=0, col="red", lty=2,lwd=2)
#     if(mm==2) mtext(heur[rr-1], side=2, padj=-2.2)
#     box(lwd=2)
#   }
# }
# dev.off()

####################################
## Figure 7                       ##
####################################

png("./MS19-470Fig7.png", type = "cairo", width = width, height = height, res = 300)
par(mfrow=c(h,M), mar = c(2.2,3,2,1))
for(hh in 1:h){
  for(mm in 1:M) {
    if(mm==1) par(mar=c(2.2,4.5,2,1)) else par(mar=c(2,2,2,1))
    ylim1 <- range(irftvarext_ext_3[,mm,,],irftvar_ext[,mm,,])
    plot.ts(irftvarext_ext_3[4,mm,,hh], ylim = ylim1, axes=FALSE,
            xlab = "", ylab = "", lty = 5, lwd=3, main=varnames_plot[mm])
    polygon(c(1:nhor,rev(1:nhor)), c(irftvarext_ext_3[1,mm,,hh], rev(irftvarext_ext_3[7,mm,,hh])),
            col = "grey85", border=NA)
    polygon(c(1:nhor,rev(1:nhor)), c(irftvarext_ext_3[2,mm,,hh], rev(irftvarext_ext_3[6,mm,,hh])),
            col = "grey65", border=NA)
    polygon(c(1:nhor,rev(1:nhor)), c(irftvarext_ext_3[3,mm,,hh], rev(irftvarext_ext_3[5,mm,,hh])), 
            col = "grey30", border=NA)
    lines(irftvarext_ext_3[4,mm,,hh], lty=1, lwd=3)
    abline(h=0, col = "black", lty=2, lwd=2)
    axis(1, at = seq(1,nhor,by=12), labels = seq(0,nhor,by=12), lwd=2)
    axis(2, lwd=2, las=2)
    if(mm==1){
      if(hh==1)  {mtext("Optimistic Credit Regime", side=2, padj=-2.2)}
      if(hh==2) {mtext("Pessimistic Credit Regime", side=2, padj=-2.2)}
    }
    box(lwd=2)
    # add baseline model
    lines(irftvar_ext[4,mm,,hh], col="#CD661D", lwd=2)
    lines(irftvar_ext[7,mm,,hh], col="#EE7621", lwd=2, lty=2)
    lines(irftvar_ext[1,mm,,hh], col="#EE7621", lwd=2, lty=2)
    lines(irftvar_ext[6,mm,,hh], col="#EE7621", lwd=2, lty=3)
    lines(irftvar_ext[2,mm,,hh], col="#EE7621", lwd=2, lty=3)
    lines(irftvar_ext[3,mm,,hh], col="#EE7621", lwd=2, lty=4)
    lines(irftvar_ext[5,mm,,hh], col="#EE7621", lwd=2, lty=4)
    if(mm == 5 && hh==2) legend("bottomright", c("Extended Model", "Baseline Model"), col=c("black", "#CD661D"), lwd=c(2,2))
  }
}
dev.off()

####################################
## Figure F1                      ##
####################################

png("./MS19-470FigF1.png", type = "cairo", width = width, height = height/2, res = 300)
par(mfrow=c(1,M), mar=c(2.2,2,2,1))
for(mm in 1:M){
  plot.ts(irfvarrob_chol[[2]][4,mm,], col = "black", ylim = range(irfvarrob_chol[[2]][,mm,]),
          axes=FALSE, xlab = "", ylab = "",
          main = varnames_plot[order[[2]][mm]],
          lty = 5, lwd = 3)
  polygon(c(1:nhor,rev(1:nhor)), c(irfvarrob_chol[[2]][1,mm,],rev(irfvarrob_chol[[2]][7,mm,])),
          col = "grey85", border=NA)
  polygon(c(1:nhor,rev(1:nhor)), c(irfvarrob_chol[[2]][2,mm,],rev(irfvarrob_chol[[2]][6,mm,])),
          col = "grey65", border=NA)
  polygon(c(1:nhor,rev(1:nhor)), c(irfvarrob_chol[[2]][3,mm,],rev(irfvarrob_chol[[2]][5,mm,])), 
          col = "grey30", border=NA)
  lines(irfvarrob_chol[[2]][4,mm,], col="black", lty=1, lwd=3)
  axis(1, at = seq(1,nhor+1,by=12), labels = seq(0,nhor,by=12), lwd=2)
  axis(2, lwd=2, las=2)
  abline(h=0, col = "black", lty=2,lwd=2)
  box(lwd=2, bty="l")
}
dev.off()

####################################
## Figure F2                      ##
####################################

png("./MS19-470FigF2.png", type = "cairo", width = width, height = height, res = 300)
par(mfrow=c(h,M), mar = c(2.2,3,2,1))
for(hh in 1:h){
  for(mm in 1:M) {
    if(mm==1) par(mar=c(2.2,4.5,2,1)) else par(mar=c(2,2,2,1))
    ylim1 <- range(irftvarrob_chol[[2]][,mm,,])
    plot.ts(irftvarrob_chol[[2]][4,mm,,hh], ylim = ylim1, axes=FALSE,
            xlab = "", ylab = "", lty = 5, lwd=3, main=varnames_plot[order[[2]][mm]])
    polygon(c(1:nhor,rev(1:nhor)), c(irftvarrob_chol[[2]][1,mm,,hh],rev(irftvarrob_chol[[2]][7,mm,,hh])),
            col = "grey85", border=NA)
    polygon(c(1:nhor,rev(1:nhor)), c(irftvarrob_chol[[2]][2,mm,,hh],rev(irftvarrob_chol[[2]][6,mm,,hh])),
            col = "grey65", border=NA)
    polygon(c(1:nhor,rev(1:nhor)), c(irftvarrob_chol[[2]][3,mm,,hh],rev(irftvarrob_chol[[2]][5,mm,,hh])), 
            col = "grey30", border=NA)
    lines(irftvarrob_chol[[2]][4,mm,,hh], lty=1, lwd=3)
    abline(h=0, col = "black", lty=2, lwd=2)
    axis(1, at = seq(1,nhor,by=12), labels = seq(0,nhor,by=12), lwd=2)
    axis(2, lwd=2, las=2)
    if(mm==1){
      if(hh==1)  {mtext("Optimistic Credit Regime", side=2, padj=-2.2)}
      if(hh==2) {mtext("Pessimistic Credit Regime", side=2, padj=-2.2)}
    }
    box(lwd=2, bty="l")
  }
}
dev.off()

####################################
## Figure F3                      ##
####################################

png("./MS19-470FigF3.png", type = "cairo", width = width, height = height, res = 300)
par(mfrow=c(h,M), mar = c(2.2,3,2,1))
for(hh in 1:h){
  for(mm in 1:M) {
    if(mm==1) par(mar=c(2.2,4.5,2,1)) else par(mar=c(2,2,2,1))
    ylim1 <- range(irftvarext_ext_7[,mm,,],irftvar_ext[,mm,,])
    plot.ts(irftvarext_ext_7[4,mm,,hh], ylim = ylim1, axes=FALSE,
            xlab = "", ylab = "", lty = 5, lwd=3, main=varnames_plot[mm])
    polygon(c(1:nhor,rev(1:nhor)), c(irftvarext_ext_7[1,mm,,hh], rev(irftvarext_ext_7[7,mm,,hh])),
            col = "grey85", border=NA)
    polygon(c(1:nhor,rev(1:nhor)), c(irftvarext_ext_7[2,mm,,hh], rev(irftvarext_ext_7[6,mm,,hh])),
            col = "grey65", border=NA)
    polygon(c(1:nhor,rev(1:nhor)), c(irftvarext_ext_7[3,mm,,hh], rev(irftvarext_ext_7[5,mm,,hh])), 
            col = "grey30", border=NA)
    lines(irftvarext_ext_7[4,mm,,hh], lty=1, lwd=3)
    abline(h=0, col = "black", lty=2, lwd=2)
    axis(1, at = seq(1,nhor,by=12), labels = seq(0,nhor,by=12), lwd=2)
    axis(2, lwd=2, las=2)
    if(mm==1){
      if(hh==1)  {mtext("Optimistic Credit Regime", side=2, padj=-2.2)}
      if(hh==2) {mtext("Pessimistic Credit Regime", side=2, padj=-2.2)}
    }
    box(lwd=2, bty="l")
    # add baseline model
    lines(irftvar_ext[4,mm,,hh], col="#CD661D", lwd=2)
    lines(irftvar_ext[7,mm,,hh], col="#EE7621", lwd=2, lty=2)
    lines(irftvar_ext[1,mm,,hh], col="#EE7621", lwd=2, lty=2)
    lines(irftvar_ext[6,mm,,hh], col="#EE7621", lwd=2, lty=3)
    lines(irftvar_ext[2,mm,,hh], col="#EE7621", lwd=2, lty=3)
    lines(irftvar_ext[3,mm,,hh], col="#EE7621", lwd=2, lty=4)
    lines(irftvar_ext[5,mm,,hh], col="#EE7621", lwd=2, lty=4)
    if(mm == 5 && hh==2) legend("bottomright", c("Extended Model", "Baseline Model"), col=c("black", "#CD661D"), lwd=c(2,2))
  }
}
dev.off()

####################################
## Figure F4                      ##
####################################

varnames_plot_int <- c("Instrument",varnames_plot)

png("./MS19-470FigF4.png", type = "cairo", width = width, height = height/2, res = 300)
par(mfrow=c(1,M+1), mar=c(2.2,3,2,1))
for(mm in 1:(M+1)){
  plot.ts(irfvar_int[4,mm,], col = "black", ylim = range(irfvar_int[,mm,]),
          axes=FALSE, xlab = "", ylab = "",
          main = varnames_plot_int[mm],
          lty = 5, lwd = 3)
  polygon(c(1:nhor,rev(1:nhor)), c(irfvar_int[1,mm,],rev(irfvar_int[7,mm,])),
          col = "grey85", border=NA)
  polygon(c(1:nhor,rev(1:nhor)), c(irfvar_int[2,mm,],rev(irfvar_int[6,mm,])),
          col = "grey65", border=NA)
  polygon(c(1:nhor,rev(1:nhor)), c(irfvar_int[3,mm,],rev(irfvar_int[5,mm,])), 
          col = "grey30", border=NA)
  lines(irfvar_int[4,mm,], col="black", lty=1, lwd=3)
  axis(1, at = seq(1,nhor+1,by=12), labels = seq(0,nhor,by=12), lwd=2)
  axis(2, lwd=2, las=2)
  abline(h=0, col = "black", lty=2,lwd=2)
  box(lwd=2, bty="l")
}
dev.off()

####################################
## Figure F5                      ##
####################################

varnames_plot_int <- c("Instrument",varnames_plot)

png("./MS19-470FigF5.png", type = "cairo", width = width, height = height, res = 300)
par(mfrow=c(h,M+1), mar = c(2.2,3,2,1))
for(hh in 1:h){
  for(mm in 1:(M+1)) {
    if(mm==1) par(mar=c(2.2,4.5,2,1)) else par(mar=c(2,2,2,1))
    ylim1 <- range(irftvar_int[,mm,,])
    plot.ts(irftvar_int[4,mm,,hh], ylim = ylim1, axes=FALSE,
            xlab = "", ylab = "", lty = 5, lwd=3, main=varnames_plot_int[mm])
    polygon(c(1:nhor,rev(1:nhor)), c(irftvar_int[1,mm,,hh],rev(irftvar_int[7,mm,,hh])),
            col = "grey80", border=NA)
    polygon(c(1:nhor,rev(1:nhor)), c(irftvar_int[2,mm,,hh],rev(irftvar_int[6,mm,,hh])),
            col = "grey60", border=NA)
    polygon(c(1:nhor,rev(1:nhor)), c(irftvar_int[3,mm,,hh],rev(irftvar_int[5,mm,,hh])), 
            col = "grey40", border=NA)
    lines(irftvar_int[4,mm,,hh], lty=1, lwd=3)
    abline(h=0, col = "black", lty=2, lwd=2)
    axis(1, at = seq(1,nhor,by=12), labels = seq(0,nhor,by=12), lwd=2)
    axis(2, lwd=2, las=2)
    if(mm==1){
      if(hh==1)  {mtext("Optimistic Credit Regime", side=2, padj=-2.2)}
      if(hh==2) {mtext("Pessimistic Credit Regime", side=2, padj=-2.2)}
    }
    box(lwd=2, bty="l")
  }
}
dev.off()

# 8) Tables

####################################
## Table E1                       ##
####################################

lapply(var_conv, function(l) round(l,3))
lapply(tvar_conv_reg1, function(l) round(l,3))
lapply(tvar_conv_reg2, function(l) round(l,3))
lapply(varext_conv_3,function(l)round(l,3))
lapply(tvarext_conv_reg1_3, function(l) round(l,3))
lapply(tvarext_conv_reg2_3, function(l) round(l,3))




