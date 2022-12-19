###############################################
### Load and Build Dataset                  ###
### The Impact of Credit Market Sentiments  ###
### Maximilian Boeck                        ###
### 30/04/2021                              ###
###############################################

library(readxl)
library(dplyr)
library(lubridate)

# load data and build data set
dataset_full <- read.csv("./data/mccracken-monthly.csv", stringsAsFactors = FALSE)
dataset_tcode <- as.numeric(dataset_full[1,-1])
dataset_full <- dataset_full[-1,]
dataset_full <- dataset_full[-722,]
dataset_full$sasdate <- as.Date(dataset_full$sasdate, format = "%m/%d/%Y")

# add shadowrate
shadowrate  <- read_xlsx("./data/WuXiaShadowRate.xlsx", col_names = TRUE, sheet = "Data")
shadowrate  <- shadowrate[,-c(4:5)]
shadowrate  <- shadowrate[-672,]
colnames(shadowrate) <- c("time","FFRWX","SR")
shadowrate$time <- as.Date(shadowrate$time, format = "%Y-%m-%d")
dataset_full <- left_join(dataset_full,shadowrate,by=c("sasdate"="time"))
dataset_full$FFRWXSR <- c(dataset_full$FFRWX[1:600],dataset_full$SR[601:683],rep(NA,38))

# compute BAA spread
dataset_full$BAAT10 <- dataset_full$BAA - dataset_full$GS10

# subset data
dataset_est   <- subset(dataset_full, sasdate <= end_sample & sasdate >= ymd(as.Date(begin_sample)) %m-% months(diff))
dataset_train <- subset(dataset_full, sasdate <= end_sample)

# NBER Recession dates
nber_m <- read.csv("./data/USRECM.csv", stringsAsFactors = FALSE)
nber_m$DATE <- as.Date(nber_m$DATE, format = "%Y-%m-%d")
nber_m <- subset(nber_m, DATE >= "1959-01-01")
nber_m <- nber_m[-722,]
nber_mm <- subset(nber_m, DATE <= end_sample)
nber_mm <- subset(nber_mm, DATE >= begin_sample)
nber_mm$diff <- c(NA,diff(nber_mm$USREC))
nber_mm <- nber_mm[-1,]
nbermat <- matrix(NA, length(time_sample), 4)
colnames(nbermat) <- c("xstart","ystart","xend","yend")
nbermat[,2] <- 0.001
nbermat[,4] <- 0.999
nbermat[which(nber_mm$diff==1),1]  <- which(nber_mm$diff==1)/length(time_sample)
nbermat[which(nber_mm$diff==1),3] <- (which(nber_mm$diff==-1)-1)/length(time_sample)
nbermat <- nbermat[-which(is.na(nbermat[,1])),]

# zoom stuff
begin_zoom_m <- as.Date("1998-01-01")
begin_zoom_lag <- as.Date("1997-12-01")
end_zoom_m   <- as.Date("2004-12-01")
time_zoom   <- seq.Date(begin_zoom_m,end_zoom_m, by="1 month")

nber_zoom_mm <- subset(nber_m, DATE <= end_zoom_m)
nber_zoom_mm <- subset(nber_zoom_mm, DATE >= begin_zoom_lag)
nber_zoom_mm$diff <- c(NA,diff(nber_zoom_mm$USREC))
nber_zoom_mm <- nber_zoom_mm[-1,]
nbermat_zoom <- matrix(NA, length(time_zoom), 4)
colnames(nbermat_zoom) <- c("xstart","ystart","xend","yend")
nbermat_zoom[,2] <- 0.001
nbermat_zoom[,4] <- 0.999
nbermat_zoom[which(nber_zoom_mm$diff==1),1] <- which(nber_zoom_mm$diff==1)/length(time_zoom)
nbermat_zoom[which(nber_zoom_mm$diff==1),3] <- (which(nber_zoom_mm$diff==-1)-1)/length(time_zoom)
nbermat_zoom <- nbermat_zoom[-which(is.na(nbermat_zoom[,1])),,drop=F]

nbermat_common <- nbermat
nbermat_common[3:6,"yend"] <- 0.5

# delete unnecessary stuff
rm(shadowrate, nber_m, nber_mm, nber_zoom_mm)
