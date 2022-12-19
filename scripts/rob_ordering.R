###############################################
### Robustness: Ordering in the VAR / TVAR  ###
### The Impact of Credit Market Sentiments  ###
### Maximilian Boeck                        ###
### 08/05/2021                              ###
###############################################

Yraw1 <- as.matrix(dataset_est[,vars])
Zraw1 <- as.matrix(dataset_est[,thrshvar])
M     <- ncol(Yraw1)

# transformations
for(mm in 1:M){
  Yraw1[,mm] <- transx(Yraw1[,mm], tcode=tcode[mm], lag=diff)
}
Yraw1 <- Yraw1[-c(1:diff),]
Zraw1 <- Zraw1[-c(1:diff),,drop=FALSE]
Yraw1 <- apply(Yraw1, 2, scale)
rownames(Yraw1)<-as.character(time_sample)

# get all permutations
orderings <- gtools::permutations(n = 5, r = M, v = 1:M)

for(oo in 1:nrow(orderings)){
  if(file.exists(paste0("./saves/modorder=",oo,"_diff=",diff,"_plag=",plag,"_draws=",draws+burnin,".rda"))){
    cat(paste0("Round: ", oo,".\n"))
    next
  }
  
  set.seed(571)
  
  order <- orderings[oo,]
  
  temp_var  <- try(bvar(Yraw = Yraw1[, order], plag = plag, nsave = draws, nburn = burnin, thin = thin, 
                        cons = TRUE, trend = FALSE, sv = FALSE, eigen = TRUE), silent=TRUE)
  temp_tvar <- try(btvar(Yraw = Yraw1[, order], plag = plag, d.min = 1, d.max = 4, Zraw = Zraw1, nsave = draws, nburn = burnin, thin = thin, 
                         cons = TRUE, trend = FALSE, sv = FALSE, eigen = TRUE), silent=TRUE)
  
  save(temp_var, temp_tvar, file=paste0("./saves/modorder=",oo,"_diff=",diff,"_plag=",plag,"_draws=",draws+burnin,".rda"))
  
  cat(paste0("Round: ", oo,".\n"))
}

# compute correlations
corr_var  <- numeric(length=nrow(orderings))
corr_tvar <- numeric(length=nrow(orderings))

load(paste0("./saves/modorder=1_diff=",diff,"_plag=",plag,"_draws=",draws+burnin,".rda"))
var_base<-temp_var
tvar_base<-temp_tvar

rowvec1 <- NULL
for(pp in 1:plag) rowvec1 <- c(rowvec1,orderings[1,]+(pp-1)*M)
rowvec1 <- c(rowvec1,66)
colvec1 <- orderings[1,]

for(rr in 1:nrow(orderings)){
  filepath <- paste0("./saves/modorder=",rr,"_diff=",diff,"_plag=",plag,"_draws=",draws+burnin,".rda")
  if(file.exists(filepath))
    load(filepath) else next
  
  var_comp<-temp_var
  tvar_comp<-temp_tvar
  
  rowvec2 <- NULL
  for(pp in 1:plag) rowvec2 <- c(rowvec2,orderings[rr,]+(pp-1)*M)
  rowvec2 <- c(rowvec2,66)
  colvec2 <- orderings[rr,]
  
  post1 <- try(apply(var_base$store$A_store,c(2,3),median), silent=TRUE)
  post2 <- try(apply(var_comp$store$A_store,c(2,3),median), silent=TRUE)
  post3 <- try(apply(tvar_base$store$A_store,c(2,3,4),median), silent=TRUE)
  post4 <- try(apply(tvar_base$store$A_store,c(2,3,4),median), silent=TRUE)
  
  if(is(post1,"try-error") | is(post2,"try-error") | is(post3,"try-error") | is(post4,"try-error")) next
  
  corr_var[rr]  <- cor(as.vector(post1[rowvec1,colvec1]), as.vector(post2[rowvec2,colvec2]))
  corr_tvar[rr] <- cor(as.vector(post3[rowvec1,colvec1,]), as.vector(post4[rowvec2,colvec2,]))
  
  print(paste0("Round: ", rr))
}

rm(Yraw1, Zraw1, M, orderings, oo, var_base, tvar_base, temp_var, temp_tvar, rowvec1, rowvec2, colvec1, colvec2,
   post1, post2, post3, post4, rr)

