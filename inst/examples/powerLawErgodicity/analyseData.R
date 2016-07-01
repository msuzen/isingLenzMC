#  Ergodic Dynamics of Ising Model 
#  R simulation functions, data analysis
#  (c) 2013, 2014, 2015, 2016 by Dr.Mehmet Suzen
#  GPLv3 or higher

rm(list=ls())

####################################################################
## Read Data
####################################################################
load("ising_ergo.RData")
library(igraph) # for power-law fit

# function to read data
get_omega_time <- function(ii, mg="metropolis") {
  data_current <- ising_ergo[[ii]]
  data_current$N
  data_current$ikBT
  data_current$H
  if(mg == "metropolis") {
    time_   <- data_current$metropolis[,1]
    omega_  <- data_current$metropolis[,2]
  }
  if(mg == "glauber") {
    time_   <- data_current$glauber[,1]
    omega_  <- data_current$glauber[,2]
  }
  ll     <- length(time) # sample by 100
  return(list(time_=time_, omega_=omega_, ikBT=data_current$ikBT))
}

####################################################################
## Compute exponents
####################################################################
# mg <- "metropolis" # run for each dynamics
mg          <- "glauber"
ll          <- length(ising_ergo)
pf_list     <- list()
time_list   <- list()
omega_list  <- list()
temp <- c()
for(i in 1:ll) {
  print(i)
  time_omega <- get_omega_time(i, mg=mg)
  temp[i]    <- time_omega$ikBT 
  time_      <- time_omega$time_
  omega_     <- time_omega$omega_
  time_list[[i]]  <- time_
  omega_list[[i]] <- omega_
  lo <- length(omega_)
  ixx <- 1:lo
  if(lo >10000) {
    ixx<-seq(1,lo,100)
  }
  pf            <- power.law.fit(omega_[ixx])
  pf_list[[i]]  <- pf
}
# alphas
alphas <- c()
xmins <- c()
ks_stats<- c()
for(i in 1:ll) {
 alphas[i]   <- pf_list[[i]]$alpha 
 xmins[i]    <- pf_list[[i]]$xmin
 ks_stats[i] <- pf_list[[i]]$KS.stat
}

if(mg == "glauber") { 
  naccept<-sapply(1:ll, function(i) length(ising_ergo[[i]]$glauber[,1]))
}
if(mg == "metropolis") { 
  naccept<-sapply(1:ll, function(i) length(ising_ergo[[i]]$metropolis[,1]))
}
df <- data.frame(temp=temp, alphas=alphas, xmins=xmins, ks_stats=ks_stats,  naccept=naccept)
write.table(df,file=paste(c(mg,".csv"),collapse = ""),sep=",")

####################################################################
## plot two curve
##   One with alpha >2.0 and one with alpha <2.0
####################################################################
ps.options(pointsize=25)
plotName <- paste(c(mg, ".eps"), collapse="")
csvName <- paste(c(mg, ".csv"), collapse="")
setEPS()
postscript(plotName)
par(mar=c(4,4,0.5,0.5), xaxs = "i", yaxs = "i")
# 
AT       <- c(0, 10, 100, 1000, 10000, 100000, 1000000)
AT.label <- expression(0, 10^1, 10^2, 10^3, 10^4, 10^5, 10^6) 
AT2       <- c(0, 0.1, 0.01, 0.001, 0.0001, 0.00001,0.000001)
AT2.label <- expression(0, 10^-1, 10^-2, 10^-3, 10^-4, 10^-5, 10^-6) 
legendText <- c()

ix_vec  <- which(alphas < 3.0)
ixx     <- ix_vec[1]
time_   <- time_list[[ixx]][-1] # omit 1
omega_  <- omega_list[[ixx]][-1]
  
plot(time_,  omega_, log="xy", type="l", lty=1, axes=FALSE, ann=FALSE,
     lwd=5,  cex.lab=5.5, cex.axis=5.5, cex.main=5.5, cex.sub=5.5) # 
title(xlab="MC Time", ylab="Inverse Rate");
axis(1, at=AT, labels=AT.label)
axis(2, at=AT2, labels=AT2.label)
box()
legendText[1] <- paste("1/kT=", sprintf('%.1f', temp[1]));
pf_list[[1]] <- pf
k <- 1 # legend
for(i in ix_vec[-1]) {
  # print(length(omega_))
  if(i%%8 == 0) { # plot at every 10
    k <- k +1
    time_      <- time_list[[i]]
    omega_     <- omega_list[[i]]
    lines(time_, omega_, type="l", lty=i, lwd=5)
    legendText[k] <- paste("1/kT=", sprintf('%.1f', temp[i]));
  }
}

legend("bottomleft", legendText, lty=1:k, bty="n", lwd=5)
dev.off()

# Scaling exponents vs. Temperature
ps.options(pointsize=25)
plotName <- paste(c(mg, "_scaling_exponents.eps"), collapse="")
#setEPS()
postscript(plotName)
plot(temp[ix_vec],  alphas[ix_vec], log="xy", type="p", lty=1, axes=FALSE, ann=FALSE,
     lwd=5,  cex.lab=5.5, cex.axis=5.5, cex.main=5.5, cex.sub=5.5) # 
title(xlab="Temperature", ylab="Scaling Exponents");
axis(1)
axis(2)
box()
dev.off()
