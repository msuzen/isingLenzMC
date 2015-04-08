#  Plot Data
#  (c) 2013, 2014 by Dr.Mehmet Suzen
#  GPLv3 or higher
#

source("isingErgodicity.R"); #

# Plots
#   Two sets: metropolis and glauber 120 plots total.
#   1. Fixed N (5 cases)  (total of 40 plots)
#      a. Fixed T vary H (4 cases)
#      b. Fixed H vary T (4 cases)
#   2. Given (T, H) vary N (20 cases)


getFileName <- function(dirData, dynamicsName, N, ikBT, H) {
  dataFileName <- paste(dirData, '/', dynamicsName, "TmMetric-n=", sprintf('%.3e', N),
                        "ikBT=",  sprintf('%.3f', ikBT),
                        "H=",  sprintf('%.3f', H), ".dat",
                        sep="")
  dataFileName
}

# Parameters
N    <- c(32, 64, 128, 256, 512)
ikBT <- c(0.5, 1.0, 1.5, 2.0)
#H    <- c(-1.0, -0.50,  0.0, 0.50, 1.0)
H    <- c(-1.0, -0.50,  0.50, 1.0)

# 1(a)  Fixed N, kBT : Vary H
dynamicsName <- 'metropolis'; # switch names for additional plots
dynamicsName <- 'glauber';
dirData      <- 'data';
dirPlot      <- 'plots/varyH';
for(i in 1:length(N)) {
  for(j in 1:length(ikBT)) {
        plotName   <- paste(dirPlot, "/", dynamicsName, "VaryH-N=", sprintf('%.3e', N[i]), sprintf('ikBT=%.3f', ikBT[j]), ".png", sep="")
        print(plotName)
        png(filename=plotName);
        legendText <- c()
        fData      <- getFileName(dirData, dynamicsName, N[i], ikBT[j], H[1])
        aa         <- readDat(fData);
        plot(aa[,1], aa[1,2]/aa[,2], log="xy", type="l", lty=1, axes=FALSE, ann=FALSE) # col, pch not need while ploting lines
        # col, pch not need while ploting lines
        mainText<- paste(dynamicsName, " N=", sprintf('%.3e', N[i]), "ikBT=",  sprintf('%.3f', ikBT[j]))
        title(main=mainText , xlab="MC Time", ylab="Inverse Effective Ergodic Convergence Rate");
        axis(1)
        axis(2)
        box()
        legendText[1] <- paste("h=", sprintf('%.3e', H[1]));
    for(k in 2:length(H)) {
         aa <- matrix()
         fData <- getFileName(dirData, dynamicsName, N[i], ikBT[j], H[k]) 
         aa <- readDat(fData);
         lines(aa[,1], aa[1,2]/aa[,2], type="l", lty=k)
         legendText[k] <- paste("h=", sprintf('%.3e', H[k]));
 
     }
        legend(1, max(aa[1,2]/aa[,2])*0.8, legendText, lty=1:5)
        dev.off();
   }
}
# 1(b)  Fixed N, H : Vary kBT
dirData      <- 'data';
dirPlot      <- 'plots/varyT';
for(i in 1:length(N)) {
  for(k in 1:length(H)) {
        plotName   <- paste(dirPlot, "/", dynamicsName, "VaryT-N=", sprintf('%.3e', N[i]), sprintf('H=%.3f', H[k]), ".png", sep="")
        print(plotName)
        png(filename=plotName);
        legendText <- c()
        fData      <- getFileName(dirData, dynamicsName, N[i], ikBT[1], H[k])
        aa         <- readDat(fData);
        plot(aa[,1], aa[1,2]/aa[,2], log="xy", type="l", lty=1, axes=FALSE, ann=FALSE) # col, pch not need while ploting lines
        # col, pch not need while ploting lines
        mainText<- paste(dynamicsName, " N=", sprintf('%.3e', N[i]), "H=",  sprintf('%.3f', H[k]))
        title(main=mainText , xlab="MC Time", ylab="Inverse Effective Ergodic Convergence Rate");
        axis(1)
        axis(2)
        box()
        legendText[1] <- paste("1/kBT=", sprintf('%.3e', ikBT[1]));
    for(j in 2:length(ikBT)) {
         aa <- matrix()
         fData <- getFileName(dirData, dynamicsName, N[i], ikBT[j], H[k])
         aa <- readDat(fData);
         lines(aa[,1], aa[1,2]/aa[,2], type="l", lty=j)
         legendText[j] <- paste("1/kBT=", sprintf('%.3e', ikBT[j]));
     }
        legend(1, max(aa[1,2]/aa[,2])*0.65, legendText, lty=1:5)
        dev.off();
   }
}

# 2) Given T,H Vary N
dirData      <- 'data';
dirPlot      <- 'plots/varyN';
for(j in 1:length(ikBT)) {
  for(k in 1:length(H)) {
        plotName   <- paste(dirPlot, "/", dynamicsName, "VaryN-T=", sprintf('%.3e', ikBT[j]), sprintf('H=%.3f', H[k]), ".png", sep="")
        print(plotName)
        png(filename=plotName);
        legendText <- c()
        fData      <- getFileName(dirData, dynamicsName, N[1], ikBT[j], H[k])
        aa         <- readDat(fData);
        plot(aa[,1], aa[1,2]/aa[,2], log="xy", type="l", lty=1, axes=FALSE, ann=FALSE) # col, pch not need while ploting lines
        # col, pch not need while ploting lines
        mainText<- paste(dynamicsName, " ikBT=", sprintf('%.3e', ikBT[j]), "H=",  sprintf('%.3f', H[k]))
        title(main=mainText , xlab="MC Time", ylab="Inverse Effective Ergodic Convergence Rate");
        axis(1)
        axis(2)
        box()
        legendText[1] <- paste("N=", sprintf('%.3e', N[1]));
    for(i in 2:length(N)) {
         aa <- matrix()
         fData <- getFileName(dirData, dynamicsName, N[i], ikBT[j], H[k])
         aa <- readDat(fData);
         lines(aa[,1], aa[1,2]/aa[,2], type="l", lty=i)
         legendText[i] <- paste("N=", sprintf('%.3e', N[i]));
     }
        legend(1, max(aa[1,2]/aa[,2])*0.65, legendText, lty=1:5)
        dev.off();
   }
}

