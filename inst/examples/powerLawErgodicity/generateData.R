#
#  Ergodic Dynamics of Ising Model 
#  R simulation functions, data generation
#  (c) 2013, 2014, 2015, 2016 by Dr.Mehmet Suzen
#  GPLv3 or higher
#

rm(list=ls())
source("functions.R"); 

#' 2015_ more data to measure the power laws.
#' 
library(matlab)

# Generate Data 
# Parameters for scaling law
N    <- c(512)
ikBT <- seq(0.5,2.0,0.025)
H    <- c(1.0)
#
J      <- 1.0
nstep  <- 2.0e6

#' sim parameter lengths
nn <- length(N)
kk <- length(ikBT)
hh <- length(H)

ising_ergo <- vector("list", nn*kk*hh)
l <- 1
tic();
for(i in 1:nn) {
  for(j in 1:kk) {
    for(k in 1:hh) {
      print(N[i])
      print(ikBT[j])
      print(H[k])
      # magnetisation metric
      ising_ergo[[l]]$N          <- N[i]
      ising_ergo[[l]]$ikBT       <- ikBT[j]
      ising_ergo[[l]]$H          <- H[k]
      ising_ergo[[l]]$metropolis <- runSim(ikBT[j], J, H[k], N[i], nstep, 1) # metropolis
      ising_ergo[[l]]$glauber    <- runSim(ikBT[j], J, H[k], N[i], nstep, 2) # glauber
      l <- l + 1
     }
   }
}
toc();

save(ising_ergo, file="ising_ergo.RData")

