#
#  Ergodic Dynamics of Ising Model 
#  R simulation functions, data generation
#  (c) 2013, 2014 by Dr.Mehmet Suzen
#  GPLv3 or higher
#

source("isingErgodicity.R"); #


# Generate Data
# Parameters
N    <- c(32, 64, 128, 256, 512)
ikBT <- c(0.5, 1.0, 1.5, 2.0)
H    <- c(-1.0, -0.50,  0.0, 0.50, 1.0)
#
J      <- 1.0
nstep  <- 500000

for(i in 1:length(N)) {
  for(j in 1:length(ikBT)) {
    for(k in 1:length(H)) {
      print(N[i])
      print(ikBT[j])
      print(H[k])
      runSim(ikBT[j], J, H[k], N[i], nstep, 1) # metropolis
      runSim(ikBT[j], J, H[k], N[i], nstep, 2) # glauber
     }
   }
}


# This is an example analysis
#
# magnetisationMetric   <- magnetisationMetricTM1D(1.0, 1.0, 1.0, 50, 500000, 1) # metropolis
# metricTM              <- magnetisationMetric[1]/magnetisationMetric # TM metric
# time                  <- 1:length(metricTM)
# Dt                    <- 1:length(metricTM)/metricTM # time evolution of the diffusion coefficient
# ff                    <- lm(metricTM ~ time)
# coeff                 <- ff$coefficients
# stdErrorsCoefficients <- coef(summary(ff))[, "Std. Error"] # metricTM == D_G t + intercept
# plot(time, metricTM)
# lines(time, ff$fitted.values)
