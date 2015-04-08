### R code from vignette source 'isingLenzMC.Rnw'

###################################################
### code chunk number 1: conf
###################################################
require(isingLenzMC)
set.seed(123456)
N <- 7 
myInitialConfig <- genConfig1D(N)
myInitialConfig 
myNextConfig    <- flipConfig1D(myInitialConfig)
myNextConfig   
# nearest neighbour energy for initial config
lattice1DenergyNN(myInitialConfig) 
# transition probability  at J=H=1/kBT=1.0
transitionProbability1D(1.0, myInitialConfig, myNextConfig, 1.0, 1.0, 1) # Metropolis


###################################################
### code chunk number 2: conf
###################################################
require(isingLenzMC)
set.seed(123456)
N <- 7 
myInitialConfig <- genConfig1D(N)
myInitialConfig
# 1 step Monte Carlo move
isStep1D(1.0, myInitialConfig, 1.0, 1.0, 1) # Metropolis


###################################################
### code chunk number 3: conf
###################################################
Tm <- transferMatrix(1.0, 1.0, 1.0)
# Free Energy
log(Tm$evalues[1]^7 + Tm$evalues[2]^7)


###################################################
### code chunk number 4: conf (eval = FALSE)
###################################################
## require(isingLenzMC)
## set.seed(123456)
## ensembleM <- 0.9934346
## N         <- 200
## x         <- genConfig1D(N)
## mcData    <- isPerform1D(1.0, x, 1.0, 1.0, 10000, ensembleM, 1)  # Metropolis


