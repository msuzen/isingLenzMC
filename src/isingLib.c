/* 
  MC Dynamics of Ising Model
    C functions
  (c) 2013 by msuzen (Dr.Mehmet Suezen)
  GPLv3 or higher
*/
#include <math.h>
#include <stdio.h>
#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>
#include <Rinternals.h>

/*
  Generate Random Configuration in 1D
    Author : msuzen
    Arguments  
    ll     : number of sites 
    myNum  : site spin configuration  [1 or -1]
     Uses R's unif_rand(); (default generator MT if not changed from called)
 */
void genConfig1D(int *ll, double *myNum) {
  double uNum=0;
  int i, n=ll[0];
  for(i=0; i<n; i++) {
    GetRNGstate();
    myNum[i] = unif_rand();
    PutRNGstate();
    if(myNum[i] >= 0.5) {
      uNum=1.0; 
     } else {
      uNum=-1.0; 
    }
    myNum[i]=uNum;
  }
}

/*
    Random single given flip configuration in 1D
    Author : msuzen
    Arguments  
    ll     : number of sites 
    myNum : site spin configuration  [1 or -1]
     Uses R's unif_rand() (default generator MT)
 */
void flipConfig1D(int *ll, double *myNum) {
  int n=ll[0], rn;
  GetRNGstate();
  rn= (int) ((double) n* unif_rand());
  PutRNGstate();
  myNum[rn] = -1.0 * myNum[rn];
}

/*
    Random flip given configuration in 1D many times
    Author : msuzen
    Arguments  
    ll     : number of sites 
    myNum  : site spin configuration  [1 or -1]
    upperF : Number of flips to perform 
     Uses R's unif_rand() (default generator MT)
 */
void flipConfig1Dmany(int *ll, double *myNum, int *upperF) {
  int i=0;
  for(i=0;i< upperF[0];i++) {
    flipConfig1D(ll, myNum);
  } 
}

/*
  Sum given vector
    Author : msuzen
    Arguments
    ll     : number of sites 
    vec    : site spin configuration 
    sum    : sum over vec to be returned
*/
void sumVec(int *ll, double *x, double *sum) {
  int i, n = ll[0];
  sum[0] = 0.0;
  for(i=0; i< n; i++) {
    sum[0] += x[i];
  }
}

/*
  Nearest-Neighbour energy in periodic boundary conditions in 1D
    Author : msuzen
    Arguments  
    ll     : number of sites 
    vec    : site spin configuration 
    energy : energy to be returned 
 */
void lattice1DenergyNN(int *ll, double *vec, double *energy) {
  int i, n = ll[0];
  energy[0] = 0.0;
  for(i=1; i< n; i++) {
    energy[0] += vec[i] * vec[i-1];
  }
  /* Additional contribution because of
     periodic boundary conditions */
  energy[0] += vec[(n-1)] * vec[0];
}

/* Total Energy in 1D
    Author : msuzen
   Arguments
   ll   : length of the configuration
   vec  : spin configuration
   J    : interaction strength
   H    : external field
   totalEnergy
 */
void totalEnergy1D(int *ll, double *vec, double *J, double *H, double *totalEnergy) {
  double energyNN[1], energySum[1];
  lattice1DenergyNN(ll, vec, energyNN);
  sumVec(ll, vec, energySum);
  totalEnergy[0] = J[0] * energyNN[0] + H[0] * energySum[0];
}

/* 
  Compute Transition probability
   Author : msuzen
   Arguments
   ikBT    : 1 over Temperature times Boltzmann constant
   ll      : n spins
   vecOrg  : Original spin states
   vecFlip : one flip spin states
   prob    : transition probability
   probSel : which transition probability to use 1 Metropolis 2 Glauber
*/       
void transitionProbability1D(double *ikBT, int *ll, double *vecOrg, double *vecFlip, double *J, double *H, double *prob, int *probSel) {
  double energyOrg[1], energyFlip[1], DeltaE;
  totalEnergy1D(ll, vecOrg, J, H, energyOrg);
  totalEnergy1D(ll, vecFlip, J, H, energyFlip);
  DeltaE  = energyOrg[0] - energyFlip[0];
  if(probSel[0] == 1) prob[0]   = fmin(1, exp(-ikBT[0]*DeltaE));  /* Metropolis */
  if(probSel[0] == 2) prob[0]   = 1.0/(1.0+exp(ikBT[0]*DeltaE)); /* Glauber    */
}

/*
  Take One importance sampling MC Step 
   Author : msuzen
   Arguments
   ikBT    : 1 over Temperature times Boltzmann constant
   ll      : n spins
   vec     : Original spin states
   prob    : transition probability
   accept  : if step accepted it will be 1 (so pass/make it zero before running this function)
   probSel : which transition probability to use 1 Metropolis 2 Glauber
*/
void isStep1D(double *ikBT, int *ll, double *vec, double *J, double *H, double *prob, int *accept, int *probSel) {
  double energyOrg[1], energyFlip[1], DeltaE, flipId, rnd;
  int rn,i;
  totalEnergy1D(ll, vec, J, H, energyOrg);
  GetRNGstate();
  rn= (int) ((double) ll[0]*unif_rand());
  PutRNGstate();
  //printf("flip id=%d\n", rn);
  vec[rn] = -1.0*vec[rn]; // flip 
  totalEnergy1D(ll, vec, J, H, energyFlip);
  vec[rn]   = -1.0*vec[rn]; // flip back
  DeltaE    = energyOrg[0] - energyFlip[0];
  if(probSel[0] == 1) prob[0]   = fmin(1, exp(-ikBT[0]*DeltaE));          /* Metropolis */
  if(probSel[0] == 2) prob[0]   = 1.0/(1.0+exp(ikBT[0]*DeltaE)); /* Glauber */
  GetRNGstate();
  rnd       = unif_rand();
  PutRNGstate();
  accept[0] = 0;
  if(prob[0] > rnd) {
    accept[0] = 1;
    vec[rn]   = -1.0*vec[rn]; //flip it accepted
  }
}

/*

#  Perform importance sampling MC on 1D
#
#   Author : msuzen
#   Arguments
#   ikBT    : 1 over Temperature times Boltzmann constant
#   ll      : number of sites 
#   vec     : 1D Spin sites on the lattice
#   J       : interaction strength
#   H       : external field
#   ensembleM : Value of the theoretical magnetization (could be thermodynamic limit value)
#   omegaM    : Fluctuating metric vector for Magnetisation (nstep length)
#   nstep     : number of MC steps requested
#   naccept   : number of MC steps accepted
#   nreject   : number of MC steps rejected
#   times     : this is 
#   probSel : which transition probability to use 1 Metropolis 2 Glauber
#
*/
void isPerform1D(double *ikBT, int *ll, double *vec, double *J, double *H, double *ensembleM, 
                 double *omegaM, int *nstep, int *naccept, int *nreject, int *times, int *probSel) {
  int i, k, accept[1];
  double prob[1];
  double timeM[1], diff, magtime;  /* time average magnetization */
  prob[0]  = 0.0;
  timeM[0] = 0.0;
  magtime  = 0.0;
   k = 0;
  for(i=0 ;i < nstep[0]; i++) {
    accept[0] = 0;
    isStep1D(ikBT, ll, vec, J, H, prob, accept, probSel);
    if(accept[0] < 1) nreject[0]++;
    if(accept[0] > 0) { 
      times[naccept[0]] = i;
      sumVec(ll, vec, timeM);
      magtime += timeM[0]/ll[0];
      diff    = (magtime/(k+1) - ensembleM[0]); 
      omegaM[k] = diff*diff;
      k++;
      naccept[0]++;
    }
   }
} 
