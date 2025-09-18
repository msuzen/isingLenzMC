# Ergodic Dynamics of Ising Model : Diffusion regimes

    
   (c) 2013, 2014, 2015, 2016, 2025   
   Süzen   
   GPL v3 

Reproducing the data and analysis for the article: 

Anomalous diffusion in convergence to effective ergodicity
M. Suezen  
[arXiv:1606.08693](https://arxiv.org/abs/1606.08693)

## Notebooks 

We provide notebooks for data generation and analysis 

* `data_generate.ipynb` : Notebook will generate the dynamic trajectories of Ising Models with external field under 
different temperatures at N=512,1024 and Metropolis/Glauber dynamics. 
Parameters reads 
```R
N    <- c(512, 1024)
ikBT <- seq(0.5,2.0,0.025)
H    <- c(1.0)
#
J      <- 1.0
nstep  <- 2.0e6
```
* `data_analysis.ipynb` : This will compute the approach the ergodicity
and power-law exponents. Ploting the results. Resulting data.frames are also stored. 

## Generated Data Outputs

Expected output is `ising_ergo.RData` file. 
This contains all the trajectories. 

## Analysis Outputs

Analysis generates eps, png and csv files. 
Plots are for approach to ergodicity and 
power-law exponents showing

```txt
glauberN1024.csv
glauberN1024.eps
glauberN1024.png
glauberN512.csv
glauberN512.eps
glauberN512.png
glauber_N1024scaling_exponents.eps
glauber_N1024scaling_exponents.png
glauber_N512scaling_exponents.eps
glauber_N512scaling_exponents.png
metropolisN1024.csv
metropolisN1024.eps
metropolisN1024.png
metropolisN512.csv
metropolisN512.eps
metropolisN512.png
metropolis_N1024scaling_exponents.eps
metropolis_N1024scaling_exponents.png
metropolis_N512scaling_exponents.eps
metropolis_N512scaling_exponents.png
```

## License

This project and all contributions are licensed under :
* All non-code  [![License: CC BY 4.0](https://i.creativecommons.org/l/by/4.0/88x31.png)](https://creativecommons.org/licenses/by/4.0/)
* Code under [![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)