# isingLenzMC: Monte Carlo for Classical Ising Model

[![CRAN Version](https://www.r-pkg.org/badges/version/isingLenzMC)](https://cran.r-project.org/package=isingLenzMC)
[![Total RStudio Cloud Downloads](https://cranlogs.r-pkg.org/badges/grand-total/isingLenzMC?color=brightgreen)](https://cran.r-project.org/package=isingLenzMC)
[![RStudio Cloud Downloads](https://cranlogs.r-pkg.org/badges/isingLenzMC?color=brightgreen)](https://cran.r-project.org/package=isingLenzMC)
[![arXiv:1606.08693](http://img.shields.io/badge/arXiv-1606.08693-B31B1B.svg)](https://arxiv.org/abs/1606.08693)
[![arXiv:1405.4497](http://img.shields.io/badge/arXiv-1405.4497-B31B1B.svg)](https://arxiv.org/abs/1405.4497)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1065942.svg)](https://doi.org/10.5281/zenodo.1065942)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17151290.svg)](https://doi.org/10.5281/zenodo.17151290)

* [Stable release on CRAN](https://CRAN.R-project.org/package=isingLenzMC)
* [Development repository](https://github.com/msuzen/isingLenzMC)

## Description

Classical Ising Model is a land mark system in statistical physics. The model explains
the physics of spin glasses and magnetic materials, and cooperative phenomenon
in general, for example phase transitions and neural networks. This package provides
utilities to simulate one dimensional Ising Model with Metropolis and Glauber Monte
Carlo with single flip dynamics in periodic boundary conditions. Utility functions
for exact solutions are provided. Such as transfer matrix for 1D. Example use cases
are as follows: Measuring effective ergodicity and power-laws in so called
functional-diffusion.

## Scientific Output

The package is used in the following peer-reviewed works:

* Effective ergodicity in single-spin-flip dynamics  
  Mehmet Suezen, [Phys. Rev. E 90, 032141](https://doi.org/10.1103/PhysRevE.90.032141)  
  [Dataset](https://doi.org/10.5281/zenodo.1065942) (2014)
* Anomalous diffusion in convergence to effective ergodicity,
  Suzen, Mehmet, [Physica Scripta](https://doi.org/10.1088/1402-4896/ae483e) (2026)
  [Dataset](https://doi.org/10.5281/zenodo.18223717)

See Bibtex entries at [CITATION](inst/CITATION).

We provide replicating results from the papers under examples directory.

* [Measuring effective ergodicity on differing temperature ranges](inst/examples/effectiveErgodicity/README.md)
* [Ergodic Dynamics of Ising Model : functional-diffusion regimes](inst/examples/powerLawErgodicity)
  * There are quite a lot of functionality here that are not yet in the main package:
    * `log-log` regression for time-power laws, with KS-test windowing.
    * Usage of distributional power-laws.
    * Parallel run utilities.
    * Finite Size Scaling (FSS) analysis.
    * Autocorrelation function.
    * Bias-corrected bootstrapping for error bars.
    * Enhance plotting utilities: Publication quality.

## Dev Scratch Notes

* Building and checks for the package release before CRAN upload.

```bash
R CMD build --compact-vignettes=both isingLenzMC
R CMD check --as-cran  isingLenzMC_<version>.tar.gz
```

* Generating `init.c`: This is a note and not needed in general. This is due to fact that we have core C library.

```R
tools::package_native_routine_registration_skeleton("isingLenzMC", character_only = FALSE)
```

## License

This project and all contributions are licensed under :

* All non-code  [![License: CC BY 4.0](https://i.creativecommons.org/l/by/4.0/88x31.png)](https://creativecommons.org/licenses/by/4.0/)
* Code under [![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
