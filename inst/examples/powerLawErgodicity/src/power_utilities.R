#
# Ergodic Dynamics of Ising Model : Diffusion regimes : Utility functions
#  (c) 2013, 2014, 2015, 2016, 2025 Süzen
#  GPL v3

require("isingLenzMC")
library(boot)

#
# List of utility methods
#
# * getEnsembleMagnetisation: Analytic solution of ensemble average.
# * isPerform1D_two: Calling C level importance sampling, returns the last config as well.
# * performIS: This is to perform 1D IsingLenz Importance Sampling
# * fit_power_law: Fits power-law given sim result list.
# * ks_distance: Compute KS-stat, D distance.
# * select_powerlaw: Select a power-law.
# * is1Dpowerlaw: Main utility to generate data/detailted plots given
#   N, beta, dynamics, and field.


getEnsembleMagnetisation <- function(n, ikBT, J, H) {
  # See IsingLenzMC's examples
  # Maxima output above
  if (n < 700) { # This is ad-hoc assumption is that we recover n -> inf solution above this
    cHB <- cosh(H * ikBT)
    sHB <- sinh(H * ikBT)
    e4B <- exp(-4 * ikBT * J)
    sq4B <- sqrt(e4B + sHB^2)
    isq4B <- ikBT * cHB * sHB / sq4B
    # Per site ensemble averaged magnetisation
    eMag <- isq4B + ikBT * sHB
    eMag <- n * eMag * (sq4B + cHB)^(n - 1)
    eMag <- eMag + n * (ikBT * sHB - isq4B) * (cHB - sq4B)^(n - 1)
    eMag <- eMag * ((sq4B + cHB)^(n) + (cHB - sq4B)^(n))^(-1)
    if (is.nan(eMag) || is.infinite(eMag)) {
      return(0.9934344)
    } # 0.9934344 infinite solution
    return(eMag / ikBT / n)
  } else {
    return(0.9934344)
  }
}


isPerform1D_two <- function(ikBT, x, J, H, nstep, ensembleM, probSel) {
  if (!is.numeric(x)) {
    stop("argument x must be numeric")
  }
  omegaM <- rep(0, nstep)
  times <- rep(1, nstep)
  naccept <- 0
  nreject <- 0
  out <- .C("isPerform1D",
    ikBT = as.double(ikBT), n = as.integer(length(x)),
    vec = as.double(x), J = as.double(J), H = as.double(H),
    ensembleM = as.double(ensembleM), omegaM = as.double(omegaM),
    nstep = as.integer(nstep), naccept = as.integer(naccept),
    nreject = as.integer(nreject), times = as.integer(times),
    probSel = as.integer(probSel)
  )
  rout <- list(
    omegaM = out$omegaM, nreject = out$nreject,
    naccept = out$naccept, times = out$times, config = out$vec
  )
  rout
}

average_magnetisation_traj <- function(N, ikBT, J, H, transP, max_mc_time, delta_mc_time, nrepeat) {
  # Average magnetisation trajectory
  ensembleM <- getEnsembleMagnetisation(N, ikBT, J, H) # ensemble average magnetisation per-site
  config <- genConfig1D(N) # Random initial condition

  nconfig <- as.integer(max_mc_time / delta_mc_time)
  mc_time <- seq(from = 1e3, by = delta_mc_time, length.out = nconfig)
  average_magnetisation_repeats <- vector("list", nrepeat)

  for (j in 1:nrepeat) {
    average_magnetisation <- vector("numeric", nconfig)
    for (i in 1:nconfig) {
      mySim <- isPerform1D_two(ikBT, config, J, H, delta_mc_time, ensembleM, transP)
      config <- mySim$config
      average_magnetisation[i] <- mean(config)
    }
    average_magnetisation_repeats[[j]] <- average_magnetisation
  }


  list(
    "N" = N, "ikBT" = ikBT, "J" = J, "H" = H, "transP" = transP, "max_mc_time" = max_mc_time,
    "delta_mc_time" = delta_mc_time, "mc_time" = mc_time, nrepeat = nrepeat,
    "average_magnetisation_repeats" = average_magnetisation_repeats
  )
}

performIS <- function(ikBT, J, H, N, nstep, transP) {
  #
  # This is to perform 1D IsingLenz Importance Sampling
  #
  # - Gamma : Rate of Thirumalai - Mountain Fluctuation Metric is returned.
  # - nstep is requested MC simulation time.
  #
  #
  # Gamma:
  #  isPerform1D C-code records only naccept values. In the remaining of the
  #  returned arraystimes`` and `omegaM` is just initialisations values
  #  This is due to random nature of accept times.
  #  Even though
  #
  # Ensemble average magnetisation per-site
  ensembleM <- getEnsembleMagnetisation(N, ikBT, J, H)
  x <- genConfig1D(N) # Random initial condition
  mcrun <- isPerform1D_two(ikBT, x, J, H, nstep, ensembleM, transP)
  range_of_accept <- 1:mcrun$naccept
  gamma <- mcrun$omegaM[range_of_accept]
  accept_times <- mcrun$times[range_of_accept]
  if (transP == 1) {
    dynamics <- "Metropolis"
  }
  if (transP == 2) {
    dynamics <- "Glauber"
  }
  list(
    "ikBT" = ikBT, "J" = J, "H" = H, "N" = N,
    "dynamics" = dynamics,
    "accept_times" = accept_times,
    "gamma" = gamma
  )
}


performIS_repeat <- function(ikBT, J, H, N, nstep, transP, nrepeat) {
  # Repeated perf
  ensembleM <- getEnsembleMagnetisation(N, ikBT, J, H)

  accept_times_repeats <- vector("list", nrepeat)
  gamma_repeats <- vector("list", nrepeat)
  for (i in 1:nrepeat) {
    x <- genConfig1D(N) # Random initial condition
    mcrun <- isPerform1D_two(ikBT, x, J, H, nstep, ensembleM, transP)
    range_of_accept <- 1:mcrun$naccept
    gamma_repeats[[i]] <- mcrun$omegaM[range_of_accept]
    accept_times_repeats[[i]] <- mcrun$times[range_of_accept]
  }

  if (transP == 1) {
    dynamics <- "Metropolis"
  }
  if (transP == 2) {
    dynamics <- "Glauber"
  }
  list(
    "ikBT" = ikBT, "J" = J, "H" = H, "N" = N,
    "dynamics" = dynamics,
    "accept_times_repeats" = accept_times_repeats,
    "gamma_repeats" = gamma_repeats
  )
}

self_time_correlation <- function(observed_series) {
    # Autocorrelation function 
    # C(t) = <A(t) A(t)>
    # C(t) = 1/N Sum_0^n A(0) A(t_n)
    observed_zero <- observed_series[1]
    mtime <- length(observed_series)
    c_t <- vector("numeric", mtime)
    c_t <- observed_zero %*% observed_series
    c_t <- cumsum(c_t) / (1:mtime)
    c_t
}

get_mean_boot95_bca <- function(data) {
    # Bias-corrected bootstraped means and 95$ CI.
    mean_boot <- function(data, indices) {
      return(mean(data[indices]))
    }
    boot_result <- boot(data = as.numeric(data), statistic = mean_boot, R = 1000)
    boot_means <- boot_result$t
    ci <- boot::boot.ci(boot_result, type = "bca")$bca[4:5]
    list("mean" = boot_result$t0, "lower" = ci[1], "upper" = ci[2], "ci" = ci)
  }




