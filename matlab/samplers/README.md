# TT-IRT (Matlab version)
Inverse Rosenblatt Transports (IRT) + MCMC sampling using Tensor Train (TT) approximation

## Samplers

This directory contains functions for producing samples or, more generally, for transforming coordinates.
Each function contains an individual description in the header of the file.

 * `tt_irt_lin.m`     IRT computed using an approximate density in TT format and linear spline interpolation
 * `tt_irt_sqr.m`     IRT computed using an approximate square root of the density in TT format and linear spline interpolation
 * `tt_irt_fourier.m` IRT computed using an approximate square root of the density in TT format and Fourier basis interpolation
 * `mcmc_prune.m`     Metropolis-Hastings rejection algorithm (assumes independent proposals)
 * `iw_prune.m`       Importance Weighting correction method
 * `qmcnodes.m`       Produces a QMC lattice on [0,1]^d
 * `tt_dirt_sample.m` Invokes the Deep Inverse Rosenblatt Transform (DIRT) of given points on the reference domain
 * `randref.m`        Produces independent random samples distributed according to a given reference measure
 * `essinv.m`         Computes *N*/*ESS*, the reciprocal Effective Sample Size for given independent samples of exact and approximate log-density functions. The Chi^2-divergence between these densities can be estimated as *N*/*ESS* - 1.
 * `tt_irt_debias.m`  **(deprecated)** a driver function for switching between MCMC and Importance Weighting (IW) correction of samples produced by `tt_irt_lin`


## The *IRT* structure

The function `tt_dirt_approx` returns a **struct** that contains TT cores of all ratio functions together with necessary auxiliary information.
This structure can then be fed to `tt_dirt_sample` to map new samples through the DIRT, or to `tt_dirt_approx` again to add new levels to the existing structure. The *IRT* structure contains the following fields.

 * *x0* A cell array of grid points used in each variable at Level 0
 * *beta* A vector of bridging parameters (e.g. tempering powers). Must be in increasing order.
 * *reference* A string specifying the reference measure
 * *crossmethod* A string specifying the TT-Cross implementation used
 * *interpolation* A string specifying the basis for interpolation and CDF inversion in IRT
 * *evalcnt* A vector of numbers of function evaluations in TT-Cross from all levels
 * *F0* A cell array of TT cores of the density function at Level 0
 * *Fprev* A tt_tensor of the density/ratio function from the last level
 * *x* A vector of grid points in the reference variables stacked together
 * *F* A cell array of TT formats of ratio functions from all levels above 0. Each TT format is in turn a cell array similarly to *F0*.

See also the header of `tt_dirt_approx.m`.


## External samplers

Running scripts `check_mcmc`, `check_qmc` in the `utils/` directory, as well as `check_svn` in the `examples/predator_prey` directory, downloads the corresponding external packages, which are stored in the `samplers/` directory.
See the main [README](https://github.com/dolgov/TT-IRT/blob/master/README.md) for references.

 * `dramcode` Delayed Rejection Adaptive Metropolis (DRAM) code (for DRAM experiments)
 * `lattice-39102-1024-1048576.3600.txt` A generating vector for the QMC lattice produced by `qmcnodes`
 * `UWerr.m` Integrated autocorrelation time (IACT) estimator as described in [[Ulli Wolff, *Monte Carlo errors with less errors*, CPC 156(2)](https://doi.org/10.1016/S0010-4655(03)00467-3)]
