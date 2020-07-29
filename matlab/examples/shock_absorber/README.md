# TT-IRT (Matlab version)
Inverse Rosenblatt Transports (IRT) + MCMC sampling using Tensor Train (TT) approximation


## Shock absorber example

This is an example of inferring parameters of a Weibull distribution of the time-to-fail from observed times, which are partially censored.
This benchmark is particularly convenient since the unnormalised posterior density can be computed explicitly at any value of the parameters.
It is also naturally extensible in the dimension by changing the number of covariates.

### Running the experiments

The files `test_shock_absorber_*.m` run a benchmark of the corresponding algorithms (see the main [README](https://github.com/dolgov/TT-IRT/blob/master/README.md)).
These will ask interactively (or take from *varargin*) the following parameters.

#### Model parameters (all tests)
 * *D* Number of extra covariates (in addition to the two parameters theta1, theta2 of the Weibull density)
 * *x* A *D* x 38 array of covariate values. The file `shock-xdata-d6.dat` contains the values used in the [[Statistics and Computing](https://doi.org/10.1007/s11222-019-09910-z)] paper.
 * *log2N* Base-2 logarithm of the number of samples *N* produced in MCMC and QMC
 * *runs* Number of runs (replicas) of the experiment

#### Single-level TT approximation parameters (`test_shock_absorber_tt.m`)

 * *n* Number of discretization points in each random variable in the posterior density function
 * *delta* Truncation threshold in TT-Cross

#### DRAM benchmark (`test_shock_absorber_dram.m`)

A burn-in of *N*/4 samples is discarded, so *log2N* should be greater than 2.



### Output variables

Each script creates certain variables that remain accessible in the main Matlab workspace.
Use the `who` command for the full list.
Some of the interesting variables are:

 * *N_cross* A vector of numbers of density evaluations in TT-Cross in each simulation
 * *ttimes_cross* CPU times of TT-Cross
 * *ttimes_invcdf* CPU times of the inverse Rosenblatt transform of MC points
 * *ttimes_dram* CPU times of DRAM
 * *num_of_rejects* Numbers of rejections in MCMC or DRAM
 * *tauint_tt* Integrated Autocorrelation Times (IACT) after MCMC
 * *tauint_dram* IACT of DRAM
 * *err_TT* Estimated TT approximation (excluding discretization) error. Needs *runs*>=8.
 * *Q_mh* A *runs* x 2 array of Quantities of Interest computed using MCMC. Each row contains [mean_quantile, quantile_of_mean], where mean_quantile is the posterior-mean 95%-quantile of the Weibull density, whereas quantile_of_mean is the 95%-quantile of the posterior-mean density.
 * *Q_iw* A *runs* x 2 array of Quantities of Interest computed using Importance Weighting seeded with QMC quadrature.
 * *Q_dram* A *runs* x 2 array of Quantities of Interest computed using DRAM


### Function files

 * `parse_shock_inputs.m`    A function for requesting/extracting model parameters
 * `shock_log_prior.m`       Log-prior (s-Normal-Gamma) function
 * `shock_log_weibull.m`     Log-likelihood (censored Weibull) function
 * `shock_quantiles.m`       A function for computing posterior quantiles
 * `shock-xdata-d6.dat`      6-dimensional covariates data used in the paper (for reproducing the experiments)
