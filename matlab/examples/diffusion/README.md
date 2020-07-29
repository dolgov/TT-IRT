# TT-IRT (Matlab version)
Inverse Rosenblatt Transports (IRT) + MCMC sampling using Tensor Train (TT) approximation


## Inverse Diffusion equation

This example benchmarks the "stylized" diffusion equation introduced first in [[Eigel, M., Pfeffer, M. & Schneider, R. Adaptive stochastic Galerkin FEM with hierarchical tensor representations. Numer. Math. 136, 765–803 (2017)](https://doi.org/10.1007/s00211-016-0850-x)].
The "stylized" stands for a synthetic Fourier expansion of the diffusion coefficient with a tunable decay rate.
This directory implements an example of the Bayesian inference of the random variables given averaged noisy observations of the solution, as described in [[Dolgov, S., Anaya-Izquierdo, K., Fox, C. et al. Approximation and sampling of multivariate probability distributions in the tensor train decomposition. Stat Comput 30, 603–625 (2020)](https://doi.org/10.1007/s11222-019-09910-z)].

### Running the experiments

The files `test_diffusion_*.m` run a benchmark of the corresponding algorithms (see the main [README](https://github.com/dolgov/TT-IRT/blob/master/README.md)).
These will ask interactively (or take from *varargin*) the following parameters.
#### PDE/Prior parameters (all tests)
 * *sigma* Variance of the affine diffusion coefficient expansion (or of **log** of the coefficient for log-uniform and log-normal fields)
 * *corr_length* Correlation length of the expansion
 * *nu* Decay rate
 * *meshlevel* Level of the spatial discretization (where *meshlevel+1* corresponds to halving the step size in each variable)

#### Inverse problem parameters (all tests)
 * *sigma_n* Variance of the observation noise
 * *m0* Number of observation points in each variable (total number of observations is m0*m0)
 * *y0* Synthetic truth value of the random variables
 * *log2N* Base-2 logarithm of the number of samples produced in QMC or MCMC
 * *runs* Number of runs (replicas) of the experiment

#### Single-level TT approximation parameters (`test_diffusion_tt.m`)

 * *ny*  Number of discretization points in each random variable in the forward PDE model
 * *rmax* Maximum TT rank in the forward problem
 * *npi* Number of discretization points in each random variable in the posterior density function
 * *delta* Truncation threshold in TT-Cross
 * *correction* Debiasing algorithm (MCMC or Importance Weighting)

#### DIRT approximation parameters (`test_diffusion_dirt.m`)

 * *ny*  Number of discretization points in each random variable in the forward PDE model
 * *rmax* Maximum TT rank in the forward problem
 * *npi* Number of discretization points in each random variable in the ratio functions on each level of DIRT
 * *rpi* TT rank of each DIRT ratio function (fixed-rank TT decompositions are used)
 * *beta* A vector of tempering powers. Should be in increasing order, and the last value should be 1.

Only MCMC debiasing is benchmarked with DIRT for simplicity.

#### DRAM benchmark (`test_diffusion_dram.m`)

Note that the burn-in of 6000 samples is discarded, so *log2N* must be greater than 12.


### Output variables

Each script creates certain variables that remain accessible in the main Matlab workspace.
Use the `who` command for the full list.
Some of the interesting variables are:

 * *ttimes_forward* A vector of CPU times of solving the forward model in each run
 * *nsolves_forward* Numbers of deterministic PDE solves in constructing the forward model surrogate
 * *ttimes_pi* CPU times of approximating the posterior density
 * *ttimes_invcdf* CPU times of (single-level) IRT
 * *ttimes_dirt* CPU times of DIRT
 * *ttimes_debias* CPU times of computing exact density values in MCMC
 * *bias* Number of rejections or bias in MCMC/Importance weighting
 * *tau_tt* Integrated Autocorrelation Times (IACT) of MCMC chain
 * *err_Pi* TT approximation (excluding discretization) error
 * *Q_tt* A *runs* x *2* matrix of quantities of interest. Each row contains *[mean_flux, exceedance_probability]*.
 * *ttimes_dram* CPU times of DRAM
 * *num_of_rejects* Number of rejections in DRAM
 * *tau_dram* IACT of DRAM chain (without burn-in)
 * *Q_dram* A *runs* x *2* matrix of quantities of interest (computed without burn-in)
 * *ttimes_qmcrat* CPU times of computing exact density values at QMC points
 * *Q_qmcrat* A *runs* x *2* matrix of quantities of interest

### Function files

 * `parse_diffusion_inputs.m`    A function for requesting/extracting model parameters
 * `build_grid_and_kle.m`        Discretization of the diffusion equation
 * `diffusion_assem_solve.m`     Solution of the diffusion equation (for ALS-Cross)
 * `diffusion_likelihood.m`      Computation of the likelihood (solves the PDE for each sample)

### External file used from `constructors/`
 * `als_cross_parametric.m`      ALS-Cross algorithm for the parametrized PDE solution, see also [http://people.bath.ac.uk/sd901/als-cross-algorithm/](http://people.bath.ac.uk/sd901/als-cross-algorithm/)
`

### External files used from `utils/`

 * `project_blockdiag_mex.c`     Faster MEX version of the reduction step in ALS-Cross.
 * `solve_blockdiag_mex.c`       Faster MEX version of the local solver in ALS-Cross.
