# TT-IRT (Matlab version)
Inverse Rosenblatt Transports (IRT) + MCMC sampling using Tensor Train (TT) approximation


## Predator-prey example

This is an 8-dimensional example of calibrating a differential equation model by observing noisy states of the trajectory.
Some of the distant variables are strongly nonlinearly correlated, which makes the posterior density difficult for a straightforward approximation.

### Running the experiments

The files `test_predator_prey_*.m` run a benchmark of the corresponding algorithms (see the main [README](https://github.com/dolgov/TT-IRT/blob/master/README.md)).
These will ask interactively (or take from *varargin*) the following parameters.

#### Model parameters (all tests)
 * *sigma_n* Observation noise variance
 * *xtrue* A vector of synthetic true values of the parameters
 * *obs_times* A vector of time points where the trajectory is observed
 * *data* A matrix of observed states. Each row contains [P(t_i), Q(t_i)] observed at *obs_times*(i). See `pp_observables.dat` for a matrix of values used in the [[DIRT](https://arxiv.org/abs/2007.06968)] paper.
 * *domain* The support interval of each variable, **relative** to *xtrue*
 * *Nsamples* Number of MCMC, DRAM or SVN samples
 * *runs* Number of runs (replicas) of the experiment (must be >1 for error estimation)

#### DIRT approximation parameters (`test_predator_prey_dirt.m`)
 * *n* Number of grid points in each variable
 * *R0* TT rank for the ratio functions at each DIRT layer (all TT ranks are set to *R0*)
 * *beta* A vector of tempering powers. Must be in increasing order, and the last element must be 1.

#### DRAM benchmark (`test_predator_prey_dram.m`)

All samples are used (no burn-in)

#### Stein Variational Newton benchmark (`test_predator_prey_svn.m`)
 * *stepsize* Step size (damp) of the Newton increment
 * *itermax* maximum number of Newton iterations
 * *initial_std* standard deviation of the initial cloud of points, **relative** to *xtrue*


### Output variables

Each script creates certain variables that remain accessible in the main Matlab workspace.
Use the `who` command for the full list.
Some of the interesting variables are:

 * *num_of_rejects* A vector of numbers of rejections in MCMC, produced in each simulation
 * *tau* Integrated Autocorrelation Time (IACT) of MCMC or DRAM
 * *tau_ess* N/ESS, the reciprocal Effective Sample Size
 * *evalcnt* Total number of function evaluations in DIRT
 * *ttimes_approx* CPU time of DIRT construction
 * *ttimes_sample* CPU time of sampling, MCMC only for DIRT, the entire algorithms for DRAM and SVN
 * *Fdist* Average Forstner-Moonen distance between posterior covariance matrices from different runs
 * *meanZ* Mean variables (1 x 8 x *runs* tensor)
 * *covZ*  Covariances (8 x 8 x *runs* tensor)

If the standard sequential *for* is used, additional variables created are:

 * *IRT* DIRT structure (contains all TT decompositions)
 * *z* MCMC/DRAM/SVN samples (in the form of a *Nsamples* x 8 matrix)
 * *lFapp* log(proposal density) values (after rejection)
 * *lFex* log(exact density) values (after rejection)



### Function files

 * `check_svn.m` A helper function to check/download Stein Variational Samplers repository
 * `forward_solve.m` A wrapper of the forward ODE model for SVN
 * `grad_mlpt.m` A wrapper for the gradient of the minus-log-posterior for SVN
 * `parse_pp_inputs.m` A helper function for parsing input parameters
 * `PP_loglikelihood.m` A log-likelihood function
 * `pp_observables.dat` Observables data used in the paper (for reproducing the experiments)
 * `PP_RHS_grad.m`  A function returning the right-hand side for the adjoint model
 * `PP_RHS.m` A function returning the right-hand side of the forward model

