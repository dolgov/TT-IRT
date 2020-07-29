# TT-IRT (Matlab version)
Inverse Rosenblatt Transports (IRT) + MCMC sampling using Tensor Train (TT) approximation


## Lorenz-d example

This is a variable-dimensional (Lorenz-40 is the canonical setup) example of inferring and quantifying uncertainty of the initial state of a chaotic dynamical system, given noisy observations of a fraction of the final state, such as even coordinates only. In addition to the nearest-neighbour coupling of the variables, there is a cyclic interaction between the first and the last variables, which makes the posterior density more difficult for straightforward tensor approximations.

### Running the experiments

The file `test_lorenz.m` runs a DIRT benchmark (see also the main [README](https://github.com/dolgov/TT-IRT/blob/master/README.md)).
It will ask interactively (or take from *varargin*) the following parameters.

#### Problem parameters
 * *d* Problem dimension. Using *d=40* reproduces the experiment in the paper, although it may take 5-10 minutes. A reduced example with *d=10* is recommended for hands-on trials.
 * *sigma_n* Standard deviation of the observation noise (must be >0)
 * *x0true* Synthetic true initial state (vector of *d* values)
 * *sigma_truth* Standard deviation of the initial state perturbationbs_times
 * *Nsamples* Length of MCMC chain to produce

#### DIRT approximation parameters
 * *n* Number of grid points in each variable
 * *a* Size of the domain in the original space (must be >0)
 * *R0* TT rank for the ratio functions at each DIRT layer (all TT ranks are set to *R0*)
 * *beta* A vector of tempering powers. Must be in increasing order, and the last element must be 1.


### Output variables

The script creates certain variables that remain accessible in the main Matlab workspace.
Use the `who` command for the full list.
Some of the interesting variables are:

 * *tau_iact* Integrated Autocorrelation Time (IACT) of MCMC
 * *tau_ess* N/ESS, the reciprocal Effective Sample Size
 * *IRT* DIRT structure (contains all TT decompositions)
 * *z* DIRT samples (in the form of a *Nsamples* x *d* matrix) (**before** rejection)
 * *lFapp* log(proposal density) values (**before** rejection)
 * *lFex* log(exact density) values (**before** rejection)



### Function files

 * `lorenz_ll.m` A log-likelihood function
 * `lorenz_rhs.m` A function returning the right-hand side of the forward model

