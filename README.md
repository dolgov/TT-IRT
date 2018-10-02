# TT-IRT
Inverse Rosenblatt Transform (Conditional Distribution) + MCMC sampling using Tensor Train approximation

Algorithms and running codes for the paper "Approximation and sampling of multivariate probability distributions in the tensor train decomposition" [[arXiv:1809.xxxx](http://arxiv.org/abs/)] by Sergey Dolgov, Karim Anaya-Izquierdo, Colin Fox and Robert Scheichl.

The toolbox consists of **matlab** and **python** parts. Moreover, the conditional distribution sampler (TT-CD) is also implemented in C with two types of integers, `matlab/tt_irt1_int64.c` and `python/tt_irt_py/tt_irt1_int32.c` which can be linked into Matlab MEX, Python CTypes or other project.

### Installation

The package is based on several projects: Tensor Train (TT) algorithms are implemented in [TT-Toolbox](https://github.com/oseledets/TT-Toolbox) (Matlab) and [ttpy](https://github.com/oseledets/ttpy) (python), the QMC lattice is produced using a generating vector from [F. Kuo website](http://web.maths.unsw.edu.au/~fkuo/lattice/index.html), and the integrated autocorrelation time (IACT) is estimated by [UWerr.m](https://www.physik.hu-berlin.de/de/com/ALPHAsoft) / [puwr.py](https://github.com/dhesse/py-uwerr).
The TT algorithms can be compared to [Delayed Rejection Adaptive Metropolis (DRAM)](http://helios.fmi.fi/~lainema/dram/) MCMC.
Installation steps are thus different in **matlab** and **python**.

**matlab/**

After changing to `matlab/` directory, the general configuration script `install` can be run in order to compile MEX files.
Moreover, this script will ask you whether the parallel *for* loops should be used for different runs of the tests.

External packages are downloaded automatically by `check_tt`, `check_mcmc` and `check_qmc` scripts when the corresponding test is run. Should the automated download fail, the scripts will tell you where you can get the necessary prerequisite manually.
Subdirectories of the external codes will be added to Matlab path.
However, the current `matlab/` directory contains all TT-IRT files and hence you don't need to have it in the path.


**python/**

First, install [ttpy](https://github.com/oseledets/ttpy) following any of the recipes listed at [https://github.com/oseledets/ttpy](https://github.com/oseledets/ttpy) (using *conda*, *pip* or directly the setup script).

The TT-CD sampler package `tt_irt` is located in `python/tt_irt_py` directory. Change there and run
```
python setup.py install [--user]
```
The installation can be checked by `verify.py` script. It will also try to download `puwr.py` for the IACT estimation.


### Scripts for running examples

All files for running experiments start with a `test_` prefix. Currently implemented are the shock absorber and the diffusion equation examples.

**matlab/**
 * `test_shock_absorber_tt.m`       Shock absorber test with TT-MH and TT-qIW methods
 * `test_shock_absorber_dram.m`     Shock absorber test with DRAM
 * `test_diffusion_tt.m`            Inverse diffusion test with TT-MH/TT-qIW
 * `test_diffusion_dram.m`          Inverse diffusion test with DRAM
 * `test_diffusion_qmcrat.m`        Inverse diffusion test with QMC ratio estimator

Each test can be run without any arguments. In this case, they will be interactively asked from the user. For batch runs, parameters can be passed in pairs of inputs ``'param_name1'``, ``param_value1``, ``'param_name2'``, ``param_value2``, and so on. For example,
```
test_shock_absorber_tt('D', 6, 'x', [], 'log2N', 14, 'delta', 0.05, 'n', 16, 'runs', 8)
```
will run the shock absorber test with all default parameters. Only a subset of arguments can be given, e.g.
```
test_shock_absorber_tt('log2N', 18, 'runs', 32)
```
will take the corresponding values from the inputs, and ask the user only for the remaining parameters.

Each test will print some statistical data (expected values, standard deviations, CPU times and so on).
Moreover, it will create all the variables in the main Matlab workspace, so they can be accessed afterwards.


**python/**
 * `test_shock_absorber_tt.py`      Shock absorber test with TT-MH

The parameters are set up statically in this file. In order to change them, go to lines 88-92 and edit
```
d = 6       # number of covariates
n = 16      # grid size
tol = 5e-2  # TT stopping threshold
log2N = 16  # log2(number of samples)
runs = 8    # number of runs
```


### Function files

Each function file contains an extended description. See e.g.
```
help('tt_irt')
```
or open the file in the editor.

**matlab/**
  * *Core*
    - `tt_irt.m`  TT-CD sampler from an approximate density in TT format
    - `tt_irt_debias.m`  MCMC and Importance Weighting (IW) correction of the TT-CD samples
    - `mcmc_prune.m`     Metropolis-Hastings rejection algorithm (assumes independent proposals)
    - `iw_prune.m`       Importance Weighting correction method
    - `amen_cross_s.m`   Enhanced TT-Cross algorithm for the TT approximation
    - `lagrange_interpolant.m`  Lagrange interpolation between two sets of points
    - `tt_sample_lagr.m`  Lagrange interpolation of a TT decomposition
    - `tracemultm.m`     Auxiliary file for conditioning of TT blocks in TT-CD (+ faster MEX version `tracemult.c`)
    - `tt_irt_mex.c`     Faster MEX version of the TT-CD sampler (uses `tt_irt1_int64.c`)
    - `qmcnodes.m`       Produces a QMC lattice on [0,1]^d
  * *Shock absorber example*
    - `parse_shock_inputs.m`    A function for requesting/extracting model parameters
    - `shock_log_prior.m`       Log-prior (s-Normal-Gamma) function
    - `shock_log_weibull.m`     Log-likelihood (censored Weibull) function
    - `shock_quantiles.m`       A function for computing posterior quantiles
    - `shock-xdata-d6.dat`      6-dimensional covariates data used in the paper (for reproducing the experiments)
  * *Inverse diffusion example*
    - `parse_diffusion_inputs.m`    A function for requesting/extracting model parameters
    - `als_cross_parametric.m`      ALS-Cross algorithm for the parametrized PDE solution, see [http://people.bath.ac.uk/sd901/als-cross-algorithm/](http://people.bath.ac.uk/sd901/als-cross-algorithm/)
    - `build_grid_and_kle.m`        Discretization of the diffusion equation
    - `diffusion_assem_solve.m`     Solution of the diffusion equation (for ALS-Cross)
    - `diffusion_likelihood.m`      Computation of the likelihood (solves the PDE for each sample)
    - `project_blockdiag_mex.c`     Faster MEX version of the reduction step in ALS-Cross.
    - `solve_blockdiag_mex.c`       Faster MEX version of the local solver in ALS-Cross.
  * *Utilities*
    - `check_mcmc.m`   Check/download/add-to-path for DRAM and UWerr
    - `check_tt.m`     Check/download/add-to-path for TT-Toolbox
    - `check_qmc.m`    Check/download for the QMC lattice generating vector
    - `install.m`      Compile MEX files, switch between sequential and parallel for loops.

**python/**
  * Core package in `tt_irt_py` directory
    - `tt_irt1_int32.c`  TT-CD sampler from an approximate density in TT format (uses 32-bit integers)
    - `tt_irt.py`        CTypes wrapper
    - `setup.py`         Compilation and installation script
  * *Shock absorber example*
    - `test_shock_absorber_tt.py`  Run script, encapsulates prior and likelihood functions
  * *Utilities*
    - `verify.py`        Check availability of ttpy and tt_irt, download puwr.py

