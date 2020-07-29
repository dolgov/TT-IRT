# TT-IRT
Inverse Rosenblatt Transports (IRT) + MCMC sampling using Tensor Train (TT) approximation

Algorithms and examples are described in "Approximation and sampling of multivariate probability distributions in the tensor train decomposition" [[Statistics and Computing](https://doi.org/10.1007/s11222-019-09910-z)] and "Deep Composition of Tensor Trains using Squared Inverse Rosenblatt Transports" [[arXiv:2007.06968](https://arxiv.org/abs/2007.06968)].

The toolbox consists of **matlab** and **python** parts. Moreover, the inverse Rosenblatt transform based on linear splines is also implemented in C with two types of integers, `matlab/utils/tt_irt1_int64.c` and `python/tt_irt_py/tt_irt1_int32.c` which can be linked into Matlab MEX, Python CTypes or other project.

### Changelog

 * `2.0` Squared (SIRT) and Deep (DIRT) Inverse Rosenblatt Transports [[arXiv:2007.06968](https://arxiv.org/abs/2007.06968)]. Functions are structured in directories. New `predator_prey` and `lorenz` examples for DIRT.
 * `1.1` improved baseline version where logarithms of density functions are used to prevent overflows. Debiasing functions return more diagnostic information. Grid work is simplified such that the sets of grid points are the same for `tt_irt` and tt_tensors. New interface of `tracemult*`, MEX version can handle complex numbers.
 * `1.0` original version used for numerical experiments in the [[Statistics and Computing](https://doi.org/10.1007/s11222-019-09910-z)] paper.

### Installation

The package is based on several projects: Tensor Train (TT) algorithms are implemented in [TT-Toolbox](https://github.com/oseledets/TT-Toolbox) (Matlab) and [ttpy](https://github.com/oseledets/ttpy) (python), the QMC lattice is produced using a generating vector from [F. Kuo website](http://web.maths.unsw.edu.au/~fkuo/lattice/index.html), and the integrated autocorrelation time (IACT) is estimated by [UWerr.m](https://www.physik.hu-berlin.de/de/com/ALPHAsoft) / [puwr.py](https://github.com/dhesse/py-uwerr).
The TT algorithms can be compared to [Delayed Rejection Adaptive Metropolis (DRAM)](http://helios.fmi.fi/~lainema/dram/) MCMC and [Stein variational Newton (SVN)](https://github.com/gianlucadetommaso/Stein-variational-samplers).
Installation steps are thus different in **matlab** and **python**.

**matlab/**

After changing to `matlab/` directory, the general configuration script `install` can be run in order to add subdirectories to the Matlab path, to test and to compile MEX files.
Moreover, this script will ask you whether the parallel *for* loops should be used for different runs of the tests.

External packages are downloaded automatically by `check_tt`, `check_mcmc` and `check_qmc` scripts in the `utils/` directory when the corresponding test is run. Should the automated download fail, the scripts will tell you where you can get the necessary prerequisite manually.
To simply add all subdirectories of TT-IRT to the path, you can run `utils/check_ttirt` script.


**python/**

First, install [ttpy](https://github.com/oseledets/ttpy) following any of the recipes listed at [https://github.com/oseledets/ttpy](https://github.com/oseledets/ttpy) (using *conda*, *pip* or directly the setup script).

The TT-CD sampler package `tt_irt` is located in `python/tt_irt_py` directory. Change there and run
```
python setup.py install [--user]
```
The installation can be checked by `verify.py` script. It will also try to download `puwr.py` for the IACT estimation.

>Warning: the Python module is significantly behind the Matlab version. In particular, SIRT and DIRT are implemented only in Matlab.

### Scripts for running examples

All files for running experiments start with a `test_` prefix. You may click on each item to open an individual directory or file.

**matlab/examples**
 * [`shock_absorber/`](https://github.com/dolgov/TT-IRT/tree/master/matlab/examples/shock_absorber)                   Shock absorber reliability example
   - `test_shock_absorber_tt.m`        TT-MH and TT-qIW methods (single-level)
   - `test_shock_absorber_dram.m`      DRAM test
 * [`diffusion/`](https://github.com/dolgov/TT-IRT/tree/master/matlab/examples/diffusion)                        Inverse diffusion example
   - `test_diffusion_tt.m`             TT-MH/TT-qIW test
   - `test_diffusion_dirt.m`           DIRT test
   - `test_diffusion_dram.m`           DRAM test
   - `test_diffusion_qmcrat.m`         QMC ratio test
 * [`predator_prey/`](https://github.com/dolgov/TT-IRT/tree/master/matlab/examples/predator_prey)                    Example of calibration of the Predator-prey model
   - `test_predator_prey_dirt.m`       DIRT test
   - `test_predator_prey_dram.m`       DRAM test
   - `test_predator_prey_svn.m`        SVN test
 * [`lorenz/`](https://github.com/dolgov/TT-IRT/tree/master/matlab/examples/lorenz)                           A Lorenz-40 identification example
   - `test_lorenz.m`                   DIRT test

Each test can be run without any arguments. In this case, they will be interactively asked from the user. For batch runs, parameters can be passed in pairs of inputs ``'param_name1'``, ``param_value1``, ``'param_name2'``, ``param_value2``, and so on. For example,
```
test_shock_absorber_tt('D', 6, 'x', [], 'log2N', 14, 'delta', 0.05, 'n', 17, 'runs', 8)
```
will run the shock absorber test with all default parameters. Only a subset of arguments can be given, e.g.
```
test_shock_absorber_tt('log2N', 18, 'runs', 32)
```
will take the corresponding values from the inputs, and ask the user only for the remaining parameters.

Each test will print some statistical data (expected values, standard deviations, CPU times and so on).
Moreover, it will create all the variables in the main Matlab workspace, so they can be accessed afterwards.


**python/**
 * Shock absorber example
   - `test_shock_absorber_tt.py` Script for running the shock absorber test with TT-MH, encapsulates prior and likelihood functions
 * Core package in `tt_irt_py` directory
   - `tt_irt1_int32.c`  TT-CD sampler from an approximate density in TT format (uses 32-bit integers)
   - `tt_irt.py`        CTypes wrapper
   - `setup.py`         Compilation and installation script
 * *Utilities*
   - `verify.py`        Check availability of ttpy and tt_irt, download puwr.py

In the shock absorber example, the parameters are set up statically in the `test_shock_absorber_tt.py` file. In order to change them, go to lines 92-96 and edit
```
d = 6       # number of covariates
n = 17      # grid size
tol = 5e-2  # TT stopping threshold
log2N = 16  # log2(number of samples)
runs = 8    # number of runs
```


### Structure of the Matlab repository

Click on each item to see a Readme for that individual directory and files therein.

 * [`constructors`](https://github.com/dolgov/TT-IRT/tree/master/matlab/constructors)   Functions for constructing TT decompositions using cross interpolation and DIRT.
 * [`examples`](https://github.com/dolgov/TT-IRT/tree/master/matlab/examples)       Scripts and auxiliary functions for running numerical experiments.
 * [`samplers`](https://github.com/dolgov/TT-IRT/tree/master/matlab/samplers)       Functions to produce samples. These include Rosenblatt transforms, MC, QMC and MCMC samplers.
 * [`utils`](https://github.com/dolgov/TT-IRT/tree/master/matlab/utils)          Helper functions including interpolations, MEX files and installation checkers.


### Further docs

Navigate into leaf directories to see their individual README files.
In addition, each function file contains an extended description. See e.g.
```
help('tt_irt_sqr')
```
or open the file in the editor.
