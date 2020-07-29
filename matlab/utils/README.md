# TT-IRT (Matlab version)
Inverse Rosenblatt Transports (IRT) + MCMC sampling using Tensor Train (TT) approximation

## Utils

This directory contains functions for carrying out "elementary" steps of IRT and other algorithms, as well as for checking the consistency of the repository and for downloading missing prerequisites.
Each function contains an individual description in the header of the file.

 * `check_mcmc.m`   Check/download/add-to-path for DRAM and UWerr
 * `check_tt.m`     Check/download/add-to-path for TT-Toolbox
 * `check_qmc.m`    Check/download for the QMC lattice generating vector
 * `check_ttirt.m`  Adds TT-IRT subdirectories to the path
 * `lagrange_interpolant.m`  Lagrange interpolation between two sets of points
 * `tt_sample_lagr.m`  Lagrange interpolation of a TT decomposition
 * `tracemultm.m`     Auxiliary file for conditioning of TT blocks in IRT
 * `tracemult.c`      Faster MEX version of the tracemult operation
 * `project_blockdiag_mex.c`     Faster MEX version of the reduction step in ALS-Cross.
 * `solve_blockdiag_mex.c`       Faster MEX version of the local solver in ALS-Cross.
 * `tt_irt_mex.c`     **(deprecated)** MEX version of IRT using linear splines (uses `tt_irt1_int64.c`)

Besides, the directory contains pre-compiled MEX files for Linux x86_64:

 * `project_blockdiag_mex.mexa64`
 * `solve_blockdiag_mex.mexa64`
 * `tracemult.mexa64`
 * `tt_irt_mex.mexa64`

However, if you use a different platform or significantly different versions of Matlab and system libraries,
you may want to recompile MEX files using the `install` script in the main directory.

