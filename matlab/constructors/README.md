# TT-IRT (Matlab version)
Inverse Rosenblatt Transports (IRT) + MCMC sampling using Tensor Train (TT) approximation

## Constructors

This directory contains functions for constructing TT decompositions (or sequences thereof). Each function contains an individual description in the header of the file.

 * `amen_cross_s.m`   Enhanced TT-Cross algorithm for the TT approximation.
 * `tt_dirt_approx.m` Construction of the Deep Inverse Rosenblatt transform (DIRT) structure. It uses `amen_cross_s.m` at each layer (by default), or [`greedy2_cross.m`](https://github.com/oseledets/TT-Toolbox/blob/master/cross/greedy2_cross.m) from the [TT-Toolbox](https://github.com/oseledets/TT-Toolbox) (see the *crossmethod* parameter).
 * `als_cross_parametric.m`      ALS-Cross algorithm for the parametrized PDE solution, see also [http://people.bath.ac.uk/sd901/als-cross-algorithm/](http://people.bath.ac.uk/sd901/als-cross-algorithm/)


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

