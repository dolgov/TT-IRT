# TT-IRT
Inverse Rosenblatt (cumulative distribution) + MCMC sampling using Tensor Train approximation

Algorithms and running codes for the paper "Approximation and sampling of multivariate probability distributions in the tensor train decomposition" [[arXiv:1809.xxxx](http://arxiv.org/abs/)] by Sergey Dolgov, Karim Anaya-Izquierdo, Colin Fox and Robert Scheichl.

The toolbox consists of **matlab** and **python** parts. Moreover, the conditional distribution sampler (TT-CD) is also implemented in C with two types of integers: `matlab/tt_irt1_int64.c` and `python/tt_irt_py/tt_irt1_int32.c` which can be linked into  Matlab MEX, Python CDLL or other project.

### Installation



### Running examples

**matlab/**

 * `test_shock_absorber_tt.m`		Shock absorber test with TT-MH and TT-qIW methods
 * `test_shock_absorber_dram.m`     Shock absorber test with DRAM
 * `test_diffusion_tt.m`            Inverse diffusion test with TT-MH/TT-qIW
 * `test_diffusion_dram.m`          Inverse diffusion test with DRAM
 * `test_diffusion_qmcrat.m`        Inverse diffusion test with QMC ratio estimator

**python/**
 * `test_shock_absorber_tt.py`      Shock absorber test with TT-MH

### Content
