from __future__ import print_function, absolute_import, division
import numpy as np
import tt
import tt.cross
from tt.cross.rectcross import rect_cross
from tt_irt import tt_irt1
import puwr
import time
import sys

# log-prior distribution
def logf_prior(theta, beta_mean, beta_var):
    d = theta.shape[1]-2
    I = theta.shape[0]

    alpha = 6.8757
    beta = 2.2932
    theta2 = theta[:,d+1]
    theta2 = np.reshape(theta2, [I,1], order='F')

    beta_mean_I = np.repeat(np.reshape(beta_mean, [1,d+1]), I, axis=0)
    beta_mean_I = np.reshape(beta_mean_I, [I, d+1], order='F')
    beta_var_I = np.repeat(np.reshape(beta_var, [1,d+1]), I, axis=0)
    beta_var_I = np.reshape(beta_var_I, [I, d+1], order='F')
    theta2_d = np.repeat(theta2, d+1, axis=1)
    theta2_d = np.reshape(theta2_d, [I, d+1], order='F')
    theta_in = theta[:, range(0,d+1)]
    theta_in = np.reshape(theta_in, [I, d+1], order='F')

    F = -((theta_in - beta_mean_I)**2)*0.5*theta2_d/beta_var_I
    F = np.sum(F, axis=1)
    F = np.reshape(F, [I, 1], order='F')
    F = F + (alpha-0.5)*np.log(theta2) - beta*theta2
    return F


# log-Likelihood
def logL_weibull(theta, x, y, c):
    d = theta.shape[1]-2
    I = theta.shape[0]
    beta = theta[:, range(0,d+1)]
    beta = np.reshape(beta, [I, d+1], order='F')
    lam = np.reshape(theta[:, d+1], [I,1], order='F')
    F = np.zeros(I, dtype=np.float64)
    F = np.reshape(F, [I,1], order='F')
    m = np.size(y)
    x = np.reshape(x, [d,m], order='F')
    beta0 = np.reshape(beta[:,0], [I,1], order='F')
    for i in range(0,m):
        logeta = beta[:,range(1,d+1)].dot(x[:,i])
        logeta = np.reshape(logeta, [I,1], order='F')
        logeta = logeta + beta0
        eta = np.exp(logeta)
        yeta = y[i]/eta
        if c[i] == 1:  # censored
            f = -(yeta**lam)
        else:
            f = np.log(lam) - logeta + (lam-1.0)*(np.log(y[i]) - logeta) - (yeta**lam)
            f = f + np.log(30000.0) # to prevent underflow

        F = F + f
    return F


# Full index function for cross
def crossfun(ind, theta, x, y, censind, beta_mean, beta_var):
    d = ind.shape[1]
    I = ind.shape[0]

    # Cross indices are integers mathematically,
    # but in rect_cross they are sometimes stored as doubles.
    # Fix that
    ind = ind.astype(np.int32)

    thetai = np.zeros(I*d)
    thetai = np.reshape(thetai, [I,d], order='F')
    for i in range(0,d):
        thetai[:,i] = theta[i][ind[:,i]]
    thetai = np.reshape(thetai, [I,d], order='F')
    # Posterior density
    F = np.exp(logf_prior(thetai, beta_mean, beta_var) + logL_weibull(thetai, x, y, censind))
    F = np.reshape(F, [I,1], order='F')
    return F


###################### MAIN ###############

d = 6
n = 16
tol = 5e-2
lvls = 16
nruns = 8

y = [6700, 6950, 7820, 8790, 9120, 9660, 9820, 11310, 11690, 11850, 11880, 12140, 12200,
     12870, 13150, 13330, 13470, 14040, 14300, 17520, 17540, 17890, 18420, 18960,
     18980, 19410, 20100, 20100, 20150, 20320, 20900, 22700, 23490, 26510, 27410,
     27490, 27890, 28100]
censind = [0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 0, 1, 0, 1, 0, 1, 1]

# Estimate the ranges of beta
beta_mean = np.zeros(d+1)
beta_mean[0] = np.log(30796.0)
beta_var = np.ones(d+1)
beta_var[0] = 0.1563

a = beta_mean - 3.0*np.sqrt(beta_var)
b = 2.0*beta_mean - a
a = np.hstack([a, 0.0])
b = np.hstack([b, 13.0])

# Vectors of grid nodes
theta = (d+2)*[None]
theta_f = (d+2)*[None]
for i in range(0,d+2):
    h = (b[i]-a[i])/n
    theta[i] = a[i]+np.arange(1,n+1)*h
    theta_f[i] = a[i]+np.arange(0,n+1)*h
    theta_f[i] = np.reshape(theta_f[i], [n+1, 1], order='F')

theta_f = np.hstack(theta_f)
theta_f = np.reshape(theta_f, [(d+2)*(n+1), 1], order='F')

# Simulate artificial covariates
#x = np.random.randn(d*np.size(y))/d;
x = np.ones(d*np.size(y))/d
x = np.reshape(x, [d, np.size(y)], order='F')

# short index function for cross
myfun = lambda ind: crossfun(ind, theta, x, y, censind, beta_mean, beta_var)

print('Running '+repr(nruns)+' tests for '+repr(d)+' covariates, producing '+repr(2**lvls)+' samples')
print('')

Q_py = nruns*[None]
tau = nruns*[None]
for irun in range(0,nruns):
    # Run cross
    f0 = tt.rand(n, d+2, 8)
    f = rect_cross.cross(myfun, f0, nswp=50, kickrank=2, rf=2, eps=tol)
    print(f)


    # Seed points
    M = 2**lvls
    q = np.random.random([M, d+2])
    q = np.reshape(q, [M, d+2], order = 'F')


    ttt = time.time()
    # Sample
    Z, Pz = tt_irt1(q, f, theta_f)
    ttimes_invcdf = time.time() - ttt
    print('Sampling time = '+repr(ttimes_invcdf))

    Z = np.reshape(Z, [M,d+2], order='F') # Proposed samples
    Pz = np.reshape(Pz, [M,1], order='F') # Approximate density

    # Exact density
    Pex = np.exp(logL_weibull(Z, x, y, censind) + logf_prior(Z, beta_mean, beta_var))

    # MCMC rejections
    num_of_rejects = 0
    for i in range(0,M-1):
        alpha = np.exp(np.log(Pex[i+1]) - np.log(Pex[i]) + np.log(Pz[i]) - np.log(Pz[i+1]))
        if alpha<np.random.random(1): # reject
            Z[i+1,:] = Z[i,:]
            Pex[i+1] = Pex[i]
            Pz[i+1] = Pz[i]
            num_of_rejects = num_of_rejects + 1

    print('rejection rate = ' + repr((num_of_rejects*100.0)/M))

    # QoI = Quantile function
    q = 0.05
    theta1 = np.exp(Z[:,0])
    theta2 = Z[:,d+1]
    q_post = theta1*((-np.log(q))**(1.0/theta2))
    # IACT -- only works in python 2.7
    if sys.version_info<(3,0):
        dumpmean, delta, tint, d_tint = puwr.tauint(np.reshape(q_post, [1, 1, M]), 0)
        tau[irun] = tint*2.0

    print('tau = ' + repr(tau[irun]))
    q_post = np.mean(q_post)
    print('Mean quantile = ' + repr(q_post))
    Q_py[irun] = q_post


print('')
print('TT Shock absorber completed. Some average values:')
print('Q_py = '+repr(np.mean(Q_py))+' +- '+repr(np.sqrt(np.sum((Q_py - np.mean(Q_py))**2)/(nruns-1))) )
print('tau = '+repr(np.mean(tau))+' +- '+repr(np.sqrt(np.sum((tau - np.mean(tau))**2)/(nruns-1))) )

