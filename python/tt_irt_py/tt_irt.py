from ctypes import cdll, c_int, c_double, POINTER
import numpy as np
import os
import glob

import tt

libfile = os.path.dirname(__file__)
libfile = os.path.join(libfile, "tt_irt1*")
libfile = glob.glob(libfile)[0]
lib = cdll.LoadLibrary(libfile)

def tt_irt1(q, f, xsf):
    """ Inverse Rosenblatt sampler, linear splines
        Inputs:
          q: seed samples from [0,1]^d (np.float64 M x d Fortran shaped)
          f: tt.tensor of the PDF (dimension d), constructed on a grid without left boundaries
          xsf: vector of grid points of all variables stacked together (np.float64, size sum(f.n+1))
        Returns:
          Z: transformed samples (np.float64 M x d Fortran shaped)
          Pz: values of the sampling density at Z (np.float64 M x 1)
    """
    c_irt1 = lib.tt_irt1
    c_irt1.restype = None
    c_irt1.argtypes = [c_int, POINTER(c_int), POINTER(c_double), POINTER(c_int), POINTER(c_double), c_int, POINTER(c_double), POINTER(c_double), POINTER(c_double)]
    #                  d           n              xsf              ttrank          ttcore            M           q                 Z                 Pz

    # Cores of f must be extracted carefully, since we might have discontinuous ps
    core = np.zeros((f.core).size, dtype=np.float64)
    ps_my = 0
    for i in range(0,f.d):
        cri = f.core[range(f.ps[i]-1,f.ps[i+1]-1)]
        core[range(ps_my,ps_my+f.r[i]*f.n[i]*f.r[i+1])] = cri
        ps_my = ps_my + f.r[i]*f.n[i]*f.r[i+1]

    d = c_int(f.d)
    n = (np.array(f.n)).ctypes.data_as(POINTER(c_int))
    xsfp = xsf.ctypes.data_as(POINTER(c_double))
    rf = (np.array(f.r)).ctypes.data_as(POINTER(c_int))
    corep = core.ctypes.data_as(POINTER(c_double))
    M = c_int(q.shape[0])
    qp = q.ctypes.data_as(POINTER(c_double))
    
    Z = np.zeros([q.shape[0], q.shape[1]], dtype=np.float64, order='F')
    Pz = np.zeros([q.shape[0]], dtype=np.float64, order='F')

    Zp = Z.ctypes.data_as(POINTER(c_double))
    Pzp = Pz.ctypes.data_as(POINTER(c_double))

    # Sampler is actually here
    c_irt1(d, n, xsfp, rf, corep, M, qp, Zp, Pzp)

    return (Z, Pz)
    
