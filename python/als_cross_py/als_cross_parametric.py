import warnings
import tt
import numpy
import time


def localcross(Y,tol, return_indices=False):
  """
  Full-pivoted cross for truncating one ket TT block instead of SVD

  :param Y: TT block
  :param tol: truncation tolerance
  :param return_indices: (optional) return indices if true

  :return: 
    if ``return_indices=True`` returns tuple u,v,I; else returns tuple u,v
  """
  # TODO fast implementation (or find in ttpy)
  if len(numpy.shape(Y)) == 2:
    n,m = numpy.shape(Y)
    b = 1
  else:
    n,m,b = numpy.shape(Y)

  minsize = min(n, m*b)
  u = numpy.zeros((n,minsize))
  v = numpy.zeros((minsize, m*b))
  res = numpy.reshape(Y, (n, m*b), order='F')
    
  I = numpy.zeros((minsize)) # return also indices

  val_max = numpy.max(numpy.abs(Y))
  r_tol = 1 # rank after truncation (0 tensor is also rank 1)
  for r in range(minsize):
    # res = numpy.reshape(res, (-1,1), order='F')
    piv = numpy.argmax(res)
    piv = numpy.unravel_index(piv, numpy.shape(res))
    val = res[piv]
    if val <= tol*val_max:
      break

    r_tol = r+1
    u[:,r] = res[:,piv[1]]
    v[r,:] = res[piv[0],:] / val
    res -= numpy.outer(u[:,r], v[r, :])
    I[r] = piv[0]

  I = I[:r_tol]
  u = u[:,:r_tol]
  v = v[:r_tol,:]

  # QR u
  u, rv = numpy.linalg.qr(u)
  v = numpy.matmul(rv, v)

  if return_indices:
    return u,v,numpy.reshape(I, (1,-1))
  else:
    return u,v



def als_cross_parametric(coeff, assem_solve_fun, tol , **varargin):
  """
  TT ALS-Cross algorithm.

  :param coeff: (d+1)-dimensional block TT format storing coefficients, with the
  first TT rank Mc corresponding to different coefficients, and the 
  entire first TT block corresponding to the deterministic part with
  Nxc degrees of freedom in the deterministic variable. The other 
  TT blocks correspond to the d (stochastic) parameters.

  :param assem_solve_fun: a handle to the function which should take the
    coefficients Ci of size ``[Mc, Nxc, r]`` and return lists ``U,A,F``, 
    respectively solutions, matrices and RHS at the given Ci. 
    Here r is the TT rank of the solution in the first TT block.
    ``U,A,F`` must be cell arrays of size 1 x r, each snapshot ``U[i]`` must be
    a column vector of size Nxu x 1, each snapshot ``F[i]`` must be a 
    column vector of size Nxa x 1, and each matrix ``A[i]`` must be a
    Nxa x Nxa matrix.

    Alternatively, if ``use_indices`` is enabled, assem_solve_fun 
    should take indices of parameters where
    the systems must be solved, in the form of r x d integer matrix,
    with elements in the k-th column ranging from 1 to n_k, the number
    of grid points in the k-th parameter. The output format is the same.

    !!! In both cases, A and F should depend near linearly on coeff !!!

  :param tol: cross truncation and stopping tolerance

  :param **varargin: Optional parameters:
    - ``Pua``: a matrix to project spatial block of solution to spat. block of matrix
      For good performance, Pua should be a full rank projector, with
      size(Pua)==[Nxa,Nxu] with Nxu>=Nxa.
      Default empty Pua assumes Nxu==Nxa.
    - ``nswp``: max number of iterations (default 5)
    - ``kickrank``: max TT rank of the residual/enrichment (default 10)
    - ``random_init``: if greater than 0, take random_init random indices at
      start; if 0 (default), take maxvol indices of coeff
    - ``used_indices``: selects the type of input for assem_solve_fun:
      ``False`` (default) assumes that the function takes values of
      the coefficients,
      ``True`` assumes that the function takes indices of parameters

  :returns:
    u: solution in TT format
    time_extern: list [time solve, time project] of cpu times
    funevals: total number of deterministic solves

  """
  # TODO paper uses c-like order, matlab fortran-like order. Problem? 

  # default values for optional parameters
  nswp = 5 # number of iterations
  kickrank = 10 # The base for enrichment rank. The actual ranks are scaled to the coeff. ranks
  Pua = []  # A matrix that maps spatial DOFS of the solution to the spatial DOFs of the matrix/rhs
  random_init = 0 # If >0, start with random indices, instead of the coefficient
  use_indices = False # whether assem_solve_fun takes indices (or values)

  # parse parameters
  for (arg, val) in varargin.items():
    match arg:
      case 'nswp':
        nswp = val
      case 'kickrank':
        kickrank = val
      case 'pua':
        Pua = val
      case 'random_init':
        random_init = val
      case 'use_indices':
        use_indices = val
      case _:
        warnings.warn('unknown argument \'' + arg + '\'', SyntaxWarning)

  # rng
  rng = numpy.random.default_rng()
  
  # all grid sizes
  Nxc = coeff.n[0]  # spatial grid size
  ny = coeff.n[1:]  # parametric grid sizes
  d = coeff.d - 1   # d is parameter dimension only
  Mc = coeff.r[0]   # number of coefficient components
  rc = coeff.r[1:]  
  ru = numpy.copy(rc)          # these will be TT ranks of the solution
  coeff_old = coeff.copy()
  coeff = tt.vector.to_list(coeff)
  C0 = coeff[0]
  coeff = coeff[1:]

  # Prepare storage for reduction/sampling matrices
  if nswp > 1:
    UAU = [None] * (d+1)  # left Galerkin reductions, matrix
    UF = [None] * (d+1)   # left Galerkin reductions, RHS
  UC = [None] * (d+1) # right cross samples of C on U-indices
  UC[-1] = numpy.ones((1,1))
  # Prepare storage for the residual
  if kickrank > 0:
    ZZ = [None] * (d+1)
    ZZ[-1] = numpy.ones((1,1))
    ZU = ZZ.copy() # ZAU from the left, ZU from the right
    ZC = ZZ.copy()
    rz = numpy.round(kickrank * rc / numpy.max(rc)).astype(int)
    rz[rz<1] = 1
    rz[-1] = 1

  xi = numpy.ones((1, random_init))
  if use_indices:
    # Initialise global indices if the user function works with them
    Ju = []

  # First, orthogonalize the coefficient.
  # We can derive its optimal indices (as an initial guess), or use random
  # ones
  v = numpy.ones((1,1))
  for i in reversed(range(d)):
    crc = numpy.reshape(coeff[i], (rc[i]*ny[i], -1), order='F')
    crc = numpy.matmul(crc, numpy.transpose(v))
    crc = numpy.reshape(crc, (rc[i], ny[i]*rc[i+1]), order='F')
    crc = numpy.transpose(crc)
    (crc,v) = numpy.linalg.qr(crc)
    # v = v[:, :rc[i]] # TODO useless?
    rc[i] = numpy.shape(crc)[1]
    ind = tt.maxvol.maxvol(crc) 
    # ind = numpy.arange(rc[i])
    crc = numpy.transpose(crc) # TODO reorder ok
    CC = crc[:, ind]
    crc = numpy.linalg.solve(CC, crc)
    v = numpy.matmul(numpy.transpose(CC), v)
    coeff[i] = numpy.reshape(crc, (rc[i], ny[i], rc[i+1]), order='F')

    if use_indices:
      Ju = numpy.hstack(
        numpy.tile(numpy.arange(ny[i]), [rc[i+1], 1]),
        numpy.repeat(Ju, ny[i], axis=0)
      ) # TODO verify correctness
      Ju = Ju[ind, :]
      
    if random_init > 0 and i > 0:
      # Random sample coeff from the right
      ind = rng.integers(ny[i], size=(random_init))
      # ind = [0]
      xi_new = numpy.zeros((rc[i], random_init))
      for k in range(random_init):
        xi_new[:, k] = (coeff[i][:, ind[k], :].reshape(rc[i],rc[i+1])) @ xi[:rc[i+1], k]
      # xi = numpy.einsum('i...j,j...->i...', coeff[i][:, ind, :], xi) # TODO verify this works
      xi = xi_new.copy()
      UC[i] = xi
      ru[i] = random_init
    else:
      UC[i] = numpy.eye(rc[i])
      ru[i] = rc[i]

    if kickrank > 0:
      # Initialize the residual
      crz = rng.standard_normal((ny[i] * rz[i+1], rz[i]))
      crz = numpy.linalg.qr(crz)[0]
      rz[i] = numpy.shape(crz)[1]
      ind = tt.maxvol.maxvol(crz)
      crz = numpy.transpose(crz)
      # Sample the coefficient and solution at Z-indices
      ZC[i] = numpy.matmul(numpy.reshape(coeff[i], (rc[i]*ny[i], rc[i+1]), order='F'), ZC[i+1])
      ZC[i] = numpy.reshape(ZC[i], (rc[i], ny[i] * rz[i+1]), order='F')
      ZC[i] = ZC[i][:, ind]
      ZU[i] = ZC[i] # in the first iteration this is the same
  
  # Init empty solution storage
  u = [None] * d
  U0 = []

  # The coefficient+rhs at sampled indices
  C0 = numpy.reshape(C0, (Mc*Nxc,-1), order='F') # size Mc*Nxc, rc1
  C0 = numpy.matmul(C0, numpy.transpose(v))
  # This is the spatial block of the coefficient, in the representation when
  # all parametric blocks contain identity matrices.
  # The coeff will not change anymore
  C0 = numpy.reshape(C0, (Mc, Nxc, rc[0]), order='F')

  coeff_new = tt.vector.from_list([C0]+coeff)
  print("err_c = "+str((coeff_old - coeff_new).norm()))
  
  # Initialise cost profilers
  time_solve = 0
  time_project = 0
  funevals = 0 # init number of det. solves

  # Start iterating
  i = 0
  max_dx = 0
  dir = 1
  swp = 1
  while swp <= nswp:
    if i == 0:
      ##### Work on spatial block ###################################
      # Previous guess (if we have one)
      Uprev = U0
      # solve deterministic problems at the U-indices
      if use_indices:
        Ci = Ju
      else:
        # Construct the coeff there
        Ci = numpy.matmul(numpy.reshape(C0, (Mc*Nxc,-1), order='F'), UC[0]) # size Mc*Nxc, rc[0]
        Ci = numpy.reshape(Ci, (Mc, Nxc, ru[0]), order='F')

      t1__uc = time.perf_counter()
      
      print("Ci")
      print(Ci)

      if swp == 1:
        U0, A0s, F0 = assem_solve_fun(Ci)
        Nxa = numpy.shape(A0s[0])[0]
        F0 = numpy.hstack(F0)
        # In the first sweep, Ci==C0, and we need the corresponding
        # matrices (A0s) and RHS (F0), since we'll use them in
        # subsequent iterations ...
      else:
        # ... where it's enough to compute the solution only
        U0 = assem_solve_fun(Ci)[0]

      time_solve += time.perf_counter() - t1__uc
      del Ci # memory saving
      funevals += ru[0]
      U0 = numpy.hstack(U0)
      Nxu = numpy.shape(U0)[0] # this, again, can differ from Nxa or Nxc
      if Nxu != Nxa and len(Pua) == 0:
        raise RuntimeError('Numbers of spatial DOFs in u and A differ, and no transformation matrix is given. Unable to reduce model')
      
      # check the error
      if len(Uprev) > 0:
        dx = numpy.linalg.norm(U0 - Uprev) / numpy.linalg.norm(U0)
      else:
        dx = 1
      max_dx = max(max_dx, dx)

      print(f'=als_cross_parametric= 0 swp={swp}, max_dx={max_dx:.3e}, max_rank = {max(ru)}')

      # Unusual ALS exit condition: after solving for the spatial block,
      # to have the non-orth center here.
      if max_dx < tol or swp > nswp: # TODO second condition useless?
        break

      max_dx = 0

      print("U0")
      print(U0)

      # Truncate U0 via full-pivoted cross approx
      # U0, v = localcross(U0, tol/numpy.sqrt(d)) # TODO find suitable ttpy function
      U0, v = numpy.linalg.qr(U0)
      if swp > 1:
        u[0] = numpy.reshape(u[0], (ru[0], ny[0]*ru[1]), order='F')
        u[0] = numpy.matmul(v, u[0])
      ru[0] = numpy.shape(U0)[1]

      if kickrank > 0:
        # Compute residual
        # Compute A * U at Z indices
        cru = numpy.linalg.multi_dot((U0, v, ZU[0]))
        if Nxa != Nxu:
          cru = numpy.matmul(Pua, cru)

        Z0 = numpy.zeros((Nxa, rz[0]))
        for j in range(rz[0]):
          crA = A0s[0] * ZC[0][0,j]
          for k in range(1,rc[0]):
            crA += A0s[k] * ZC[0][k,j]

          Z0[:,j] = numpy.matmul(crA, cru[:,j])
        
        Z0 -= numpy.matmul(F0, ZC[0])
        Z0 = numpy.linalg.qr(Z0)[0]
        rz[0] = numpy.shape(Z0)[1]
        if Nxa != Nxu:
          cru = numpy.hstack((U0, numpy.matmul(numpy.conjugate(numpy.transpose(Pua)), Z0))) # TODO why conjugate
        else:
          cru = numpy.hstack((U0, Z0))

        # QR U0
        U0, v = numpy.linalg.qr(cru)
        v = v[:,:ru[0]]
        if swp > 1:
          u[0] = numpy.reshape(u[0], (ru[0], ny[0]*ru[1]), order='F')
          u[0] = numpy.matmul(v, u[0])

        ru[0] = numpy.shape(U0)[1]

      # Project the model onto the solution basis U0
      # UAU is very large here
      # Need to run some loops to save mem
      t1__uc = time.perf_counter()
      UAU_new = [None] * rc[0]
      Uprev = U0
      if Nxa != Nxu:
        Uprev = numpy.matmul(Pua, U0)
      
      for j in range(rc[0]):
        UAU_new[j] = numpy.linalg.multi_dot((numpy.conjugate(numpy.transpose(Uprev)), A0s[j], Uprev)) # TODO why conjugate
        UAU_new[j] = numpy.reshape(UAU_new[j], (-1,1), order='F')
      
      if nswp == 1:
        # we don't need to save UAU projections for all blocks, if we
        # will never iterate back
        UAU = numpy.hstack(UAU_new)
        UF = numpy.matmul(numpy.conjugate(numpy.transpose(Uprev)), F0) # TODO why conjugate
      else:
        UAU[0] = numpy.hstack(UAU_new)
        UF[0] = numpy.matmul(numpy.conjugate(numpy.transpose(Uprev)), F0) # TODO why conjugate

      time_project += time.perf_counter() - t1__uc
      del UAU_new

      # Project onto residual
      if kickrank > 0:
        # Project onto residual basis
        ZU_new = [None] * rc[0]
        for j in range(rc[0]):
          ZU_new[j] = numpy.linalg.multi_dot((numpy.conjugate(numpy.transpose(Z0)), A0s[j], Uprev)) # TODO why conjugate
          ZU_new[j] = numpy.reshape(ZU_new[j], (-1,1), order='F')
        
        ZU[0] = numpy.hstack(ZU_new)
        ZC[0] = numpy.matmul(numpy.conjugate(numpy.transpose(Z0)), F0) # TODO why conjugate
        del ZU_new
      
      # Save some memory if we only have 1 iteration
      if nswp == 1:
        del A0s
        del F0

    else:  ##### End i == 0, Loop for reduced system ################
      # TODO indices are shifted by -1, reconsider
      # Solve block-diagonal reduced system
      crC = numpy.reshape(coeff[i-1], (rc[i-1]*ny[i-1], rc[i]), order='F')
      crC = numpy.matmul(crC, UC[i]) # now these are indices from the right, and Galerkin from the left
      crC = numpy.reshape(crC, (rc[i-1], ny[i-1]* ru[i]), order='F')
      if nswp == 1:
        UAUi = UAU # cell-free local reduction from the left
        UFi = UF
      else:
        UAUi = UAU[i-1]
        UFi = UF[i-1]
      
      crF = numpy.matmul(UFi, crC)
      try:
        raise NotImplementedError('No fast solver available') # TODO implement
      except:
        crA = numpy.reshape(UAUi, (ru[i-1]*ru[i-1], rc[i-1]), order='F')
        crC = numpy.reshape(crC, (rc[i-1], ny[i-1]*ru[i]), order='F')
        cru = [None] * (ny[i-1] * ru[i]) # TODO do these with slicing of ndarray?
        for j in range(len(cru)):
          Ai = numpy.matmul(crA, crC[:,j])
          Ai = numpy.reshape(Ai, (ru[i-1], ru[i-1]))
          cru[j] = numpy.linalg.solve(Ai, crF[:,j]).reshape(-1,1)
        
        cru = numpy.hstack(cru) 

      cru = numpy.reshape(cru, (ru[i-1]*ny[i-1],ru[i]), order='C')
      print("cru shape = "+str(cru.shape[0]) + " "+str(cru.shape[1]))

      # error check
      if u[i-1] is None:
        dx = 1
      else:
        dx = numpy.linalg.norm(cru - numpy.reshape(u[i-1], (ru[i-1]*ny[i-1], ru[i]), order='F')) / numpy.linalg.norm(cru)
      max_dx = max(max_dx, dx)
      u[i-1] = numpy.reshape(cru, (ru[i-1], ny[i-1], ru[i])) 
      # if we're in the d-th block (i == d-1), we're done

      if i < d and dir > 0: ##### left-right sweep ################
        # Truncate cru with cross
        # cru, v = localcross(cru, tol/numpy.sqrt(d))
        cru, v = numpy.linalg.qr(cru)
        rv = v
        if kickrank > 0:
          # Update the residual and enrichment
          crC = numpy.reshape(coeff[i-1],(rc[i-1]*ny[i-1], rc[i]), order='F')
          crC = numpy.matmul(crC, ZC[i]) # now these are indices from the right, and Galerkin from the left
          crC = numpy.reshape(crC, (rc[i-1], ny[i-1]*rz[i]), order='F')
          Uprev = numpy.linalg.multi_dot((cru, rv, ZU[i]))
          Uprev = numpy.reshape(Uprev, (ru[i-1],ny[i-1]*rz[i]), order='F')
          # first enrichment
          crA = numpy.reshape(UAUi, (ru[i-1]*ru[i-1], rc[i-1]), order='F')
          crz = [None] * (ny[i-1]*rz[i])
          for j in range(len(crz)):
            Ai = numpy.matmul(crA, crC[:,j])
            Ai = numpy.reshape(Ai, (ru[i-1], ru[i-1]))
            crz[j] = numpy.matmul(Ai, Uprev[:,j])
            crz[j] = numpy.reshape(crz[j], (-1,1))
          
          crz = numpy.hstack(crz) # size ru1 x n*rz2
          crz -= numpy.matmul(UFi, crC)
          crz = numpy.reshape(crz, (ru[i-1]*ny[i-1],rz[i]), order='F')
          v = numpy.hstack((cru, crz))
          cru = v # TODO useless ?
          # QR u
          cru, v = numpy.linalg.qr(cru)
          v = v[:, :numpy.shape(rv)[0]]
          rv = numpy.matmul(v, rv)
          # Now the residual itself
          crA = numpy.reshape(ZU[i-1], (rz[i-1]*ru[i-1],rc[i-1]), order='F')
          crz = [None] * (ny[i-1]*rz[i])
          for j in range(len(crz)):
            Ai = numpy.matmul(crA, crC[:,j])
            Ai = numpy.reshape(Ai, (rz[i-1],ru[i-1]), order='F')
            crz[j] = numpy.matmul(Ai, Uprev[:,j])
            crz[j] = numpy.reshape(crz[j], (-1,1))

          crz = numpy.hstack(crz) # size rz1 x n*rz2
          crz -= numpy.matmul(ZC[i-1], crC)
          crz = numpy.reshape(crz, (rz[i-1]*ny[i-1], rz[i]), order='F')
        
        # cast the non-orth factor to the next block
        if swp > 1:
          u[i] = numpy.reshape(u[i], (ru[i], ny[i]*ru[i+1]), order='F')
          u[i] = numpy.matmul(rv, u[i])
          u[i] = numpy.reshape(u[i], (-1, ny[i], ru[i+1]), order='F')
        
        ru[i] = numpy.shape(cru)[1]
        u[i-1] = numpy.reshape(cru, (ru[i-1], ny[i-1], ru[i]), order='F')

        # projection from the left -- Galerkin
        # with matrix
        crC = numpy.reshape(coeff[i-1], (rc[i-1],ny[i-1],rc[i]), order='F')
        cru = numpy.reshape(cru, (ru[i-1], ny[i-1], ru[i]), order='F')
        try:
          raise NotImplementedError("No fast projection available")
        except:
          UAU_new = numpy.zeros((ru[i], ru[i]*rc[i]))
          UAUi = numpy.reshape(UAUi, (ru[i-1], ru[i-1]*rc[i-1]), order='F')
          crC = numpy.transpose(crC, (0,2,1))
          cru = numpy.transpose(cru, [0,2,1])
          for j in range(ny[i-1]):
            v = cru[:,:,j]
            crA = numpy.matmul(numpy.conjugate(numpy.transpose(v)), UAUi) # TODO why conj
            crA = numpy.reshape(crA, (ru[i]*ru[i-1], rc[i-1]), order='F')
            crA = numpy.matmul(crA, crC[:,:,j])
            crA = numpy.reshape(crA, (ru[i], ru[i-1]*rc[i]), order='F')
            crA = numpy.transpose(crA)
            crA = numpy.reshape(crA, (ru[i-1], ru[i]*rc[i]), order='F')
            crA = numpy.matmul(numpy.conjugate(numpy.transpose(v)), crA) # TODO why conj
            crA = numpy.reshape(crA, (ru[i]*rc[i], ru[i]), order='F')
            crA = numpy.transpose(crA)
            UAU_new += crA

          cru = numpy.transpose(cru, (0,2,1))
        
        # Save some mem
        if nswp == 1:
          UAU = UAU_new
        else:
          UAU[i] = UAU_new
        
        del UAU_new
        del UAUi

        # with RHS
        crC = numpy.reshape(coeff[i-1], (rc[i-1], ny[i-1]*rc[i]), order='F')
        UFi = numpy.matmul(UFi, crC)
        UFi = numpy.reshape(UFi, (ru[i-1]*ny[i-1], rc[i]), order='F')
        cru = numpy.reshape(cru, (ru[i-1]*ny[i-1], ru[i]), order='F')
        UFi = numpy.matmul(numpy.conjugate(numpy.transpose(cru)), UFi) # TODO why conj
        
        if nswp == 1:
          UF = UFi
        else:
          UF[i] = UFi
        
        del UFi

        # Projections with Z
        if kickrank > 0:
          crz = numpy.linalg.qr(crz)[0]
          rz[i] = numpy.shape(crz)[1]
          # with matrix
          crC = numpy.reshape(crC, (rc[i-1], ny[i-1], rc[i]), order='F')
          cru = numpy.reshape(cru, (ru[i-1], ny[i-1], ru[i]), order='F')
          crz = numpy.reshape(crz, (rz[i-1], ny[i-1], rz[i]), order='F')
          ZU[i] = numpy.zeros((rz[i], ru[i]*rc[i]))
          ZU[i-1] = numpy.reshape(ZU[i-1], (rz[i-1], ru[i-1]*rc[i-1]), order='F')
          crC = numpy.transpose(crC, (0,2,1))
          cru = numpy.transpose(cru, (0,2,1))
          crz = numpy.transpose(crz, (0,2,1))
          for j in range(ny[i-1]):
            v = crz[:,:,j]
            crA = numpy.matmul(numpy.conjugate(numpy.transpose(v)), ZU[i-1]) # TODO why conj
            crA = numpy.reshape(crA, (rz[i]*ru[i-1], rc[i-1]), order='F')
            crA = numpy.matmul(crA, crC[:,:,j])
            crA = numpy.reshape(crA, (rz[i], ru[i-1]*rc[i]), order='F')
            crA = numpy.transpose(crA)
            crA = numpy.reshape(crA, (ru[i-1], rz[i]*rc[i]), order='F')
            v = cru[:,:,j]
            crA = numpy.matmul(numpy.conjugate(numpy.transpose(v)), crA) # TODO why conj
            crA = numpy.reshape(crA, (ru[i]*rc[i], rz[i]), order='F')
            crA = numpy.transpose(crA)
            ZU[i] += crA

          crz = numpy.transpose(crz, (0,2,1))
          crz = numpy.reshape(crz, (rz[i-1]*ny[i-1], rz[i]), order='F')
          # with RHS
          crC = numpy.reshape(coeff[i-1], (rc[i-1],ny[i-1]*rc[i]), order='F')
          ZC[i] = numpy.matmul(ZC[i-1], crC)
          ZC[i] = numpy.reshape(ZC[i], (rz[i-1]*ny[i-1], rc[i]), order='F')
          ZC[i] = numpy.matmul(numpy.conjugate(numpy.transpose(crz)), ZC[i]) # TODO why conj
          if nswp == 1:
            ZC[i-1] = None
            ZU[i-1] = None

      elif dir < 0: ##### right-left sweep
        cru = numpy.reshape(cru, (ru[i-1], ny[i-1]*ru[i]), order='F')
        # truncate cru with cross
        v, cru = localcross(cru, tol/numpy.sqrt(d))
        # now cru is not orthogonal
        rv = v # TODO v useless, asign directly?
        if kickrank > 0:
          # Update the residual and enrichment
          # enrichment first
          crC = numpy.reshape(coeff[i-1], (rc[i-1]*ny[i-1], rc[i]), order='F')
          crC = numpy.matmul(crC, UC[i])
          crC = numpy.reshape(crC, (rc[i-1], ny[i-1]*ru[i]), order='F')
          Uprev = numpy.matmul(rv, cru) # TODO ok to skip reshape?
          Uprev = numpy.reshape(Uprev, (ru[i-1],ny[i-1]*ru[i]), order='F')
          crA = numpy.reshape(ZU[i-1], (rz[i-1]*ru[i-1], rc[i-1]), order='F')
          crz = [None] * (ny[i-1]*ru[i])
          for j in range(len(crz)):
            Ai = numpy.matmul(crA, crC[:,j])
            Ai = numpy.reshape(Ai, (rz[i-1], ru[i-1]), order='F')
            crz[j] = numpy.matmul(Ai, Uprev[:,j])
            crz[j] = numpy.reshape(crz[j], (-1,1))

          crz = numpy.hstack(crz)
          crz -= numpy.matmul(ZC[i-1], crC)
          # crz = numpy.reshape(crz, (rz[i-1]*ny[i-1],ru[i]), order='F') # TODO ok to skip reshape
          crz = numpy.reshape(crz, (rz[i-1], ny[i-1]*ru[i]), order='F')  # TODO crz should already have this shape
          v = numpy.vstack((cru, crz))
          # now the residual itself
          crC = numpy.reshape(coeff[i-1], (rc[i-1]*ny[i-1], rc[i]), order='F')
          crC = numpy.matmul(crC, ZC[i]) # now these are indices from the right, and Galerkin from the left
          crC = numpy.reshape(crC, (rc[i-1], ny[i-1]*rz[i]), order='F')
          # Uprev = numpy.matmul(rv, cru) # TODO Uprev should already have this content 
          Uprev = numpy.reshape(Uprev, (ru[i-1]*ny[i-1], ru[i]), order='F')
          Uprev = numpy.matmul(Uprev, ZU[i])
          Uprev = numpy.reshape(Uprev, (ru[i-1], ny[i-1]*rz[i]), order='F')
          crz = [None] * (ny[i-1]*rz[i])
          for j in range(len(crz)):
            Ai = numpy.matmul(crA, crC[:,j])
            Ai = numpy.reshape(Ai, (rz[i-1], ru[i-1]), order='F')
            crz[j] = numpy.matmul(Ai, Uprev[:,j])
            crz[j] = numpy.reshape(crz[j],(-1,1))
          
          crz = numpy.hstack(crz)
          crz -= numpy.matmul(ZC[i-1], crC)
          crz = numpy.reshape(crz, (rz[i-1]*ny[i-1],rz[i]), order='F')
          cru = v # TODO why not asign cru directly above?

        ru[i-1] = numpy.shape(rv)[1]
        # QR u
        cru, v = numpy.linalg.qr(numpy.transpose(cru))
        v = v[:, :ru[i-1]]
        rv = numpy.matmul(rv, numpy.transpose(v))
        # maxvol to determine new indices
        ind = tt.maxvol.maxvol(cru)
        cru = numpy.transpose(cru)
        UU = cru[:, ind]
        cru = numpy.linalg.solve(UU, cru)
        rv = numpy.matmul(rv, UU)
        ru[i-1] = numpy.shape(rv)[1]
        # Cast non-orth factor to the next block
        if i > 1:
          u[i-2] = numpy.reshape(u[i-2], (ru[i-2]*ny[i-2],-1), order='F')
          u[i-2] = numpy.matmul(u[i-2], rv)
          u[i-2] = numpy.reshape(u[i-2], (ru[i-2], ny[i-2], -1), order='F')
        else:
          U0 = numpy.matmul(U0, rv)

        ru[i-1] = numpy.shape(cru)[0]
        u[i-1] = numpy.reshape(cru, (ru[i-1], ny[i-1], ru[i]), order='F')

        if use_indices:
          Ju = numpy.hstack(
            numpy.tile(numpy.arange(ny[i]), [ru[i+1], 1]),
            numpy.repeat(Ju, ny[i], axis=0)
          ) # TODO verify correctness
          Ju = Ju[ind, :]

        # Projection from the right -- sample C on U indices
        UC[i-1] = numpy.reshape(coeff[i-1], (rc[i-1]*ny[i-1], rc[i]), order='F')
        UC[i-1] = numpy.matmul(UC[i-1], UC[i])
        UC[i-1] = numpy.reshape(UC[i-1], (rc[i-1], ny[i-1]*ru[i]), order='F')
        UC[i-1] = UC[i-1][:, ind]

        # reductions with Z
        if kickrank > 0:
          # QR and maxvol Z
          crz = numpy.reshape(crz, (rz[i-1], ny[i-1]*rz[i]), order='F')
          crz = numpy.linalg.qr(numpy.transpose(crz))[0]
          rz[i-1] = numpy.shape(crz)[1]
          ind = tt.maxvol.maxvol(crz)
          # Sample C and U
          ZC[i-1] = numpy.reshape(coeff[i-1], (rc[i-1]*ny[i-1], rc[i]), order='F')
          ZC[i-1] = numpy.matmul(ZC[i-1], ZC[i])
          ZC[i-1] = numpy.reshape(ZC[i-1], (rc[i-1], ny[i-1]*rz[i]), order='F')
          ZC[i-1] = ZC[i-1][:, ind]
          ZU[i-1] = numpy.reshape(u[i-1], (ru[i-1]*ny[i-1], ru[i]), order='F')
          ZU[i-1] = numpy.matmul(ZU[i-1], ZU[i])
          ZU[i-1] = numpy.reshape(ZU[i-1], (ru[i-1], ny[i-1]*rz[i]), order='F')
          ZU[i-1] = ZU[i-1][:, ind]

      print(
        f'=als_cross_parametric= swp={swp} ({dir}), i={i}, dx={dx:.3e}, rank =[{ru[i-1]}, {ru[i]}]'
      )

    i += dir

    if dir > 0 and i == d+1 and swp == nswp:
      break # Last block & last sweep
    elif dir > 0 and i == d and swp < nswp:
      # turn at the right end
      print(f'=als_cross_parametric= fwd swp={swp}, max_dx={max_dx:.3e}, max_rank ={max(ru)}')
      dir = -1
      swp += 1
      max_dx = 0
      if use_indices:
        Ju = []
    elif i == 0 and dir < 0:
      # turn at the left end
      dir = 1
      swp += 1
  
  U0 = numpy.reshape(U0, (1, Nxu, ru[0]), order='F')
  u = [U0] + u
  u = tt.tensor.from_list(u)     

  return (u,[time_solve, time_project],funevals)