#!/usr/bin/env python
"""
Functions for solving 2D diffusion equations of a simple type
(constant coefficient):

      u_t = a*(u_xx + u_yy) + f(x,t)    on  (0,Lx)x(0,Ly)

with boundary conditions u=0 on x=0,Lx and y=0,Ly for t in (0,T].
Initial condition: u(x,y,0)=I(x,y).

The following naming convention of variables are used.

===== ==========================================================
Name  Description
===== ==========================================================
Fx     The dimensionless number a*dt/dx**2, which implicitly
      together with dt specifies the mesh in x.
Fy     The dimensionless number a*dt/dy**2, which implicitly
      together with dt specifies the mesh in y.
Nx    Number of mesh cells in x direction.
Ny    Number of mesh cells in y direction.
dt    Desired time step. dx is computed from dt and F.
T     The stop time for the simulation.
I     Initial condition (Python function of x and y).
a     Variable coefficient (constant).
Lx     Length of the domain ([0,Lx]).
Ly     Length of the domain ([0,Ly]).
x     Mesh points in x.
y     Mesh points in y.
t     Mesh points in time.
n     Index counter in time.
u     Unknown at current/new time level.
u_n   u at the previous time level.
dx    Constant mesh spacing in x.
dy    Constant mesh spacing in y.
dt    Constant mesh spacing in t.
===== ==========================================================

The mesh points are numbered as (0,0), (1,0), (2,0),
..., (Nx,0), (0,1), (1,1), ..., (Nx,1), ..., (0,Ny), (1,Ny), ...(Nx,Ny).
2D-index i,j maps to a single index k = j*(Nx+1) + i, where i,j is the
node ID and k is the corresponding location in the solution array u (or u1).

f can be specified as None or 0, resulting in f=0.

user_action: function of (u, x, y, t, n) called at each time
level (x and y are one-dimensional coordinate vectors).
This function allows the calling code to plot the solution,
compute errors, etc.
"""
import sys
import numpy as np

def solver_dense(
    I, a, f, Lx, Ly, Nx, Ny, dt, T, theta=0.5,
    U_0x=0, U_0y=0, U_Lx=0, U_Ly=0, user_action=None):
    """
    Full solver for the model problem using the theta-rule
    difference approximation in time. Dense matrix. Gaussian solve.
    """
    import time; t0 = time.clock()  # for measuring CPU time

    x = np.linspace(0, Lx, Nx+1)       # mesh points in x dir
    y = np.linspace(0, Ly, Ny+1)       # mesh points in y dir
    dx = x[1] - x[0]
    dy = y[1] - y[0]

    dt = float(dt)                    # avoid integer division
    Nt = int(round(T/float(dt)))
    t = np.linspace(0, Nt*dt, Nt+1)   # mesh points in time

    # Mesh Fourier numbers in each direction
    Fx = a*dt/dx**2
    Fy = a*dt/dy**2

    # Allow f to be None or 0
    if f is None or f == 0:
        f = lambda x, y, t: np.zeros((x.size, y.size)) \
            if isinstance(x, np.ndarray) else 0

    u   = np.zeros((Nx+1, Ny+1))      # unknown u at new time level
    u_n = np.zeros((Nx+1, Ny+1))      # u at the previous time level

    Ix = range(0, Nx+1)
    Iy = range(0, Ny+1)
    It = range(0, Nt+1)

    # Make U_0x, U_0y, U_Lx and U_Ly functions if they are float/int
    if isinstance(U_0x, (float,int)):
        _U_0x = float(U_0x)  # Make copy of U_0x
        U_0x = lambda t: _U_0x
    if isinstance(U_0y, (float,int)):
        _U_0y = float(U_0y)  # Make copy of U_0y
        U_0y = lambda t: _U_0y
    if isinstance(U_Lx, (float,int)):
        _U_Lx = float(U_Lx)  # Make copy of U_Lx
        U_Lx = lambda t: _U_Lx
    if isinstance(U_Ly, (float,int)):
        _U_Ly = float(U_Ly)  # Make copy of U_Ly
        U_Ly = lambda t: _U_Ly

    # Load initial condition into u_n
    for i in Ix:
        for j in Iy:
            u_n[i,j] = I(x[i], y[j])

    # Two-dim coordinate arrays for vectorized function evaluations
    # in the user_action function
    xv = x[:,np.newaxis]
    yv = y[np.newaxis,:]

    if user_action is not None:
        user_action(u_n, x, xv, y, yv, t, 0)

    # Data structures for the linear system
    N = (Nx+1)*(Ny+1)  # no of unknowns
    A = np.zeros((N, N))
    b = np.zeros(N)

    # Fill in dense matrix A, mesh line by line
    m = lambda i, j: j*(Nx+1) + i

    # Equation corresponding to mesh point (i,j) has number
    # j*(Nx+1)+i and will contribute to row j*(Nx+1)+i
    # in the matrix.

    # Equations corresponding to j=0, i=0,1,... (u known)
    j = 0
    for i in Ix:
        p = m(i,j);  A[p, p] = 1
    # Loop over all internal mesh points in y diretion
    # and all mesh points in x direction
    for j in Iy[1:-1]:
        i = 0;  p = m(i,j);  A[p, p] = 1   # boundary
        for i in Ix[1:-1]:                 # interior points
            p = m(i,j)
            A[p, m(i,j-1)] = - theta*Fy
            A[p, m(i-1,j)] = - theta*Fx
            A[p, p]        = 1 + 2*theta*(Fx+Fy)
            A[p, m(i+1,j)] = - theta*Fx
            A[p, m(i,j+1)] = - theta*Fy
        i = Nx;  p = m(i,j);  A[p, p] = 1  # boundary
    # Equations corresponding to j=Ny, i=0,1,... (u known)
    j = Ny
    for i in Ix:
        p = m(i,j);  A[p, p] = 1

    # Time loop
    import scipy.linalg
    for n in It[0:-1]:
        # Compute b
        j = 0
        for i in Ix:
            p = m(i,j);  b[p] = U_0y(t[n+1])  # boundary
        for j in Iy[1:-1]:
            i = 0;  p = p = m(i,j);  b[p] = U_0x(t[n+1])  # boundary
            for i in Ix[1:-1]:
                p = m(i,j)                                # interior
                b[p] = u_n[i,j] + \
                  (1-theta)*(
                  Fx*(u_n[i+1,j] - 2*u_n[i,j] + u_n[i-1,j]) +\
                  Fy*(u_n[i,j+1] - 2*u_n[i,j] + u_n[i,j-1]))\
                    + theta*dt*f(i*dx,j*dy,(n+1)*dt) + \
                  (1-theta)*dt*f(i*dx,j*dy,n*dt)
            i = Nx;  p = m(i,j);  b[p] = U_Lx(t[n+1])     # boundary
        j = Ny
        for i in Ix:
            p = m(i,j);  b[p] = U_Ly(t[n+1])  # boundary
        #print b

        # Solve matrix system A*c = b
        # (the solve function always returns a new object so we
        # do not bother with inserting the solution in-place
        # with c[:] = ...)
        c = scipy.linalg.solve(A, b)

        # Fill u with vector c
        for i in Ix:
            for j in Iy:
                u[i,j] = c[m(i,j)]

        if user_action is not None:
            user_action(u, x, xv, y, yv, t, n+1)

        # Update u_n before next step
        u_n, u = u, u_n

    t1 = time.clock()

    return t, t1-t0

import scipy.sparse
import scipy.sparse.linalg

def solver_sparse(
    I, a, f, Lx, Ly, Nx, Ny, dt, T, theta=0.5,
    U_0x=0, U_0y=0, U_Lx=0, U_Ly=0, user_action=None,
    method='direct', CG_prec='ILU', CG_tol=1E-5):
    """
    Full solver for the model problem using the theta-rule
    difference approximation in time. Sparse matrix with
    dedicated Gaussian elimination algorithm (method='direct')
    or ILU preconditioned Conjugate Gradients (method='CG' with
    tolerance CG_tol and preconditioner CG_prec ('ILU' or None)).
    """
    import time; t0 = time.clock()  # for measuring CPU time

    x = np.linspace(0, Lx, Nx+1)       # mesh points in x dir
    y = np.linspace(0, Ly, Ny+1)       # mesh points in y dir
    dx = x[1] - x[0]
    dy = y[1] - y[0]

    dt = float(dt)                  # avoid integer division
    Nt = int(round(T/float(dt)))
    t = np.linspace(0, Nt*dt, Nt+1) # mesh points in time

    # Mesh Fourier numbers in each direction
    Fx = a*dt/dx**2
    Fy = a*dt/dy**2

    # Allow f to be None or 0
    if f is None or f == 0:
        f = lambda x, y, t: np.zeros((x.size, y.size)) \
            if isinstance(x, np.ndarray) else 0

    u   = np.zeros((Nx+1, Ny+1))    # unknown u at new time level
    u_n = np.zeros((Nx+1, Ny+1))    # u at the previous time level

    Ix = range(0, Nx+1)
    Iy = range(0, Ny+1)
    It = range(0, Nt+1)

    # Make U_0x, U_0y, U_Lx and U_Ly functions if they are float/int
    if isinstance(U_0x, (float,int)):
        _U_0x = float(U_0x)  # Make copy of U_0x
        U_0x = lambda t: _U_0x
    if isinstance(U_0y, (float,int)):
        _U_0y = float(U_0y)  # Make copy of U_0y
        U_0y = lambda t: _U_0y
    if isinstance(U_Lx, (float,int)):
        _U_Lx = float(U_Lx)  # Make copy of U_Lx
        U_Lx = lambda t: _U_Lx
    if isinstance(U_Ly, (float,int)):
        _U_Ly = float(U_Ly)  # Make copy of U_Ly
        U_Ly = lambda t: _U_Ly

    # Load initial condition into u_n
    for i in Ix:
        for j in Iy:
            u_n[i,j] = I(x[i], y[j])

    # Two-dim coordinate arrays for vectorized function evaluations
    xv = x[:,np.newaxis]
    yv = y[np.newaxis,:]

    if user_action is not None:
        user_action(u_n, x, xv, y, yv, t, 0)

    N = (Nx+1)*(Ny+1)
    main   = np.zeros(N)            # diagonal
    lower  = np.zeros(N-1)          # subdiagonal
    upper  = np.zeros(N-1)          # superdiagonal
    lower2 = np.zeros(N-(Nx+1))     # lower diagonal
    upper2 = np.zeros(N-(Nx+1))     # upper diagonal
    b      = np.zeros(N)            # right-hand side

    # Precompute sparse matrix
    lower_offset = 1
    lower2_offset = Nx+1

    m = lambda i, j: j*(Nx+1) + i
    j = 0; main[m(0,j):m(Nx+1,j)] = 1  # j=0 boundary line
    for j in Iy[1:-1]:             # Interior mesh lines j=1,...,Ny-1
        i = 0;   main[m(i,j)] = 1  # Boundary
        i = Nx;  main[m(i,j)] = 1  # Boundary
        # Interior i points: i=1,...,N_x-1
        lower2[m(1,j)-lower2_offset:m(Nx,j)-lower2_offset] = - theta*Fy
        lower[m(1,j)-lower_offset:m(Nx,j)-lower_offset] = - theta*Fx
        main[m(1,j):m(Nx,j)] = 1 + 2*theta*(Fx+Fy)
        upper[m(1,j):m(Nx,j)] = - theta*Fx
        upper2[m(1,j):m(Nx,j)] = - theta*Fy
    j = Ny; main[m(0,j):m(Nx+1,j)] = 1  # Boundary line

    A = scipy.sparse.diags(
        diagonals=[main, lower, upper, lower2, upper2],
        offsets=[0, -lower_offset, lower_offset,
                 -lower2_offset, lower2_offset],
        shape=(N, N), format='csc')
    #print A.todense()   # Check that A is correct

    if method == 'CG':
        if CG_prec == 'ILU':
            # Find ILU preconditioner (constant in time)
            A_ilu = scipy.sparse.linalg.spilu(A)  # SuperLU defaults
            M = scipy.sparse.linalg.LinearOperator(
                shape=(N, N), matvec=A_ilu.solve)
        else:
            M = None
        CG_iter = []  # No of CG iterations at time level n

    # Time loop
    for n in It[0:-1]:
        """
        # Compute b, scalar version
        j = 0
        for i in Ix:
            p = m(i,j);  b[p] = U_0y(t[n+1])          # Boundary
        for j in Iy[1:-1]:
            i = 0;  p = m(i,j);  b[p] = U_0x(t[n+1])  # Boundary
            for i in Ix[1:-1]:
                p = m(i,j)                            # Interior
                b[p] = u_n[i,j] + \
                  (1-theta)*(
                  Fx*(u_n[i+1,j] - 2*u_n[i,j] + u_n[i-1,j]) +\
                  Fy*(u_n[i,j+1] - 2*u_n[i,j] + u_n[i,j-1]))\
                    + theta*dt*f(i*dx,j*dy,(n+1)*dt) + \
                  (1-theta)*dt*f(i*dx,j*dy,n*dt)
            i = Nx;  p = m(i,j);  b[p] = U_Lx(t[n+1]) # Boundary
        j = Ny
        for i in Ix:
            p = m(i,j);  b[p] = U_Ly(t[n+1])          # Boundary
        #print b
        """
        # Compute b, vectorized version

        # Precompute f in array so we can make slices
        f_a_np1 = f(xv, yv, t[n+1])
        f_a_n   = f(xv, yv, t[n])

        j = 0; b[m(0,j):m(Nx+1,j)] = U_0y(t[n+1])     # Boundary
        for j in Iy[1:-1]:
            i = 0;   p = m(i,j);  b[p] = U_0x(t[n+1]) # Boundary
            i = Nx;  p = m(i,j);  b[p] = U_Lx(t[n+1]) # Boundary
            imin = Ix[1]
            imax = Ix[-1]  # for slice, max i index is Ix[-1]-1
            b[m(imin,j):m(imax,j)] = u_n[imin:imax,j] + \
                  (1-theta)*(Fx*(
              u_n[imin+1:imax+1,j] -
            2*u_n[imin:imax,j] +
              u_n[imin-1:imax-1,j]) +
                             Fy*(
              u_n[imin:imax,j+1] -
            2*u_n[imin:imax,j] +
              u_n[imin:imax,j-1])) + \
                theta*dt*f_a_np1[imin:imax,j] + \
              (1-theta)*dt*f_a_n[imin:imax,j]
        j = Ny;  b[m(0,j):m(Nx+1,j)] = U_Ly(t[n+1]) # Boundary

        # Solve matrix system A*c = b
        if method == 'direct':
            c = scipy.sparse.linalg.spsolve(A, b)
        elif method == 'CG':
            x0 = u_n.T.reshape(N)  # Start vector is u_n
            CG_iter.append(0)

            def CG_callback(c_k):
                """Trick to count the no of iterations in CG."""
                CG_iter[-1] += 1

            c, info = scipy.sparse.linalg.cg(
                A, b, x0=x0, tol=CG_tol, maxiter=N, M=M,
                callback=CG_callback)
            '''
            if info > 0:
                print 'CG: tolerance %g not achieved within %d iterations' \
                      % (CG_tol, info)
            elif info < 0:
                print 'CG breakdown'
            else:
                print 'CG converged in %d iterations (tol=%g)' \
                      % (CG_iter[-1], CG_tol)
            '''
        # Fill u with vector c
        #for j in Iy:  # vectorize y lines
        #    u[0:Nx+1,j] = c[m(0,j):m(Nx+1,j)]
        u[:,:] = c.reshape(Ny+1,Nx+1).T

        if user_action is not None:
            user_action(u, x, xv, y, yv, t, n+1)

        # Update u_n before next step
        u_n, u = u, u_n

    t1 = time.clock()

    return t, t1-t0


def omega_optimal(Lx, Ly, Nx, Ny):
    """Return the optimal omega for SOR according to 2D formula."""
    dx = Lx/float(Nx)
    dy = Ly/float(Ny)
    rho = (np.cos(np.pi/Nx) + (dx/dy)**2*np.cos(np.pi/Ny))/\
          (1 + (dx/dy)**2)
    omega = 2.0/(1 + np.sqrt(1-rho**2))
    return omega


def solver_classic_iterative(
    I, a, f, Lx, Ly, Nx, Ny, dt, T, theta=0.5,
    U_0x=0, U_0y=0, U_Lx=0, U_Ly=0, user_action=None,
    version='vectorized', iteration='Jacobi',
    omega=1.0, max_iter=100, tol=1E-4):
    """
    Full solver for the model problem using the theta-rule
    difference approximation in time. Jacobi or SOR iteration.
    (omega='optimal' applies an omega according a formula.)
    """
    import time; t0 = time.clock()     # for measuring CPU time

    x = np.linspace(0, Lx, Nx+1)       # mesh points in x dir
    y = np.linspace(0, Ly, Ny+1)       # mesh points in y dir
    dx = x[1] - x[0]
    dy = y[1] - y[0]

    dt = float(dt)                    # avoid integer division
    Nt = int(round(T/float(dt)))
    t = np.linspace(0, Nt*dt, Nt+1)   # mesh points in time

    # Mesh Fourier numbers in each direction
    Fx = a*dt/dx**2
    Fy = a*dt/dy**2

    # Allow f to be None or 0
    if f is None or f == 0:
        f = lambda x, y, t: np.zeros((x.size, y.size)) \
            if isinstance(x, np.ndarray) else 0

    if version == 'vectorized' and iteration == 'SOR':
        if (Nx % 2) != 0 or (Ny % 2) != 0:
            raise ValueError(
                'Vectorized SOR requires even Nx and Ny (%dx%d)'
                % (Nx, Ny))
    if version == 'SOR':
        if omega == 'optimal':
            omega = omega_optimal(Lx, Ly, Nx, Ny)

    u   = np.zeros((Nx+1, Ny+1))      # unknown u at new time level
    u_n = np.zeros((Nx+1, Ny+1))      # u at the previous time level
    u_  = np.zeros((Nx+1, Ny+1))      # most recent approx to u
    if version == 'vectorized':
        u_new = np.zeros((Nx+1, Ny+1))  # help array

    Ix = range(0, Nx+1)
    Iy = range(0, Ny+1)
    It = range(0, Nt+1)

    # Make U_0x, U_0y, U_Lx and U_Ly functions if they are float/int
    if isinstance(U_0x, (float,int)):
        _U_0x = float(U_0x)  # Make copy of U_0x
        U_0x = lambda t: _U_0x
    if isinstance(U_0y, (float,int)):
        _U_0y = float(U_0y)  # Make copy of U_0y
        U_0y = lambda t: _U_0y
    if isinstance(U_Lx, (float,int)):
        _U_Lx = float(U_Lx)  # Make copy of U_Lx
        U_Lx = lambda t: _U_Lx
    if isinstance(U_Ly, (float,int)):
        _U_Ly = float(U_Ly)  # Make copy of U_Ly
        U_Ly = lambda t: _U_Ly

    # Load initial condition into u_n
    for i in Ix:
        for j in Iy:
            u_n[i,j] = I(x[i], y[j])

    # Two-dim coordinate arrays for vectorized function evaluations
    # in the user_action function
    xv = x[:,np.newaxis]
    yv = y[np.newaxis,:]

    if user_action is not None:
        user_action(u_n, x, xv, y, yv, t, 0)

    # Time loop
    import scipy.linalg
    for n in It[0:-1]:
        # Solve linear system by Jacobi or SOR iteration at time level n+1
        u_[:,:] = u_n  # Start value
        converged = False
        r = 0
        while not converged:
            if version == 'scalar':
                if iteration == 'Jacobi':
                    u__ = u_
                elif iteration == 'SOR':
                    u__ = u
                j = 0
                for i in Ix:
                    u[i,j] = U_0y(t[n+1])  # Boundary
                for j in Iy[1:-1]:
                    i = 0;   u[i,j] = U_0x(t[n+1])  # Boundary
                    i = Nx;  u[i,j] = U_Lx(t[n+1])  # Boundary
                    for i in Ix[1:-1]:
                        u_new = 1.0/(1.0 + 2*theta*(Fx + Fy))*(theta*(
                            Fx*(u_[i+1,j] + u__[i-1,j]) +
                            Fy*(u_[i,j+1] + u__[i,j-1])) + \
                        u_n[i,j] + (1-theta)*(
                          Fx*(
                        u_n[i+1,j] - 2*u_n[i,j] + u_n[i-1,j]) +
                          Fy*(
                        u_n[i,j+1] - 2*u_n[i,j] + u_n[i,j-1]))\
                          + theta*dt*f(i*dx,j*dy,(n+1)*dt) + \
                        (1-theta)*dt*f(i*dx,j*dy,n*dt))
                        u[i,j] = omega*u_new + (1-omega)*u_[i,j]
                j = Ny
                for i in Ix:
                    u[i,j] = U_Ly(t[n+1])  # boundary
            elif version == 'vectorized':
                j = 0;  u[:,j] = U_0y(t[n+1])  # boundary
                i = 0;  u[i,:] = U_0x(t[n+1])  # boundary
                i = Nx; u[i,:] = U_Lx(t[n+1])  # boundary
                j = Ny; u[:,j] = U_Ly(t[n+1])  # boundary
                # Internal points
                f_a_np1 = f(xv, yv, t[n+1])
                f_a_n   = f(xv, yv, t[n])
                def update(u_, u_n, ic, im1, ip1, jc, jm1, jp1):
                    #print '''
#ic:  %s
#im1: %s
#ip1: %s
#jc:  %s
#jm1: %s
#jp1: %s
#''' % (range(u_.shape[0])[ic],range(u_.shape[0])[im1],range(u_.shape[0])[ip1],
#       range(u_.shape[1])[ic],range(u_.shape[1])[im1],range(u_.shape[1])[ip1])
                    return \
                    1.0/(1.0 + 2*theta*(Fx + Fy))*(theta*(
                        Fx*(u_[ip1,jc] + u_[im1,jc]) +
                        Fy*(u_[ic,jp1] + u_[ic,jm1])) +\
                    u_n[ic,jc] + (1-theta)*(
                      Fx*(u_n[ip1,jc] - 2*u_n[ic,jc] + u_n[im1,jc]) +\
                      Fy*(u_n[ic,jp1] - 2*u_n[ic,jc] + u_n[ic,jm1]))+\
                      theta*dt*f_a_np1[ic,jc] + \
                      (1-theta)*dt*f_a_n[ic,jc])

                if iteration == 'Jacobi':
                    ic  = jc  = slice(1,-1)
                    im1 = jm1 = slice(0,-2)
                    ip1 = jp1 = slice(2,None)
                    u_new[ic,jc] = update(
                        u_, u_n, ic, im1, ip1, jc, jm1, jp1)
                    u[ic,jc] = omega*u_new[ic,jc] + (1-omega)*u_[ic,jc]
                elif iteration == 'SOR':
                    u_new[:,:] = u_
                    # Red points
                    ic  = slice(1,-1,2)
                    im1 = slice(0,-2,2)
                    ip1 = slice(2,None,2)
                    jc  = slice(1,-1,2)
                    jm1 = slice(0,-2,2)
                    jp1 = slice(2,None,2)
                    u_new[ic,jc] = update(
                        u_new, u_n, ic, im1, ip1, jc, jm1, jp1)

                    ic  = slice(2,-1,2)
                    im1 = slice(1,-2,2)
                    ip1 = slice(3,None,2)
                    jc  = slice(2,-1,2)
                    jm1 = slice(1,-2,2)
                    jp1 = slice(3,None,2)
                    u_new[ic,jc] = update(
                        u_new, u_n, ic, im1, ip1, jc, jm1, jp1)

                    # Black points
                    ic  = slice(2,-1,2)
                    im1 = slice(1,-2,2)
                    ip1 = slice(3,None,2)
                    jc  = slice(1,-1,2)
                    jm1 = slice(0,-2,2)
                    jp1 = slice(2,None,2)
                    u_new[ic,jc] = update(
                        u_new, u_n, ic, im1, ip1, jc, jm1, jp1)

                    ic  = slice(1,-1,2)
                    im1 = slice(0,-2,2)
                    ip1 = slice(2,None,2)
                    jc  = slice(2,-1,2)
                    jm1 = slice(1,-2,2)
                    jp1 = slice(3,None,2)
                    u_new[ic,jc] = update(
                        u_new, u_n, ic, im1, ip1, jc, jm1, jp1)

                    # Relax
                    c = slice(1,-1)
                    u[c,c] = omega*u_new[c,c] + (1-omega)*u_[c,c]

            r += 1
            converged = np.abs(u-u_).max() < tol or r >= max_iter
            #print r, np.abs(u-u_).max(), np.sqrt(dx*dy*np.sum((u-u_)**2))
            u_[:,:] = u

        print 't=%.2f: %s %s (omega=%g) finished in %d iterations' % \
              (t[n+1], version, iteration, omega, r)

        if user_action is not None:
            user_action(u, x, xv, y, yv, t, n+1)

        # Update u_n before next step
        u_n, u = u, u_n

    t1 = time.clock()

    return t, t1-t0

def quadratic(theta, Nx, Ny):
    """Exact discrete solution of the scheme."""

    def u_exact(x, y, t):
        return 5*t*x*(Lx-x)*y*(Ly-y)
    def I(x, y):
        return u_exact(x, y, 0)
    def f(x, y, t):
        return 5*x*(Lx-x)*y*(Ly-y) + 10*a*t*(y*(Ly-y)+x*(Lx-x))

    # Use rectangle to detect errors in switching i and j in scheme
    Lx = 0.75
    Ly = 1.5
    a = 3.5
    dt = 0.5
    T = 2

    def assert_no_error(u, x, xv, y, yv, t, n):
        """Assert zero error at all mesh points."""
        u_e = u_exact(xv, yv, t[n])
        diff = abs(u - u_e).max()
        tol = 1E-12
        msg = 'diff=%g, step %d, time=%g' % (diff, n, t[n])
        print msg
        assert diff < tol, msg

    print '\ntesting dense matrix'
    t, cpu = solver_dense(
        I, a, f, Lx, Ly, Nx, Ny,
        dt, T, theta, user_action=assert_no_error)

    print '\ntesting sparse matrix'
    t, cpu = solver_sparse(
        I, a, f, Lx, Ly, Nx, Ny,
        dt, T, theta, user_action=assert_no_error,
        method='direct')

    def assert_small_error(u, x, xv, y, yv, t, n):
        """Assert small error at all mesh points for iterative methods."""
        u_e = u_exact(xv, yv, t[n])
        diff = abs(u - u_e).max()
        tol = 1E-12
        tol = 1E-4
        msg = 'diff=%g, step %d, time=%g' % (diff, n, t[n])
        print msg
        assert diff < tol, msg

    tol = 1E-5  # Tolerance in iterative methods
    for iteration in 'Jacobi', 'SOR':
        for version in 'scalar', 'vectorized':
            for theta in 1, 0.5:
                print '\ntesting %s, %s version, theta=%g, tol=%g' % \
                      (iteration, version, theta, tol)
                t, cpu = solver_classic_iterative(
                    I=I, a=a, f=f, Lx=Lx, Ly=Ly, Nx=Nx, Ny=Ny,
                    dt=dt, T=T, theta=theta,
                    U_0x=0, U_0y=0, U_Lx=0, U_Ly=0,
                    user_action=assert_small_error,
                    version=version, iteration=iteration,
                    omega=1.0, max_iter=100, tol=tol)

    print '\ntesting CG+ILU, theta=%g, tol=%g' % (theta, tol)
    solver_sparse(
        I, a, f, Lx, Ly, Nx, Ny, dt, T, theta=0.5,
        user_action=assert_small_error,
        method='CG', CG_prec='ILU', CG_tol=tol)

    print '\ntesting CG, theta=%g, tol=%g' % (theta, tol)
    solver_sparse(
        I, a, f, Lx, Ly, Nx, Ny, dt, T, theta=0.5,
        user_action=assert_small_error,
        method='CG', CG_prec=None, CG_tol=tol)

    return t, cpu

def test_quadratic():
    # For each of the three schemes (theta = 1, 0.5, 0), a series of
    # meshes are tested (Nx > Ny and Nx < Ny)
    for theta in [1, 0.5, 0]:
        for Nx in range(2, 6, 2):
            for Ny in range(2, 6, 2):
                print '\n*** testing for %dx%d mesh' % (Nx, Ny)
                quadratic(theta, Nx, Ny)

def demo_classic_iterative(
    tol=1E-4, iteration='Jacobi',
    version='vectorized', theta=0.5,
    Nx=10, Ny=10):
    Lx = 2.0
    Ly = 1.0
    a = 1.5

    u_exact = lambda x, y, t: \
              np.exp(-a*np.pi**2*(Lx**(-2) + Ly**(-2))*t)*\
              np.sin(np.pi*x/Lx)*np.sin(np.pi*y/Ly)
    I = lambda x, y: u_exact(x, y, 0)
    f = lambda x, y, t: 0 if isinstance(x, (float,int)) else \
        np.zeros((Nx+1,Ny+1))
    dt = 0.2
    dt = 0.05
    T = 0.5

    def examine(u, x, xv, y, yv, t, n):
        # Expected error in amplitude
        dx = x[1] - x[0];  dy = y[1] - y[0];  dt = t[1] - t[0]
        Fx = a*dt/dx**2;  Fy = a*dt/dy**2
        kx = np.pi/Lx;    ky = np.pi/Ly
        px = kx*dx/2;     py = ky*dy/2
        if theta == 1:
            A_d = (1 + 4*Fx*np.sin(px)**2 + 4*Fy*np.sin(py)**2)**(-n)
        else:
            A_d = ((1 - 2*Fx*np.sin(px)**2 - 2*Fy*np.sin(py)**2)/\
                   (1 + 2*Fx*np.sin(px)**2 + 2*Fy*np.sin(py)**2))**n
        A_e  = np.exp(-a*np.pi**2*(Lx**(-2) + Ly**(-2))*t[n])
        A_diff = abs(A_e - A_d)
        u_diff = abs(u_exact(xv, yv, t[n]).max() - u.max())
        print 'Max u: %.2E' % u.max(), \
              'error in u: %.2E' % u_diff, 'ampl.: %.2E' % A_diff, \
              'iter: %.2E' % abs(u_diff - A_diff)

    solver_classic_iterative(
        I=I, a=a, f=f, Lx=Lx, Ly=Ly, Nx=Nx, Ny=Ny,
        dt=dt, T=T, theta=theta,
        U_0x=0, U_0y=0, U_Lx=0, U_Ly=0, user_action=examine,
        #version='vectorized', iteration='Jacobi',
        version=version, iteration=iteration,
        omega=1.0, max_iter=300, tol=tol)


def convergence_rates(theta, num_experiments=6):
    """
    Compute convergence rates. The error measure is the L2 norm
    of the error mesh function in space and time:
    E^2 = dt*sum(dx*dy*sum((u_e - u)**2)).
    """
    # ...for FE and BE, error truncation analysis suggests
    # that the error E goes like C*h**1, where h = dt = dx**2
    # (note that dx = dy is chosen).
    # ...for CN, we similarly have that E goes like
    # C*h**2, where h = dt**2 = dx**2 in this case.

    r_FE_BE_expected = 1
    r_CN_expected = 2
    E_values  = []
    h_values = []
    Lx = 2.0; Ly = 2.0    # Use square domain
    a = 3.5
    T = 1
    kx = np.pi/Lx
    ky = np.pi/Ly
    #p = a*(kx**2 + ky**2)  # f=0, slower approach to asymptotic range
    p = 1

    def u_exact(x, y, t):
        return np.exp(-p*t)*np.sin(kx*x)*np.sin(ky*y)
    def I(x, y):
        return u_exact(x, y, 0)
    def f(x, y, t):
        return (-p + a*(kx**2 + ky**2))*np.exp(-p*t)*np.sin(kx*x)*np.sin(ky*y)

    def add_error_contribution(u, x, xv, y, yv, t, n):
        u_e = u_exact(xv, yv, t[n])
        E2_sum['err'] += np.sum((u_e - u)**2)
        if n == 0:
            # Store away dx, dy, dt in the dict for later use
            E2_sum['dx'] = x[1] - x[0]
            E2_sum['dy'] = y[1] - y[0]
            E2_sum['dt'] = t[1] - t[0]

    def compute_E(h):
        dx = E2_sum['dx']
        dy = E2_sum['dy']
        dt = E2_sum['dt']
        sum_xyt = E2_sum['err']
        E = np.sqrt(dt*dx*dy*sum_xyt)  # L2 norm of error mesh func
        E_values.append(E)
        h_values.append(h)

    def assert_conv_rates():
        r = [np.log(E_values[i+1]/E_values[i])/
             np.log(h_values[i+1]/h_values[i])
             for i in range(0, num_experiments-2, 1)]
        tol = 0.5
        if theta == 0:  # i.e., FE
            diff = abs(r_FE_BE_expected - r[-1])
            msg = 'Forward Euler. r = 1 expected, got=%g' % r[-1]
        elif theta == 1:  # i.e., BE
            diff = abs(r_FE_BE_expected - r[-1])
            msg = 'Backward Euler. r = 1 expected, got=%g' % r[-1]
        else:  # theta == 0.5, i.e, CN
            diff = abs(r_CN_expected - r[-1])
            msg = 'Crank-Nicolson. r = 2 expected, got=%g' % r[-1]
            #print msg
        print 'theta: %g' % theta
        print 'r: ', r
        assert diff < tol, msg

    tol = 1E-5  # Tolerance in iterative methods
    for method in 'direct', 'CG':
        print '\ntesting convergence rate, theta=%g, method=%s' % (theta, method)
        for i in range(num_experiments):
            # Want to do E2_sum += ... in local functions (closures), but
            # a stndard variable E2_sum is reported undefined. Trick: use
            # a dictionary or list instead.
            E2_sum = {'err' : 0}

            N = 2**(i+1)
            Nx = N; Ny = N   # We want dx=dy, so Nx=Ny if Lx=Ly
            dx = float(Lx)/N
            # Find single discretization parameter h = dt and its relation
            # to dt, dx and dy (E=C*h^r)
            if theta == 1:
                dt = h = dx**2
            elif theta == 0:
                h = dx**2
                K = 1./(4*a)
                dt = K*h
            elif theta == 0.5:
                dt = h = dx
            if method == 'direct':
                solver_sparse(
                    I, a, f, Lx, Ly, Nx, Ny, dt, T, theta=0.5,
                    user_action=add_error_contribution,
                    method=method)
            else:
                solver_sparse(
                    I, a, f, Lx, Ly, Nx, Ny, dt, T, theta=0.5,
                    user_action=add_error_contribution,
                    method=method, CG_prec='ILU', CG_tol=tol)
            compute_E(h)
            print 'Experiment no:%d, %d unknowns' % (i+1, (N+1)**2), E_values[-1]
        assert_conv_rates()


def convergence_rates0(theta, num_experiments=10):
    # Sveins version
    """
    Compute convergence rates. The error measure is the L2 norm
    of the error mesh function in space and time:
    E^2 = dt*sum(dx*dy*sum((u_e - u)**2)).
    """
    # ...for FE and BE, error truncation analysis suggests
    # that the error E goes like C*h**1, where h = dt = dx**2
    # (note that dx = dy is chosen).
    # ...for CN, we similarly have that E goes like
    # C*h**2, where h = dt**2 = dx**2 in this case.

    r_FE_BE_expected = 1
    r_CN_expected = 2
    E_values  = []
    dt_values = []
    Lx = 2.0; Ly = 2.0    # Use square domain
    a = 3.5
    T = 2
    p = 1  # hpl: easier to let p match q and s, so f=0
    q = 2*np.pi/Lx
    s = 2*np.pi/Ly  # parameters for exact solution, loop later...?

    # hpl: easier to have p,q,s as "global" variables in the function,
    # they are in I and f anyway...
    def u_exact(p, q, s, x, y, t):
        return np.exp(-p*t)*np.sin(q*x)*np.sin(s*y)
    def I(x, y):
        return u_exact(p, q, s, x, y, 0)
    def f(x, y, t):
        return np.exp(-p*t)*(-p + a*(q**2 + s**2))*np.sin(s*y)*np.sin(q*x)

    # Define boundary conditions (functions of space and time)
    def U_0x(t):
        return 0
    def U_0y(t):
        return 0
    def U_Lx(t):
        return 0
    def U_Ly(t):
        return 0

    def assert_correct_convergence_rate(u, x, xv, y, yv, t, n):
        u_e = u_exact(p, q, s, xv, yv, t[n])
        E2_sum['err'] += np.sum((u_e - u)**2)

        if t[n] == T:  # hpl: dangerous comparison...
            dx = x[1] - x[0]
            dt = t[1] - t[0]
            E = np.sqrt(dt*dx*E2_sum['err'])  # error, 1 simulation, t = [0,T]
            E_values.append(E)
            dt_values.append(dt)
            if counter['i'] == num_experiments:  # i.e., all num. exp. finished
                print '...all experiments finished'
                r = [np.log(E_values[i+1]/E_values[i])/
                     np.log(dt_values[i+1]/dt_values[i])
                     for i in range(0, num_experiments-2, 1)]
                tol = 0.5
                if theta == 0:  # i.e., FE
                    diff = abs(r_FE_BE_expected - r[-1])
                    msg = 'Forward Euler. r = 1 expected, got=%g' % r[-1]
                elif theta == 1:  # i.e., BE
                    diff = abs(r_FE_BE_expected - r[-1])
                    msg = 'Backward Euler. r = 1 expected, got=%g' % r[-1]
                else:  # theta == 0.5, i.e, CN
                    diff = abs(r_CN_expected - r[-1])
                    msg = 'Crank-Nicolson. r = 2 expected, got=%g' % r[-1]
                #print msg
                print 'theta: %g' % theta
                print 'r: ', r
                assert diff < tol, msg

    print '\ntesting convergence rate, sparse matrix, CG, ILU'
    tol = 1E-5  # Tolerance in iterative methods
    counter = {'i' : 0}   # initialize
    for i in range(num_experiments):
        # Want to do E2_sum += ... in local functions (closures), but
        # a stndard variable E2_sum is reported undefined. Trick: use
        # a dictionary or list instead.
        E2_sum = {'err' : 0}
        counter['i'] += 1
        N = 2**(i+1)
        Nx = N; Ny = N
        if theta == 0 or theta == 1:  # i.e., FE or BE
            dt = (float(Lx)/N)**2	# i.e., choose dt = dx**2
        else:			      # theta == 0.5, i.e., CN
            dt = float(Lx)/N            # i.e., choose dt = dx
        solver_sparse(
            I, a, f, Lx, Ly, Nx, Ny, dt, T, theta=0.5,
            user_action=assert_correct_convergence_rate,
            method='CG', CG_prec='ILU', CG_tol=tol)
        print 'Experiment no:%d' % (i+1)


def test_convergence_rate():
    # For each solver, we find the conv rate on a square domain:
    # For each of the three schemes (theta = 1, 0.5, 0), a series of
    # finer and finer meshes are tested, computing the conv.rate r
    # in each case to find the limiting value when dt and dx --> 0.

    for theta in [1, 0.5, 0]:
        convergence_rates(theta, num_experiments=6)


def efficiency():
    """Measure the efficiency of iterative methods."""
    # JUST A SKETCH!
    cpu = {}

    # Find a more advanced example and use large Nx, Ny
    def u_exact(x, y, t):
        return 5*t*x*(Lx-x)*y*(Ly-y)
    def I(x, y):
        return u_exact(x, y, 0)
    def f(x, y, t):
        return 5*x*(Lx-x)*y*(Ly-y) + 10*a*t*(y*(Ly-y)+x*(Lx-x))

    # Use rectangle to detect errors in switching i and j in scheme
    Lx = 0.75
    Ly = 1.5
    a = 3.5
    dt = 0.5
    T = 2

    print '\ntesting sparse matrix LU solver'
    t, cpu = solver_sparse(
        I, a, f, Lx, Ly, Nx, Ny,
        dt, T, theta, user_action=None,
        method='direct')

    theta = 0.5
    tol = 1E-5  # Tolerance in iterative methods
    # Testing Jacobi and Gauss-Seidel
    for iteration in 'Jacobi', 'SOR':
        for version in 'scalar', 'vectorized':
            print '\ntesting %s, %s version, theta=%g, tol=%g' % \
                  (iteration, version, theta, tol)
            t, cpu_ = solver_classic_iterative(
                I=I, a=a, f=f, Lx=Lx, Ly=Ly, Nx=Nx, Ny=Ny,
                dt=dt, T=T, theta=theta,
                U_0x=0, U_0y=0, U_Lx=0, U_Ly=0,
                user_action=None,
                version=version, iteration=iteration,
                omega=1.0, max_iter=100, tol=tol)
            cpu[iteration+'_'+version] = cpu_

    for omega in 'optimal', 1.2, 1.5:
        print '\ntesting SOR, omega:', omega
        t, cpu_ = solver_classic_iterative(
            I=I, a=a, f=f, Lx=Lx, Ly=Ly, Nx=Nx, Ny=Ny,
            dt=dt, T=T, theta=theta,
            U_0x=0, U_0y=0, U_Lx=0, U_Ly=0,
            user_action=None,
            version=version, iteration=iteration,
            omega=1.0, max_iter=100, tol=tol)
        cpu['SOR(omega=%g)' % omega] = cpu_

    print '\ntesting CG+ILU, theta=%g, tol=%g' % (theta, tol)
    t, cpu_ = solver_sparse(
        I, a, f, Lx, Ly, Nx, Ny, dt, T, theta=0.5,
        user_action=None,
        method='CG', CG_prec='ILU', CG_tol=tol)
    cpu['CG+ILU'] = cpu_

    print '\ntesting CG, theta=%g, tol=%g' % (theta, tol)
    t, cpu_ = solver_sparse(
        I, a, f, Lx, Ly, Nx, Ny, dt, T, theta=0.5,
        user_action=None,
        method='CG', CG_prec=None, CG_tol=tol)
    cpu['CG'] = cpu_

    return t, cpu

if __name__ == '__main__':
    #test_quadratic()
    #demo_classic_iterative(
    #    iteration='Jacobi', theta=0.5, tol=1E-4, Nx=20, Ny=20)
    test_convergence_rate()
