#!/usr/bin/env python
"""
Functions for solving 2D diffusion equations of a simple type
(constant coefficient):

      u_t = a*(u_xx + u_yy) + f(x,t)    on  (0,Lx)x(0,Ly)

with boundary conditions u=0 on x=0,Lx and y=0,Ly for t in (0,T].
Initial condition: u(x,0)=I(x).

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
u_1   u at the previous time level.
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
        f = lambda x, y, t: 0

    u   = np.zeros((Nx+1, Ny+1))      # unknown u at new time level
    u_1 = np.zeros((Nx+1, Ny+1))      # u at the previous time level

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

    # Load initial condition into u_1
    for i in Ix:
        for j in Iy:
            u_1[i,j] = I(x[i], y[j])

    # Two-dim coordinate arrays for vectorized function evaluations
    # in the user_action function
    xv = x[:,np.newaxis]
    yv = y[np.newaxis,:]

    if user_action is not None:
        user_action(u_1, x, xv, y, yv, t, 0)

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
                b[p] = u_1[i,j] + \
                  (1-theta)*(
                  Fx*(u_1[i+1,j] - 2*u_1[i,j] + u_1[i-1,j]) +\
                  Fy*(u_1[i,j+1] - 2*u_1[i,j] + u_1[i,j-1]))\
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

        # Update u_1 before next step
        u_1, u = u, u_1

    t1 = time.clock()
    # Return u_1 as solution since we set u_1=u above
    return t, t1-t0

import scipy.sparse
import scipy.sparse.linalg

def solver_sparse(
    I, a, f, Lx, Ly, Nx, Ny, dt, T, theta=0.5,
    U_0x=0, U_0y=0, U_Lx=0, U_Ly=0, user_action=None):
    """
    Full solver for the model problem using the theta-rule
    difference approximation in time. Sparse matrix with
    dedicated Gaussian elimination algorithm.
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
        f = lambda x, y, t: 0

    u   = np.zeros((Nx+1, Ny+1))    # unknown u at new time level
    u_1 = np.zeros((Nx+1, Ny+1))    # u at the previous time level

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

    # Load initial condition into u_1
    for i in Ix:
        for j in Iy:
            u_1[i,j] = I(x[i], y[j])

    # Two-dim coordinate arrays for vectorized function evaluations
    xv = x[:,np.newaxis]
    yv = y[np.newaxis,:]

    if user_action is not None:
        user_action(u_1, x, xv, y, yv, t, 0)

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
        shape=(N, N), format='csr')
    #print A.todense()   # Check that A is correct

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
                b[p] = u_1[i,j] + \
                  (1-theta)*(
                  Fx*(u_1[i+1,j] - 2*u_1[i,j] + u_1[i-1,j]) +\
                  Fy*(u_1[i,j+1] - 2*u_1[i,j] + u_1[i,j-1]))\
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
            b[m(imin,j):m(imax,j)] = u_1[imin:imax,j] + \
                  (1-theta)*(Fx*(
              u_1[imin+1:imax+1,j] -
            2*u_1[imin:imax,j] +
              u_1[imin-1:imax-1,j]) +
                             Fy*(
              u_1[imin:imax,j+1] -
            2*u_1[imin:imax,j] +
              u_1[imin:imax,j-1])) + \
                theta*dt*f_a_np1[imin:imax,j] + \
              (1-theta)*dt*f_a_n[imin:imax,j]
        j = Ny;  b[m(0,j):m(Nx+1,j)] = U_Ly(t[n+1]) # Boundary

        # Solve matrix system A*c = b
        c = scipy.sparse.linalg.spsolve(A, b)

        # Fill u with vector c
        #for j in Iy:  # vectorize y lines
        #    u[0:Nx+1,j] = c[m(0,j):m(Nx+1,j)]
        u[:,:] = c.reshape(Ny+1,Nx+1).T

        if user_action is not None:
            user_action(u, x, xv, y, yv, t, n+1)

        # Update u_1 before next step
        u_1, u = u, u_1

    t1 = time.clock()
    # Return u_1 as solution since we set u_1=u above
    return t, t1-t0


def solver_jacobi(
    I, a, f, Lx, Ly, Nx, Ny, dt, T, theta=0.5,
    U_0x=0, U_0y=0, U_Lx=0, U_Ly=0, user_action=None,
    version='vectorized', iteration='Jacobi',
    omega=1.0, max_iter=100, tol=1E-4):
    """
    Full solver for the model problem using the theta-rule
    difference approximation in time. Jacobi or SOR iteration.
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
        f = lambda x, y, t: 0

    u   = np.zeros((Nx+1, Ny+1))      # unknown u at new time level
    u_1 = np.zeros((Nx+1, Ny+1))      # u at the previous time level
    u_  = np.zeros((Nx+1, Ny+1))      # most recent approx to u

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

    # Load initial condition into u_1
    for i in Ix:
        for j in Iy:
            u_1[i,j] = I(x[i], y[j])

    # Two-dim coordinate arrays for vectorized function evaluations
    # in the user_action function
    xv = x[:,np.newaxis]
    yv = y[np.newaxis,:]

    if user_action is not None:
        user_action(u_1, x, xv, y, yv, t, 0)

    # Time loop
    import scipy.linalg
    for n in It[0:-1]:
        # Solve linear system by Jacobi iteration at time level n+1
        u_[:,:] = u_1  # Start value
        converged = False
        r = 0
        while not converged:
            if iteration == 'Jacobi':
                u__ = u_
            elif iteration == 'SOR':
                u__ = u
            if version == 'scalar':
                j = 0
                for i in Ix:
                    u[i,j] = U_0y(t[n+1])  # Boundary
                for j in Iy[1:-1]:
                    i = 0;   u[i,j] = U_0x(t[n+1])  # Boundary
                    i = Nx;  u[i,j] = U_Lx(t[n+1])  # Boundary
                    for i in Ix[1:-1]:
                        u_new = 1.0/(1.0 + 2*Fx + 2*Fy)*(theta*(
                            Fx*(u_[i+1,j] + u__[i-1,j]) +
                            Fy*(u_[i,j+1] + u__[i,j-1])) + \
                        u_1[i,j] + \
                        (1-theta)*(Fx*(
                        u_1[i+1,j] - 2*u_1[i,j] + u_1[i-1,j]) +
                          Fy*(
                        u_1[i,j+1] - 2*u_1[i,j] + u_1[i,j-1]))\
                          + theta*dt*f(i*dx,j*dy,(n+1)*dt) + \
                        (1-theta)*dt*f(i*dx,j*dy,n*dt))
                        u[i,j] = omega*u_new + (1-omega)*u_[i,j]
                j = Ny
                for i in Ix:
                    u[i,j] = U_Ly(t[n+1])  # boundary
            elif version == 'vectorized':
                if iteration == 'SOR':
                    raise NotImplementedError
                j = 0;  u[:,j] = U_0y(t[n+1])  # boundary
                i = 0;  u[i,:] = U_0x(t[n+1])  # boundary
                i = Nx; u[i,:] = U_Lx(t[n+1])  # boundary
                j = Ny; u[:,j] = U_Ly(t[n+1])  # boundary
                # Internal points
                f_a_np1 = f(xv, yv, t[n+1])
                f_a_n   = f(xv, yv, t[n])
                u_new = 1.0/(1.0 + 2*Fx + 2*Fy)*(theta*(Fx*(
                  u_[2:,1:-1] + u_[:-2,1:-1]) +
                    Fy*(
                  u_[1:-1,2:] + u_[1:-1,:-2])) +\
                u_1[1:-1,1:-1] + \
                  (1-theta)*(Fx*(
                  u_1[2:,1:-1] - 2*u_1[1:-1,1:-1] + u_1[:-2,1:-1]) +\
                    Fy*(
                  u_1[1:-1,2:] - 2*u_1[1:-1,1:-1] + u_1[1:-1,:-2]))\
                  + theta*dt*f_a_np1[1:-1,1:-1] + \
                  (1-theta)*dt*f_a_n[1:-1,1:-1])
                u[1:-1,1:-1] = omega*u_new + (1-omega)*u_[1:-1,1:-1]
            r += 1
            converged = np.abs(u-u_).max() < tol or r >= max_iter
            #print r, np.abs(u-u_).max(), np.sqrt(dx*dy*np.sum((u-u_)**2))
            u_[:,:] = u

        print 't=%.2f: Jacobi finished in %d iterations' % \
              (t[n+1], r)

        if user_action is not None:
            user_action(u, x, xv, y, yv, t, n+1)

        # Update u_1 before next step
        u_1, u = u, u_1

    t1 = time.clock()
    # Return u_1 as solution since we set u_1=u above
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

    print 'testing dense matrix'
    t, cpu = solver_dense(
        I, a, f, Lx, Ly, Nx, Ny,
        dt, T, theta, user_action=assert_no_error)

    print 'testing sparse matrix'
    t, cpu = solver_sparse(
        I, a, f, Lx, Ly, Nx, Ny,
        dt, T, theta, user_action=assert_no_error)

    return t, cpu

def test_quadratic():
    # For each of the three schemes (theta = 1, 0.5, 0), a series of
    # meshes are tested (Nx > Ny and Nx < Ny)
    for theta in [1, 0.5, 0]:
        for Nx in range(2, 6, 2):
            for Ny in range(2, 6, 2):
                print 'testing for %dx%d mesh' % (Nx, Ny)
                quadratic(theta, Nx, Ny)

def test_Jacobi():
    Nx = 10
    Ny = 10
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
    T = 0.5

    def examine(u, x, xv, y, yv, t, n):
        # Expected error in amplitude
        dx = x[1] - x[0];  dy = y[1] - y[0];  dt = t[1] - t[0]
        Fx = a*dt/dx**2;  Fy = a*dt/dy**2
        kx = np.pi/Lx;  ky = np.pi/Ly
        px = kx*dx/2;  py = ky*dy/2
        A_BE = (1 + 4*Fx*np.sin(px)**2 + 4*Fy*np.sin(py)**2)**(-n)
        A_e  = np.exp(-a*np.pi**2*(Lx**(-2) + Ly**(-2))*t[n])
        A_diff = abs(A_e - A_BE)
        print 'Max u:', u.max(), \
              'error:', abs(u_exact(xv, yv, t[n]).max() - u.max()), A_diff

    solver_jacobi(
        I=I, a=a, f=f, Lx=Lx, Ly=Ly, Nx=Nx, Ny=Ny, dt=dt, T=T, theta=1,
        U_0x=0, U_0y=0, U_Lx=0, U_Ly=0, user_action=examine,
        version='vectorized', omega=1.0, max_iter=100, tol=1E-4)

if __name__ == '__main__':
    #test_quadratic()
    test_Jacobi()
