#!/usr/bin/env python
"""
Function for solving 3D diffusion equations of a simple type
(constant coefficient) with incomplete LU factorization and 
Conjug. grad. method.

      u_t = a*(u_xx + u_yy + u_zz) + f(x,y,z,t)    on  (0,Lx)x(0,Ly)x(0,Lz)

with boundary conditions u=0 on x=0,Lx and y=0,Ly and z=0,Lz for t in (0,T].
Initial condition: u(x,y,z,0)=I(x,y,z).

The following naming convention of variables are used.

===== ==========================================================
Name  Description
===== ==========================================================
Fx     The dimensionless number a*dt/dx**2, which implicitly
      together with dt specifies the mesh in x.
Fy     The dimensionless number a*dt/dy**2, which implicitly
      together with dt specifies the mesh in y.
Fz     The dimensionless number a*dt/dz**2, which implicitly
      together with dt specifies the mesh in z.
Nx    Number of mesh cells in x direction.
Ny    Number of mesh cells in y direction.
Nz    Number of mesh cells in z direction.
dt    Desired time step. dx is computed from dt and F.
T     The stop time for the simulation.
I     Initial condition (Python function of x and y).
a     Variable coefficient (constant).
Lx     Length of the domain ([0,Lx]).
Ly     Length of the domain ([0,Ly]).
Lz     Length of the domain ([0,Lz]).
x     Mesh points in x.
y     Mesh points in y.
z     Mesh points in z.
t     Mesh points in time.
n     Index counter in time.
u     Unknown at current/new time level.
u_n   u at the previous time level.
dx    Constant mesh spacing in x.
dy    Constant mesh spacing in y.
dz    Constant mesh spacing in z.
dt    Constant mesh spacing in t.
===== ==========================================================

The mesh points are numbered as (0,0,0), (1,0,0), (2,0,0),
..., (Nx,0,0), (0,1,0), (1,1,0), ..., (Nx,1,0), ..., (0,Ny,0), 
(1,Ny,0), ...(Nx,Ny,0), (0,0,1), (1,0,1), ...(Nx,0,1), (0,1,1), 
(1,1,1), ...(Nx, Ny, Nz). 3D-index i,j,k maps to a single index 
s = k*(Nx+1)*(Ny+1) + j*(Nx+1) + i, where i,j,k is the node ID 
and s is the corresponding location in the solution array c when 
solving Ac = b.

f can be specified as None or 0, resulting in f=0.

user_action: function of (u, x, y, z, t, n) called at each time
level (x, y and z are one-dimensional coordinate vectors).
This function allows the calling code to plot the solution,
compute errors, etc.
"""
import sys
import numpy as np
import scipy.sparse
import scipy.linalg
import scipy.sparse.linalg

def solver_sparse_CG(
    I, a, f, Lx, Ly, Lz, Nx, Ny, Nz, dt, T, theta=0.5,
    U_0x=0, U_0y=0, U_0z=0, U_Lx=0, U_Ly=0, U_Lz=0, user_action=None):
    """
    Full solver for the model problem using the theta-rule
    difference approximation in time. Sparse matrix with ILU
    preconditioning and CG solve.
    """
    import time; t0 = time.clock()  # for measuring CPU time

    x = np.linspace(0, Lx, Nx+1)       # mesh points in x dir
    y = np.linspace(0, Ly, Ny+1)       # mesh points in y dir
    z = np.linspace(0, Lz, Nz+1)       # mesh points in z dir
    dx = x[1] - x[0]
    dy = y[1] - y[0]
    dz = z[1] - z[0]

    dt = float(dt)                  # avoid integer division
    Nt = int(round(T/float(dt)))
    t = np.linspace(0, Nt*dt, Nt+1) # mesh points in time

    # Mesh Fourier numbers in each direction
    Fx = a*dt/dx**2
    Fy = a*dt/dy**2
    Fz = a*dt/dz**2

    # Allow f to be None or 0
    if f is None or f == 0:
        f = lambda x, y, z, t: 0

    # unknown u at new time level
    u   = np.zeros((Nx+1, Ny+1, Nz+1)) 
    # u at the previous time level
    u_n = np.zeros((Nx+1, Ny+1, Nz+1)) 

    Ix = range(0, Nx+1)
    Iy = range(0, Ny+1)
    Iz = range(0, Nz+1)
    It = range(0, Nt+1)

    # Make U_0x, U_0y, U_0z, U_Lx, U_Ly, U_Lz 
    # functions if they are float/int
    if isinstance(U_0x, (float,int)):
        _U_0x = float(U_0x)  # make copy of U_0x
        U_0x = lambda t: _U_0x
    if isinstance(U_0y, (float,int)):
        _U_0y = float(U_0y)  # make copy of U_0y
        U_0y = lambda t: _U_0y
    if isinstance(U_0z, (float,int)):
        _U_0z = float(U_0z)  # make copy of U_0z
        U_0z = lambda t: _U_0z
    if isinstance(U_Lx, (float,int)):
        _U_Lx = float(U_Lx)  # make copy of U_Lx
        U_Lx = lambda t: _U_Lx
    if isinstance(U_Ly, (float,int)):
        _U_Ly = float(U_Ly)  # make copy of U_Ly
        U_Ly = lambda t: _U_Ly
    if isinstance(U_Lz, (float,int)):
        _U_Lz = float(U_Lz)  # make copy of U_Lz
        U_Lz = lambda t: _U_Lz

    # Load initial condition into u_n
    for i in Ix:
        for j in Iy:
            for k in Iz:
                u_n[i,j,k] = I(x[i], y[j], z[k])

    # 3D coordinate arrays for vectorized function evaluations
    xv = x[:,np.newaxis,np.newaxis]
    yv = y[np.newaxis,:,np.newaxis]
    zv = z[np.newaxis,np.newaxis,:]

    if user_action is not None:
        user_action(u_n, x, xv, y, yv, z, zv, t, 0)

    N = (Nx+1)*(Ny+1)*(Nz+1)
    main   = np.zeros(N)            	# diagonal
    lower  = np.zeros(N-1)          	# subdiagonal
    upper  = np.zeros(N-1)          	# superdiagonal
    lower2 = np.zeros(N-(Nx+1))     	# lower diagonal
    upper2 = np.zeros(N-(Nx+1))     	# upper diagonal
    lower3 = np.zeros(N-(Nx+1)*(Ny+1))  # lower diagonal
    upper3 = np.zeros(N-(Nx+1)*(Ny+1))  # upper diagonal
    b      = np.zeros(N)            	# right-hand side

    # Precompute sparse matrix
    lower_offset = 1
    lower2_offset = Nx+1
    lower3_offset = (Nx+1)*(Ny+1)

    m = lambda i, j, k: k*(Nx+1)*(Ny+1) + j*(Nx+1) + i
    k = 0; main[m(0,0,k):m(Nx+1,Ny+1,k)] = 1  # k=0 boundary layer
    for k in Iz[1:-1]:         # interior mesh layers k=1,...,Nz-1
        j = 0; main[m(0,j,k):m(Nx+1,j,k)] = 1  # j=0 boundary line
        for j in Iy[1:-1]:        # interior mesh lines j=1,...,Ny-1
            i = 0;   main[m(i,j,k)] = 1  # boundary node
            i = Nx;  main[m(i,j,k)] = 1  # boundary node
            # Interior i points: i=1,...,N_x-1
            lower3[m(1,j,k)-lower3_offset:m(Nx,j,k)-lower3_offset] = - theta*Fz
            lower2[m(1,j,k)-lower2_offset:m(Nx,j,k)-lower2_offset] = - theta*Fy
            lower[m(1,j,k)-lower_offset:m(Nx,j,k)-lower_offset] = - theta*Fx
            main[m(1,j,k):m(Nx,j,k)] = 1 + 2*theta*(Fx+Fy+Fz)
            upper[m(1,j,k):m(Nx,j,k)]  = - theta*Fx
            upper2[m(1,j,k):m(Nx,j,k)] = - theta*Fy
            upper3[m(1,j,k):m(Nx,j,k)] = - theta*Fz
        j = Ny; main[m(0,j,k):m(Nx+1,j,k)] = 1  # boundary line
    k = Nz; main[m(0,0,k):m(Nx+1,Ny+1,k)] = 1  # boundary layer

    A = scipy.sparse.diags(
        diagonals=[main, lower, upper, lower2, upper2, lower3, upper3],
        offsets=[0, -lower_offset, lower_offset,
                 -lower2_offset, lower2_offset,
                 -lower3_offset, lower3_offset],
        shape=(N, N), format='csc')
    #print A.todense()   # Check that A is correct

    # Find preconditioner for A (stays constant the whole time interval)
    A_ilu = scipy.sparse.linalg.spilu(A)   
    M = scipy.sparse.linalg.LinearOperator(shape=(N, N), matvec=A_ilu.solve)

    # Time loop
    c = None		# initialize solution vector (Ac = b)
    for n in It[0:-1]:
 
        # Compute b, scalar version
        '''
        k = 0                 # k=0 boundary layer
        for j in Iy:
            for i in Ix:
                p = m(i,j,k);  b[p] = U_0z(t[n+1])

        for k in Iz[1:-1]:    # interior mesh layers k=1,...,Nz-1
            j = 0           # boundary mesh line
            for i in Ix:
                p = m(i,j,k);  b[p] = U_0y(t[n+1])  

            for j in Iy[1:-1]:   # interior mesh lines j=1,...,Ny-1
                i = 0;  p = m(i,j,k);  b[p] = U_0x(t[n+1])  # boundary node

                for i in Ix[1:-1]:           # interior nodes
                    p = m(i,j,k)                   
                    b[p] = u_n[i,j,k] + \
                      (1-theta)*(
                      Fx*(u_n[i+1,j,k] - 2*u_n[i,j,k] + u_n[i-1,j,k]) +\
                      Fy*(u_n[i,j+1,k] - 2*u_n[i,j,k] + u_n[i,j-1,k]) +
                      Fz*(u_n[i,j,k+1] - 2*u_n[i,j,k] + u_n[i,j,k-1]))\
                        + theta*dt*f(i*dx,j*dy,k*dz,(n+1)*dt) + \
                      (1-theta)*dt*f(i*dx,j*dy,k*dz,n*dt)
                i = Nx;  p = m(i,j,k);  b[p] = U_Lx(t[n+1]) # boundary node

            j = Ny          # boundary mesh line
            for i in Ix:
                p = m(i,j,k);  b[p] = U_Ly(t[n+1])          

        k = Nz                 # k=Nz boundary layer
        for j in Iy:
            for i in Ix:
                p = m(i,j,k);  b[p] = U_Lz(t[n+1])

        #print b
        '''
        # Compute b, vectorized version

        # Precompute f in array so we can make slices
        f_a_np1 = f(xv, yv, zv, t[n+1])
        f_a_n   = f(xv, yv, zv, t[n])

        k = 0; b[m(0,0,k):m(Nx+1,Ny+1,k)] = U_0z(t[n+1])  # k=0 boundary layer
        for k in Iz[1:-1]:         # interior mesh layers k=1,...,Nz-1
            j = 0; b[m(0,j,k):m(Nx+1,j,k)] = U_0y(t[n+1])  # j=0, boundary mesh line
            for j in Iy[1:-1]:        # interior mesh lines j=1,...,Ny-1
                i = 0;   p = m(i,j,k);  b[p] = U_0x(t[n+1]) # boundary node
                i = Nx;  p = m(i,j,k);  b[p] = U_Lx(t[n+1]) # boundary node
                # Interior i points: i=1,...,N_x-1
                imin = Ix[1]
                imax = Ix[-1]  # for slice, max i index is Ix[-1]-1
                b[m(imin,j,k):m(imax,j,k)] = u_n[imin:imax,j,k] + \
                      (1-theta)*(Fx*(
                  u_n[imin+1:imax+1,j,k] -
                2*u_n[imin:imax,j,k] +
                  u_n[imin-1:imax-1,j,k]) +
                                 Fy*(
                  u_n[imin:imax,j+1,k] -
                2*u_n[imin:imax,j,k] +
                  u_n[imin:imax,j-1,k]) + \
                                 Fz*(
                  u_n[imin:imax,j,k+1] -
                2*u_n[imin:imax,j,k] +
                  u_n[imin:imax,j,k-1])) + \
                    theta*dt*f_a_np1[imin:imax,j,k] + \
                  (1-theta)*dt*f_a_n[imin:imax,j,k]
            j = Ny;  b[m(0,j,k):m(Nx+1,j,k)] = U_Ly(t[n+1]) # j=Ny, boundary mesh line
        k = Nz; b[m(0,0,k):m(Nx+1,Ny+1,k)] = U_Lz(t[n+1])  # k=Nz boundary layer
        
        # Solve matrix system A*c = b (use previous sol as start vector x0)
        c, info = scipy.sparse.linalg.cg(A, b, x0=c, tol=1e-14, maxiter=N, M=M)

        if info > 0:
            print 'CG: tolerance not achieved within %d iterations' % info
        elif info < 0:
            print 'CG breakdown'

        # Fill u with vector c
        #for k in Iz:
            #for j in Iy:
            #    u[0:Nx+1,j,k] = c[m(0,j,k):m(Nx+1,j,k)]
        u[:,:,:] = c.reshape(Nz+1,Ny+1,Nx+1).T

        if user_action is not None:
            user_action(u, x, xv, y, yv, z, zv, t, n+1)

        # Update u_n before next step
        u_n, u = u, u_n

    t1 = time.clock()

    return t, t1-t0


def quadratic(theta, Nx, Ny, Nz):
    """Exact discrete solution of the scheme."""

    def u_exact(x, y, z, t):
        return 5*t*x*(Lx-x)*y*(Ly-y)*z*(Lz-z)
    def I(x, y, z):
        return u_exact(x, y, z, 0)
    def f(x, y, z, t):
        return 5*x*(Lx-x)*y*(Ly-y)*z*(Lz-z) + 10*a*t*(
          y*(Ly-y)*z*(Lz-z)+x*(Lx-x)*z*(Lz-z)+x*(Lx-x)*y*(Ly-y))

    # Use rectangular box (cuboid) to detect errors in switching 
    # i, j and k in scheme
    Lx = 0.75; Ly = 1.5; Lz = 2.0
    a = 3.5
    dt = 0.5
    T = 2

    def assert_no_error(u, x, xv, y, yv, z, zv, t, n):
        """Assert zero error at all mesh points."""
        u_e = u_exact(xv, yv, zv, t[n])
        diff = abs(u - u_e).max()
        tol = 1E-12
        msg = 'diff=%g, step %d, time=%g' % (diff, n, t[n])
        print msg
        assert diff < tol, msg

    print 'testing sparse matrix, ILU and CG, theta=%g' % theta
    t, cpu = solver_sparse_CG(
        I, a, f, Lx, Ly, Lz, Nx, Ny, Nz,
        dt, T, theta, user_action=assert_no_error)

    return t, cpu

def test_quadratic():
    # For each of the three schemes (theta = 1, 0.5, 0), a series of
    # meshes are tested, where Nc, c = x,y,z, is successively 
    # the largest and smallest among the three Nc values.
    for theta in [1, 0.5, 0]:
        for Nx in range(2, 6, 2):
            for Ny in range(2, 6, 2):
                for Nz in range(2, 6, 2):
                    print 'testing for %dx%dx%d mesh' % (Nx, Ny, Nz)
                    quadratic(theta, Nx, Ny, Nz)

if __name__ == '__main__':
    test_quadratic()
