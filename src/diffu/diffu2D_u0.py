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

from numpy import *

def solver_theta_2D_dense(I, a, f, Lx, Ly, Nx, Ny, dt, T, theta=0.5,
            U_0x=0, U_0y=0, U_Lx=0, U_Ly=0, user_action=None):
    """
    Full solver for the model problem using the theta-rule
    difference approximation in time. Dense matrix. Gaussian solve.
    """
    
    import time; t0 = time.clock()          # for measuring CPU time

    x = linspace(0, Lx, Nx+1)  # mesh points in x dir
    y = linspace(0, Ly, Ny+1)  # mesh points in y dir
    dx = x[1] - x[0]
    dy = y[1] - y[0]

    xv = x[:,newaxis]          # for vectorized function evaluations
    yv = y[newaxis,:]

    dt = float(dt)                # avoid integer division
    Nt = int(round(T/float(dt)))
    t = linspace(0, Nt*dt, Nt+1)    # mesh points in time

    Fx = a*dt/dx**2;  Fy = a*dt/dy**2   # mesh Fourier numbers (each direction)
    
    # Allow f to be None or 0
    if f is None or f == 0:
        f = lambda x, y, t: 0 

    u = zeros((Nx+1)*(Ny+1))        # unknown u at new time level
    u_1 = zeros((Nx+1)*(Ny+1))      # u at the previous time level

    Ix = range(0, Nx+1)
    Iy = range(0, Ny+1)
    It = range(0, Nt+1)

    # Make U_0x, U_0y, U_Lx and U_Ly functions if they are float/int
    if isinstance(U_0x, (float,int)):
        _U_0x = float(U_0x)  # make copy of U_0x
        U_0x = lambda t: _U_0x
    if isinstance(U_0y, (float,int)):
        _U_0y = float(U_0y)  # make copy of U_0y
        U_0y = lambda t: _U_0y
    if isinstance(U_Lx, (float,int)):
        _U_Lx = float(U_Lx)  # make copy of U_Lx
        U_Lx = lambda t: _U_Lx
    if isinstance(U_Ly, (float,int)):
        _U_Ly = float(U_Ly)  # make copy of U_Ly
        U_Ly = lambda t: _U_Ly

    # Load initial condition into u_1
    for i in Ix:
        for j in Iy:
            k = j*(Ix[-1]+1) + i
            u_1[k] = I(x[i], y[j])

    if user_action is not None:
        user_action(u_1, x, xv, y, yv, t, 0)

    # Data structures for the linear system
    A = zeros(((Ix[-1]+1)*(Iy[-1]+1), (Ix[-1]+1)*(Iy[-1]+1)))
    b = zeros((Ix[-1]+1)*(Iy[-1]+1))
    
    # fill in dense matrix A, line by line
    for i in Ix:     # Nx+1 "first" boundary nodes
        A[i, i] = 1
    for j in Iy[1:-1]:  # loop over "groups" of Nx+1 nodes
        A[j*(Ix[-1]+1), j*(Ix[-1]+1)] = 1   # "first" node is a boundary node
        for i in Ix[1:-1]:
            A[j*(Ix[-1]+1)+i, j*(Ix[-1]+1)+i-(Ix[-1]+1)] = - theta*Fy
            A[j*(Ix[-1]+1)+i, j*(Ix[-1]+1)+i-1] = - theta*Fx
            A[j*(Ix[-1]+1)+i, j*(Ix[-1]+1)+i] = 1 + 2*theta*(Fx+Fy)
            A[j*(Ix[-1]+1)+i, j*(Ix[-1]+1)+i+1] = - theta*Fx
            A[j*(Ix[-1]+1)+i, j*(Ix[-1]+1)+i+(Ix[-1]+1)] = - theta*Fy
        A[j*(Ix[-1]+1)+Ix[-1], j*(Ix[-1]+1)+Ix[-1]] = 1     # "last" node is a boundary node
    for i in range(Iy[-1]*(Ix[-1]+1), Iy[-1]*(Ix[-1]+1)+(Ix[-1]+1)):  # Nx+1 "last" boundary nodes
        A[i, i] = 1
    #print A
       
    import scipy.linalg
    
    # Time loop
    for n in It[0:-1]:
        # compute b and solve linear system
        for i in Ix:     # Nx+1 "first" boundary nodes
            b[i] = U_0y(t[n+1])
        for j in Iy[1:-1]:  # loop over "groups" of Nx+1 nodes
            b[j*(Ix[-1]+1)] = U_0x(t[n+1])  # "first" node is a boundary node
            for i in Ix[1:-1]:
                b[j*(Ix[-1]+1)+i] = u_1[j*(Ix[-1]+1)+i] +\
                    (1-theta)*(Fx*(  u_1[j*(Ix[-1]+1)+(i+1)] -\
                                   2*u_1[j*(Ix[-1]+1)+(i)]   +\
                                     u_1[j*(Ix[-1]+1)+(i-1)])   +\
                               Fy*(  u_1[(j+1)*(Ix[-1]+1)+i] -\
                                   2*u_1[j*(Ix[-1]+1)+i]     +\
                                     u_1[(j-1)*(Ix[-1]+1)+i]))  +\
                    theta*dt*f(i*dx,j*dy,(n+1)*dt) + (1-theta)*dt*f(i*dx,j*dy,n*dt)
            b[j*(Ix[-1]+1)+Ix[-1]] = U_Lx(t[n+1])  # "last" node is a boundary node
        for i in range(Iy[-1]*(Ix[-1]+1), Iy[-1]*(Ix[-1]+1)+(Ix[-1]+1)):  # Nx+1 "last" boundary nodes
            b[i] = U_Ly(t[n+1])
        #print b
        # Solve matrix system A*u = b
        u[:] = scipy.linalg.solve(A, b)
        
        if user_action is not None:
            user_action(u, x, xv, y, yv, t, n+1)
        
        # Update u_1 before next step
        u_1, u = u, u_1
    
    t1 = time.clock()
    # Return u_1 as solution since we set u_1=u above
    return t, t1-t0


import scipy.sparse
import scipy.sparse.linalg

def solver_theta_2D_sparse(I, a, f, Lx, Ly, Nx, Ny, dt, T, theta=0.5,
                    U_0x=0, U_0y=0, U_Lx=0, U_Ly=0, user_action=None):
    """
    Full solver for the model problem using the theta-rule
    difference approximation in time. Sparse matrix with dedicated Gaussian 
    solve. 
    """
    
    import time; t0 = time.clock()          # for measuring CPU time

    x = linspace(0, Lx, Nx+1)  # mesh points in x dir
    y = linspace(0, Ly, Ny+1)  # mesh points in y dir
    dx = x[1] - x[0]
    dy = y[1] - y[0]

    xv = x[:,newaxis]          # for vectorized function evaluations
    yv = y[newaxis,:]

    dt = float(dt)                # avoid integer division
    Nt = int(round(T/float(dt)))
    t = linspace(0, Nt*dt, Nt+1)    # mesh points in time

    Fx = a*dt/dx**2;  Fy = a*dt/dy**2   # mesh Fourier numbers (each direction)
    
    # Allow f to be None or 0
    if f is None or f == 0:
        f = lambda x, y, t: 0 

    u = zeros((Nx+1)*(Ny+1))        # unknown u at new time level
    u_1 = zeros((Nx+1)*(Ny+1))      # u at the previous time level

    Ix = range(0, Nx+1)
    Iy = range(0, Ny+1)
    It = range(0, Nt+1)

    # Make U_0x, U_0y, U_Lx and U_Ly functions if they are float/int
    if isinstance(U_0x, (float,int)):
        _U_0x = float(U_0x)  # make copy of U_0x
        U_0x = lambda t: _U_0x
    if isinstance(U_0y, (float,int)):
        _U_0y = float(U_0y)  # make copy of U_0y
        U_0y = lambda t: _U_0y
    if isinstance(U_Lx, (float,int)):
        _U_Lx = float(U_Lx)  # make copy of U_Lx
        U_Lx = lambda t: _U_Lx
    if isinstance(U_Ly, (float,int)):
        _U_Ly = float(U_Ly)  # make copy of U_Ly
        U_Ly = lambda t: _U_Ly

    # Load initial condition into u_1
    for i in Ix:
        for j in Iy:
            k = j*(Ix[-1]+1) + i
            u_1[k] = I(x[i], y[j])

    if user_action is not None:
        user_action(u_1, x, xv, y, yv, t, 0)

    main = zeros((Ix[-1]+1)*(Iy[-1]+1))
    lower = zeros((Ix[-1]+1)*(Iy[-1]+1)-1)
    upper = zeros((Ix[-1]+1)*(Iy[-1]+1)-1)
    lower2 = zeros((Ix[-1]+1)*Iy[-1])
    upper2 = zeros((Ix[-1]+1)*Iy[-1])
    b = zeros((Ix[-1]+1)*(Iy[-1]+1))
       
    # Precompute sparse matrix
    main[0:Ix[-1]+1] = 1                        # Nx+1 "first" boundary nodes
    for j in Iy[1:-1]:         # loop over "groups" of Nx+1 nodes
        main[j*(Ix[-1]+1)] = 1              # "first" node is a boundary node
        lower2[(j-1)*(Ix[-1]+1)+1:(j-1)*(Ix[-1]+1)+Ix[-1]] = - theta*Fy
        lower[j*(Ix[-1]+1):j*(Ix[-1]+1)+Ix[-1]-1] = - theta*Fx
        main[j*(Ix[-1]+1)+1:j*(Ix[-1]+1)+Ix[-1]] = 1 + 2*theta*(Fx+Fy)
        upper[j*(Ix[-1]+1)+1:j*(Ix[-1]+1)+Ix[-1]] = - theta*Fx
        upper2[j*(Ix[-1]+1)+1:j*(Ix[-1]+1)+Ix[-1]] = - theta*Fy
        main[j*(Ix[-1]+1)+Ix[-1]] = 1        # "last" node is a boundary node
    main[Iy[-1]*(Ix[-1]+1):Iy[-1]*(Ix[-1]+1)+(Ix[-1]+1)] = 1        # Nx+1 "last" boundary nodes
        
    A = scipy.sparse.diags(diagonals = [main, lower, upper, lower2, upper2],
        offsets=[0, -1, 1, -(Ix[-1]+1), (Ix[-1]+1)], shape=((Ix[-1]+1)*(Iy[-1]+1), (Ix[-1]+1)*(Iy[-1]+1)),
                 format='csr')   
    #print A.todense()   # Check that A is correct
        
    # Time loop
    for n in It[0:-1]:
        # compute b and solve linear system
        for i in Ix:     # Nx+1 "first" boundary nodes
            b[i] = U_0y(t[n+1])
        for j in Iy[1:-1]:  # loop over "groups" of Nx+1 nodes
            b[j*(Ix[-1]+1)] = U_0x(t[n+1])    # "first" node is a boundary node
            for i in Ix[1:-1]:
                b[j*(Ix[-1]+1)+i] = u_1[j*(Ix[-1]+1)+i] +\
                    (1-theta)*(Fx*(  u_1[j*(Ix[-1]+1)+(i+1)] -\
                                   2*u_1[j*(Ix[-1]+1)+(i)]   +\
                                     u_1[j*(Ix[-1]+1)+(i-1)])   +\
                               Fy*(  u_1[(j+1)*(Ix[-1]+1)+i] -\
                                   2*u_1[j*(Ix[-1]+1)+i]     +\
                                     u_1[(j-1)*(Ix[-1]+1)+i]))  +\
                    theta*dt*f(i*dx,j*dy,(n+1)*dt) + (1-theta)*dt*f(i*dx,j*dy,n*dt)
            b[j*(Ix[-1]+1)+Ix[-1]] = U_Lx(t[n+1])   # "last" node is a boundary node
        for i in range(Iy[-1]*(Ix[-1]+1), Iy[-1]*(Ix[-1]+1)+(Ix[-1]+1)):  # Nx+1 "last" boundary nodes
            b[i] = U_Ly(t[n+1])
        #print b
        # Solve matrix system A*u = b        
        u[:] = scipy.sparse.linalg.spsolve(A, b)
      
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

    a = 3.5
    Lx = Ly = 1.5
    dt = 0.5
    T = 2

    def assert_no_error(u, x, xv, y, yv, t, n):
        u_e = zeros(len(x)*len(y))
        for i in range(0, Nx+1):        
            for j in range(0, Ny+1):
                k = j*(Nx+1) + i
                u_e[k] = u_exact(x[i], y[j], t[n])
         
        # just a preliminary print
        for i in range(0, len(u), 1):
            print 'u_e = %g , u = %g' % (u_e[i], u[i])
               
        diff = abs(u - u_e).max()
        tol = 1E-12
        msg = 'diff=%g, step %d, time=%g' % (diff, n, t[n])
        assert diff < tol, msg

    t, cpu = solver_theta_2D_dense(I, a, f, Lx, Ly, Nx, Ny,
                        dt, T, theta, user_action=assert_no_error)
                                                
    t, cpu = solver_theta_2D_sparse(I, a, f, Lx, Ly, Nx, Ny,
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
    
if __name__ == '__main__':
    test_quadratic()


