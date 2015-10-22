#!/usr/bin/env python
# As v1, but using scipy.sparse.diags instead of spdiags
"""
Functions for solving a 1D diffusion equations of simplest types
(constant coefficient, no source term):

      u_t = a*u_xx on (0,L)

with boundary conditions u=0 on x=0,L, for t in (0,T].
Initial condition: u(x,0)=I(x).

The following naming convention of variables are used.

===== ==========================================================
Name  Description
===== ==========================================================
F     The dimensionless number a*dt/dx**2, which implicitly
      together with dt specifies the mesh in space.
dt    Desired time step. dx is computed from dt and F.
T     The stop time for the simulation.
I     Initial condition (Python function of x).
a     Variable coefficient (constant).
L     Length of the domain ([0,L]).
x     Mesh points in space.
t     Mesh points in time.
n     Index counter in time.
u     Unknown at current/new time level.
u_1   u at the previous time level.
dx    Constant mesh spacing in x.
dt    Constant mesh spacing in t.
===== ==========================================================

user_action is a function of (u, x, t, n), u[i] is the solution at
spatial mesh point x[i] at time t[n], where the calling code
can add visualization, error computations, data analysis,
store solutions, etc.
"""
import numpy as np
import time

def solver_FE_simple(I, a, f, L, dt, F, T):
    """
    Simplest expression of the computational algorithm
    using the Forward Euler method and explicit Python loops.
    f must be a Python function of x and t. If None, a
    default f=0 is used.
    """
    import time
    t0 = time.clock()

    dt = float(dt)                # avoid integer division
    Nt = int(round(T/dt))
    t = np.linspace(0, T, Nt+1)   # mesh points in time
    dx = np.sqrt(a*dt/F)
    Nx = int(round(L/dx))
    x = np.linspace(0, L, Nx+1)   # mesh points in space

    if f is None:
        f = lambda x, t: 0

    u   = np.zeros(Nx+1)
    u_1 = np.zeros(Nx+1)

    # Set initial condition u(x,0) = I(x)
    for i in range(0, Nx+1):
        u_1[i] = I(x[i])

    for n in range(0, Nt):
        # Compute u at inner mesh points
        for i in range(1, Nx):
            u[i] = u_1[i] + F*(u_1[i-1] - 2*u_1[i] + u_1[i+1]) + \
                   dt*f(x[i], t[n])

        # Insert boundary conditions
        u[0] = 0;  u[Nx] = 0

        # Switch variables before next step
        u_1, u = u, u_1

    t1 = time.clock()
    # Return u_1 as u since we set u_1=u above
    return u_1, x, t, t1-t0


def solver_FE(I, a, f, L, dt, F, T,
              user_action=None, version='scalar'):
    """
    Vectorized implementation of solver_FE_simple.
    If version='vectorized', f must be a vectorized
    function of x and t (if f is None, a default version
    f=0 is made).
    """
    import time
    t0 = time.clock()

    dt = float(dt)                # avoid integer division
    Nt = int(round(T/dt))
    t = np.linspace(0, T, Nt+1)   # mesh points in time
    dx = np.sqrt(a*dt/F)
    Nx = int(round(L/dx))
    x = np.linspace(0, L, Nx+1)   # mesh points in space

    if f is None:
        if version == 'scalar':
            f = lambda x, t: 0
        else:
            f = lambda x, t: np.zeros(len(x))

    u   = np.zeros(Nx+1)   # solution array
    u_1 = np.zeros(Nx+1)   # solution at t-dt

    # Set initial condition
    for i in range(0,Nx+1):
        u_1[i] = I(x[i])

    if user_action is not None:
        user_action(u_1, x, t, 0)

    for n in range(0, Nt):
        # Update all inner points
        if version == 'scalar':
            for i in range(1, Nx):
                u[i] = u_1[i] + F*(u_1[i-1] - 2*u_1[i] + u_1[i+1]) + \
                       dt*f(x[i], t[n])

        elif version == 'vectorized':
            u[1:Nx] = u_1[1:Nx] +  \
                      F*(u_1[0:Nx-1] - 2*u_1[1:Nx] + u_1[2:Nx+1]) + \
                      dt*f(x[1:Nx], t[n])
        else:
            raise ValueError('version=%s' % version)

        # Insert boundary conditions
        u[0] = 0;  u[Nx] = 0
        if user_action is not None:
            user_action(u, x, t, n+1)

        # Update u_1 before next step
        #u_1[:]= u
        u_1, u = u, u_1

    t1 = time.clock()
    # Return u_1 as solution since we set u_1=u above
    return u_1, x, t, t1-t0


def solver_BE_simple(I, a, f, L, dt, F, T):
    """
    Simplest expression of the computational algorithm
    for the Backward Euler method, using explicit Python loops
    and a dense matrix format for the coefficient matrix.
    f must be a Python function of x and t. If None, a
    default f=0 is used.
    """
    import time
    t0 = time.clock()

    dt = float(dt)                # avoid integer division
    Nt = int(round(T/dt))
    t = np.linspace(0, T, Nt+1)   # mesh points in time
    dx = np.sqrt(a*dt/F)
    Nx = int(round(L/dx))
    x = np.linspace(0, L, Nx+1)   # mesh points in space

    if f is None:
        f = lambda x, t: 0

    u   = np.zeros(Nx+1)
    u_1 = np.zeros(Nx+1)

    # Data structures for the linear system
    A = np.zeros((Nx+1, Nx+1))
    b = np.zeros(Nx+1)

    for i in range(1, Nx):
        A[i,i-1] = -F
        A[i,i+1] = -F
        A[i,i] = 1 + 2*F
    A[0,0] = A[Nx,Nx] = 1

    # Set initial condition u(x,0) = I(x)
    for i in range(0, Nx+1):
        u_1[i] = I(x[i])

    for n in range(0, Nt):
        # Compute b and solve linear system
        for i in range(1, Nx):
            b[i] = u_1[i] + dt*f(x[i], t[n+1])
        b[0] = b[Nx] = 0
        u[:] = np.linalg.solve(A, b)

        # Update u_1 before next step
        #u_1[:]= u
        u_1, u = u, u_1

    t1 = time.clock()
    # Return u_1 as solution since we set u_1=u above
    return u_1, x, t, t1-t0


import scipy.sparse
import scipy.sparse.linalg

def solver_BE(I, a, f, L, dt, F, T, user_action=None):
    """
    Vectorized implementation of solver_BE_simple using also
    a sparse (tridiagonal) matrix for efficiency.
    f must be a vectorized function of x and t (or None, then
    an f=0 function is used).
    """
    import time
    t0 = time.clock()

    dt = float(dt)                # avoid integer division
    Nt = int(round(T/dt))
    t = np.linspace(0, T, Nt+1)   # mesh points in time
    dx = np.sqrt(a*dt/F)
    Nx = int(round(L/dx))
    x = np.linspace(0, L, Nx+1)   # mesh points in space

    if f is None:
        f = lambda x, t: np.zeros(len(x))

    u   = np.zeros(Nx+1)   # solution array at t[n+1]
    u_1 = np.zeros(Nx+1)   # solution at t[n]

    # Representation of sparse matrix and right-hand side
    diagonal = np.zeros(Nx+1)
    lower    = np.zeros(Nx)
    upper    = np.zeros(Nx)
    b        = np.zeros(Nx+1)

    # Precompute sparse matrix
    diagonal[:] = 1 + 2*F
    lower[:] = -F  #1
    upper[:] = -F  #1
    # Insert boundary conditions
    diagonal[0] = 1
    upper[0] = 0
    diagonal[Nx] = 1
    lower[-1] = 0

    A = scipy.sparse.diags(
        diagonals=[diagonal, lower, upper],
        offsets=[0, -1, 1], shape=(Nx+1, Nx+1),
        format='csr')
    print A.todense()

    # Set initial condition
    for i in range(0,Nx+1):
        u_1[i] = I(x[i])

    if user_action is not None:
        user_action(u_1, x, t, 0)

    for n in range(0, Nt):
        b = u_1 + dt*f(x, t[n+1])
        b[0] = b[-1] = 0.0  # boundary conditions
        u[:] = scipy.sparse.linalg.spsolve(A, b)

        if user_action is not None:
            user_action(u, x, t, n+1)

        # Update u_1 before next step
        #u_1[:] = u
        u_1, u = u, u_1

    t1 = time.clock()
    # Return u_1 as solution since we set u_1=u above
    return u_1, x, t, t1-t0


def solver_theta(I, a, f, L, dt, F, T, theta=0.5, U_0=0, U_L=0,
                 user_action=None):
    """
    Full solver for the model problem using the theta-rule
    difference approximation in time.
    The boundary conditions U_0 and U_L can be numbers, functions
    of time, or None. None means Neumann condition du/dn=0.
    Vectorized implementation and sparse (tridiagonal)
    coefficient matrix.
    """
    import time
    t0 = time.clock()

    dt = float(dt)                # avoid integer division
    Nt = int(round(T/dt))
    t = np.linspace(0, T, Nt+1)   # mesh points in time
    dx = np.sqrt(a*dt/F)
    Nx = int(round(L/dx))
    x = np.linspace(0, L, Nx+1)   # mesh points in space

    if f is None:
        f = lambda x, t: np.zeros(len(x))

    u   = np.zeros(Nx+1)   # solution array at t[n+1]
    u_1 = np.zeros(Nx+1)   # solution at t[n]

    # Make U_0 and U_L functions if they are float/int
    if isinstance(U_0, (float,int)):
        _U_0 = float(U_0)  # make copy of U_0
        U_0 = lambda t: _U_0
    if isinstance(U_L, (float,int)):
        _U_L = float(U_L)  # make copy of U_L
        U_L = lambda t: _U_L

    # Representation of sparse matrix and right-hand side
    diagonal = np.zeros(Nx+1)
    lower    = np.zeros(Nx)
    upper    = np.zeros(Nx)
    b        = np.zeros(Nx+1)

    # Precompute sparse matrix (scipy format)
    Fl = F*theta
    Fr = F*(1-theta)
    diagonal[:] = 1 + 2*Fl
    lower[:] = -Fl  #1
    upper[:] = -Fl  #1
    # Insert boundary conditions
    diagonal[0] = 1
    upper[0] = 0
    diagonal[Nx] = 1
    lower[-1] = 0

    diags = [0, -1, 1]
    A = scipy.sparse.diags(
        diagonals=[diagonal, lower, upper],
        offsets=[0, -1, 1], shape=(Nx+1, Nx+1),
        format='csr')
    #print A.todense()

    # Set initial condition
    for i in range(0,Nx+1):
        u_1[i] = I(x[i])

    if user_action is not None:
        user_action(u_1, x, t, 0)

    # Time loop
    for n in range(0, Nt):
        b[1:-1] = u_1[1:-1] + Fr*(u_1[:-2] - 2*u_1[1:-1] + u_1[2:]) \
                  + theta*dt*f(x[1:-1], t[n+1]) \
                  + (1-theta)*dt*f(x[1:-1], t[n])
        # Assign boundary conditions at this new t[n+1] step
        b[0] = U_0(t[n+1]); b[-1] = U_L(t[n+1])
        # Solve matrix system A*u = b
        u[:] = scipy.sparse.linalg.spsolve(A, b)

        if user_action is not None:
            user_action(u, x, t, n+1)

        # Update u_1 before next step
        u_1, u = u, u_1

    t1 = time.clock()
    # Return u_1 as solution since we set u_1=u above
    return u_1, x, t, t1-t0


def viz(I, a, f, L, dt, F, T, umin, umax,
        solver='solver_FE', animate=True, framefiles=True):

    solutions = []
    from scitools.std import plot, savefig

    def plot_u(u, x, t, n):
        if n == 0:
            # Store x and t first in solutions
            solutions.append(x)
            solutions.append(t)
        solutions.append(u.copy())
        plot(x, u, 'r-', axis=[0, L, umin, umax], title='t=%f' % t[n])
        if framefiles:
            savefig('tmp_frame%04d.png' % n)
            if n in [0, 2, 5, 10, 25, 50, 75, 100, 250, 500]: savefig('tmp_frame%04d.pdf' % n)
        if t[n] == 0:
            time.sleep(2)
        elif not framefiles:
            # It takes time to write files so pause is needed
            # for screen only animation
            time.sleep(0.2)

    user_action = plot_u if animate else lambda u,x,t,n: None

    u, x, t, cpu = eval(solver)(I, a, f, L, dt, F, T,
                                user_action=user_action)
    return solutions, cpu


def plug(solver='solver_FE', F=0.5, dt=0.0002):
    """Plug profile as initial condition."""
    L = 1.
    a = 1
    T = 0.1

    def I(x):
        return 0 if abs(x-L/2.0) > 0.1 else 1

    u, cpu = viz(I, a, None, L, dt, F, T, umin=-0.1, umax=1.1,
                 solver=solver, animate=True, framefiles=True)
    return u

def gaussian(solver='solver_FE', F=0.5, dt=0.0002, sigma=0.1):
    """Gaussian profile as initial condition."""
    L = 1.
    a = 1
    T = 0.1

    def I(x):
        return np.exp(-0.5*((x-L/2.0)**2)/sigma**2)

    u, cpu = viz(I, a, None, L, dt, F, T, umin=-0.1, umax=1.1,
                 solver=solver, animate=True, framefiles=True)
    return u


def expsin(scheme='FE', F=0.5, m=3):
    L = 10.0
    a = 1.0
    T = 1.2

    def exact(x, t):
        return np.exp(-m**2*pi**2*a/L**2*t)*np.sin(m*pi/L*x)

    def I(x):
        return exact(x, 0)

    Nx = 80
    viz(I, a, None, L, Nx, F, T, -1, 1, scheme=scheme,
        animate=True, framefiles=True)

    # Convergence study
    def action(u, x, t, n):
        e = abs(u - exact(x, t[n])).max()
        errors.append(e)

    errors = []
    Nx_values = [10, 20, 40, 80, 160]
    for Nx in Nx_values:
        dx = L/Nx; dt = F/a*dx**2   # find dt from Nx
        eval('solver_'+scheme)(I, a, None, L, dt, F, T, user_action=action)
        print dt, errors[-1]


def expcos(F=0.5, m=3):
    L = 10.0
    a = 1.0
    T = 1.2

    def exact(x, t):
        return exp(-m**2*pi**2*a/L**2*t)*cos(m*pi/L*x)

    def I(x):
        return exact(x, 0)

    # Convergence study
    def action(u, x, t, n):
        e = abs(u - exact(x, t[n])).max()
        errors.append(e)

    errors = []
    Nx_values = [10, 20, 40, 80, 160]
    for Nx in Nx_values:
        dx = L/Nx; dt = F/a*dx**2
        solver_theta(I, a, None, L, dt, F, T, user_action=action,
                     U_0=exact(0,t), U_L=exact(L,t))
        print dt, errors[-1]


def test_solver_FE():
    """Test the FE solvers."""
    import sympy as sym
    x, t, a, L = sym.symbols('x t a L')
    u = x*(L-x)*5*t

    def pde(u):
        return sym.diff(u, t) - a*sym.diff(u, x, x)

    f = sym.simplify(pde(u))
    print 'source', f, 'if u is', u

    # Convert f and u to Python functions. This needs
    # the parameters in u (a and L) to be numbers, and
    # L and a must be substituted by these numbers in the
    # sympy expressions.
    a = 0.5
    L = 1.5
    u_exact = sym.lambdify(
        [x, t], u.subs('L', L).subs('a', a), modules='numpy')
    f = sym.lambdify(
        [x, t], f.subs('L', L).subs('a', a), modules='numpy')
    I = lambda x: u_exact(x, 0)

    dx = L/3  # 3 cells
    F = 0.5
    dt = F*dx**2

    u, x, t, cpu = solver_FE_simple(
        I=I, a=a, f=f, L=L, dt=dt, F=F, T=2)
    u_e = u_exact(x, t[-1])
    diff = abs(u_e - u).max()
    tol = 1E-14
    assert diff < tol, 'max diff solver_FE_simple: %g' % diff

    u, x, t, cpu = solver_FE(
        I=I, a=a, f=f, L=L, dt=dt, F=F, T=2,
        user_action=None, version='scalar')
    u_e = u_exact(x, t[-1])
    diff = abs(u_e - u).max()
    tol = 1E-14
    assert diff < tol, 'max diff solver_FE, scalar: %g' % diff

    u, x, t, cpu = solver_FE(
        I=I, a=a, f=f, L=L, dt=dt, F=F, T=2,
        user_action=None, version='vectorized')
    u_e = u_exact(x, t[-1])
    diff = abs(u_e - u).max()
    tol = 1E-14
    assert diff < tol, 'max diff solver_FE, vectorized: %g' % diff

def test_solvers():
    """Test all solvers in one function."""
    def u_exact(x, t):
        return (L-x)**2*5*t  # fulfills BC at x=0 and x=L

    def I(x):
        return u_exact(x, 0)

    def f(x, t):
        return (L-x)**2*t - a*2*5*t

    a = 3.5
    L = 1.5
    import functools
    s = functools.partial  # object for calling a function w/args
    solvers = [
        s(solver_FE_simple, I=I, a=a, f=f, L=L, Nx=3, F=0.5, T=2),
        s(solver_FE,        I=I, a=a, f=f, L=L, Nx=3, F=0.5, T=2,
          user_action=None, version='scalar'),
        s(solver_FE,        I=I, a=a, f=f, L=L, Nx=3, F=0.5, T=2,
          user_action=None, version='vectorized'),
        s(solver_BE_simple, I=I, a=a, f=f, L=L, Nx=3, F=0.5, T=2),
        s(solver_BE,        I=I, a=a, f=f, L=L, Nx=3, F=0.5, T=2,
          user_action=None),
        s(solver_theta,     I=I, a=a, f=f, L=L, Nx=3, F=0.5, T=2,
          theta=0, u_L=0, u_R=0, user_action=None),
        ]
    for solver in solvers:
        u, x, t, cpu = solver()
        u_e = u_exact(x, t[-1])
        diff = abs(u_e - u).max()
        tol = 1E-14
        assert diff < tol, 'max diff: %g' % diff

if __name__ == '__main__':
    import sys
    if len(sys.argv) < 2:
        print """Usage %s function arg1 arg2 arg3 ...""" % sys.argv[0]
        sys.exit(0)
    cmd = '%s(%s)' % (sys.argv[1], ', '.join(sys.argv[2:]))
    print cmd
    eval(cmd)
