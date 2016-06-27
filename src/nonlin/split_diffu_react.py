import sys, time
import scitools.std as plt
import scipy.sparse
import scipy.sparse.linalg
import numpy as np

def diffusion_FE(I, a, f, L, dt, F, T, user_action=None):
    """Diffusion solver, Forward Euler method."""
    Nt = int(round(T/float(dt)))
    t = np.linspace(0, Nt*dt, Nt+1)   # Mesh points in time
    dx = np.sqrt(a*dt/F)
    Nx = int(round(L/dx))
    x = np.linspace(0, L, Nx+1)       # Mesh points in space
    # Make sure dx and dt are compatible with x and t
    dx = x[1] - x[0]
    dt = t[1] - t[0]

    u   = np.zeros(Nx+1)
    u_n = np.zeros(Nx+1)

    # Set initial condition u(x,0) = I(x)
    for i in range(0, Nx+1):
        u_n[i] = I(x[i])
    if user_action is not None:
        user_action(u_n, x, t, 0)

    for n in range(0, Nt):
        # Compute u at inner mesh points
        u[1:Nx] = u_n[1:Nx] +  \
                  F*(u_n[0:Nx-1] - 2*u_n[1:Nx] + u_n[2:Nx+1]) + \
                  f(u[1:Nx], t[n])

        # Insert boundary conditions
        u[0] = 0;  u[Nx] = 0

        if user_action is not None:
            user_action(u, x, t, 0)

        u_n, u = u, u_n


def solver_theta(I, a, f, L, dt, F, T, theta=0.5, u_L=0, u_R=0,
                 user_action=None):
    """
    Full solver for the model problem using the theta-rule
    difference approximation in time (no restriction on F,
    i.e., the time step when theta >= 0.5).
    Vectorized implementation and sparse (tridiagonal)
    coefficient matrix.
    """
    import time;  t0 = time.clock()  # for measuring the CPU time

    Nt = int(round(T/float(dt)))
    t = np.linspace(0, Nt*dt, Nt+1)   # Mesh points in time
    dx = np.sqrt(a*dt/F)
    Nx = int(round(L/dx))
    x = np.linspace(0, L, Nx+1)       # Mesh points in space
    # Make sure dx and dt are compatible with x and t
    dx = x[1] - x[0]
    dt = t[1] - t[0]

    u   = np.zeros(Nx+1)   # solution array at t[n+1]
    u_n = np.zeros(Nx+1)   # solution at t[n]

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
        u_n[i] = I(x[i])

    if user_action is not None:
        user_action(u_n, x, t, 0)

    # Time loop
    for n in range(0, Nt):
        b[1:-1] = u_n[1:-1] + \
                  Fr*(u_n[:-2] - 2*u_n[1:-1] + u_n[2:]) + \
                  dt*theta*f(x[1:-1], t[n+1]) + \
                  dt*(1-theta)*f(x[1:-1], t[n])
        b[0] = u_L; b[-1] = u_R  # boundary conditions
        u[:] = scipy.sparse.linalg.spsolve(A, b)

        if user_action is not None:
            user_action(u, x, t, n+1)

        # Update u_n before next step
        u_n, u = u, u_n

    t1 = time.clock()
    return t1-t0

def reaction_FE(I, b, dt, T, user_action=None):
    """Reaction solver, Forward Euler method."""
    # Use linear reaction

    def f(u, t):
        return -b*u

    Nt = int(round(T/float(dt)))
    t = np.linspace(0, Nt*dt, Nt+1)   # Mesh points in time
    u = np.zeros(Nt+1)
    u[0] = I
    for n in range(Nt):
        u[n+1] = u[n] - dt*f(u[n], t[n])
        if user_action is not None:
            user_action(u, t[n+1])

def ordinary_splitting():
    # couple solvers, make it possible to run many smaller steps in ODE solver
    # call diffusion part with f=0

def Strang_splitting():
    # couple solvers
    # call diffusion part with f=0

def exact():
    # Simple call to the diffusion part with f=-lambda u: b*u

# It might happen that FE is superior without splitting, but we could
# try Backward Euler with large time steps and show that it becomes
# better with splitting (? - if the ODE needs much finer steps...).
# An extreme is two steps, BE will be ok and correct at infinity anyway,
# while the reaction part will be very, very coarse in the compound solver,
# hopefully better with splitting and more steps for the ODE.

def test_solvers():
    def u_exact(x, t):
        return x*(L-x)*5*t  # fulfills BC at x=0 and x=L

    def I(x):
        return u_exact(x, 0)

    def f(x, t):
        return 5*x*(L-x) + 10*a*t

    a = 3.5
    L = 1.5
    Nx = 4
    F = 0.5
    # Compute dt from Nx and F
    dx = L/Nx;  dt = F/a*dx**2

    def compare(u, x, t, n):      # user_action function
        """Compare exact and computed solution."""
        u_e = u_exact(x, t[n])
        diff = abs(u_e - u).max()
        tol = 1E-14
        print 'max diff:', diff
        assert diff < tol, 'max diff: %g' % diff

    diffusion_FE(I=I, a=a, L=L, dt=dt, F=F, T=0.2,
                 user_action=compare)

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print """Usage %s function arg1 arg2 arg3 ...""" % sys.argv[0]
        sys.exit(0)
    cmd = '%s(%s)' % (sys.argv[1], ', '.join(sys.argv[2:]))
    print cmd
    eval(cmd)
