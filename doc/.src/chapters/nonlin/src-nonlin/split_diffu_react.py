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

    # Allow f to be None or 0
    if f is None or f == 0:
        f = lambda x, t: np.zeros((x.size)) \
            if isinstance(x, np.ndarray) else 0

    # Set initial condition  u = I or u(x,0) = I(x)
    if isinstance(I, np.ndarray):   # I is an array
        u_n = np.copy(I)
    else:                           # I is a function
        for i in range(0, Nx+1):
            u_n[i] = I(x[i])
        
    if user_action is not None:
        user_action(u_n, x, t, 0)

    for n in range(0, Nt):
        # Compute u at inner mesh points at t_n+1
        u[1:Nx] = u_n[1:Nx] +  \
                  F*(u_n[0:Nx-1] - 2*u_n[1:Nx] + u_n[2:Nx+1]) + \
                  dt*f(u[1:Nx], t[n])

        # Insert boundary conditions
        u[0] = 0;  u[Nx] = 0

        if user_action is not None:
            user_action(u, x, t, n+1)

        u_n, u = u, u_n
    # return u_n, i.e. u (because of ref. swapping)
    return u_n


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

def reaction_FE(I, a, b, L, Nx, local_dt, F, local_T, user_action=None):
    """Reaction solver, Forward Euler method."""

    def f(u, t):
        return -b*u

    local_Nt = int(round(local_T/float(local_dt)))    
    x = np.linspace(0, L, Nx+1)       # Mesh points in space    
    t = np.linspace(0, local_Nt*local_dt, local_Nt+1)  # Mesh points in time
    u = np.zeros(Nx+1)
    u = np.copy(I)      # u is array with values distr in space
    for n in range(local_Nt):
        u = u + local_dt*f(u, t[n])        
        if user_action is not None:
            user_action(u, x, t, n+1)
    return u


def ordinary_splitting(I, a, b, L, dt, dt_Rfactor, F, T, user_action=None):
    '''1st order scheme, i.e. Forward Euler is enough for both
    the diffusion and the reaction part. The time step dt is given for the
    diffusion step, while the time step for the reaction part is
    found as dt/dt_Rfactor, where dt_Rfactor >= 1.
    '''
    Nt = int(round(T/float(dt)))
    dx = np.sqrt(a*dt/F)
    Nx = int(round(L/dx))
    x = np.linspace(0, L, Nx+1)       # Mesh points in space
    u = np.zeros(Nx+1)
    
    # Set initial condition u(x,0) = I(x)
    for i in range(0, Nx+1):
        u[i] = I(x[i])

    # In the following loop, each time step is "covered twice", 
    # first for diffusion, then for reaction
    for n in range(0, Nt):    
        # Diffusion step (one time step dt)
        u_s = diffusion_FE(I=u, a=a, f=0, L=L, dt=dt, F=F, T=dt,
                           user_action=None)
        # Reaction step (potentially many steps within dt)
        u = reaction_FE(I=u_s, a=a, b=b, L=L, Nx=Nx, local_dt=dt/dt_Rfactor, 
                        F=F, local_T=dt, user_action=user_action)
                        

def Strang_splitting(I, a, b, L, dt, dt_Rfactor, F, T, user_action=None):
    # couple solvers
    # call diffusion part with f=0
    #S: To achieve 2nd order, we should use 2nd order schemes?!
    #S: first diff (in loop), then react (in nested loop)   
    print 'hello from Strang_splitting'

def exact():
    # Simple call to the diffusion part with f=-lambda u: b*u
    print 'hello from exact'

# It might happen that FE is superior without splitting, but we could
# try Backward Euler with large time steps and show that it becomes
# better with splitting (? - if the ODE needs much finer steps...).
# An extreme is two steps, BE will be ok and correct at infinity anyway,
# while the reaction part will be very, very coarse in the compound solver,
# hopefully better with splitting and more steps for the ODE.

def test_solvers():
    
    def u_exact(x, t):
        '''exact sol. to: du/dt = a*d^2u/dx^2 - b*u'''
        return np.exp(-(a*k**2 + b)*t) * np.sin(k*x)

    def I(x):
        return u_exact(x, 0)

    def f(u, t):
        return -b*u  

    a = 3.5
    b = 1
    L = 1.5
    k = np.pi/L
    Nx = 4
    F = 0.5
    # Compute dt from Nx and F
    dx = L/Nx;  dt = F/a*dx**2
    
    def compare(u, x, t, n):      # user_action function
        """Compare exact and (approximate) computed solution."""
        print 't:', t[n]
        print 'u:', u
        u_e = u_exact(u, t[n])
        print 'u_e:', u_e
        diff = abs(u_e - u).max()
        print 'max diff:', diff

    #diffusion_FE(I=I, a=a, f=f, L=L, dt=dt, F=F, T=0.2,
    #             user_action=compare)
    ordinary_splitting(I=I, a=a, b=b, L=L, dt=dt, 
                       dt_Rfactor=10, F=F, T=0.2,
                       user_action=compare)

def convergence_rates(
    u_exact,                 # Python function for exact solution
    I, f, a, b, L,           # physical parameters
    dx, dt_Rfactor, dt0, num_meshes, T):  # numerical parameters
    """
    Halve the time step and estimate convergence rates 
    for num_meshes simulations.
    """
    # First define an appropriate user action function
    global error
    error = 0  # error computed in the user action function

    def compute_error(u, x, t, n):
        global error  # must be global to be altered here
        # (otherwise error is a local variable, different
        # from error defined in the parent function)
        if n == 0:
            error = 0
        else:
            print 'u:', u
            print 'u_exact:', u_exact(x, t[n])
            error = max(error, np.abs(u - u_exact(x, t[n])).max())

    # Run finer and finer resolutions and compute true errors
    E = []
    h = []  # dt
    dt = dt0
    for i in range(num_meshes):
        F = a*dt/dx**2
        #ordinary_splitting(I=I, a=a, b=b, L=L, dt=dt, 
        #               dt_Rfactor=dt_Rfactor, F=F, T=T,
        #               user_action=compute_error)
        # Testing: ------------------------ should give r=1
        diffusion_FE(I=I, a=a, f=f, L=L, dt=dt, F=F, T=T,
                           user_action=compute_error)
        # -----------------------------------
        
        # error is computed in the final call to compute_error
        E.append(error)
        #h.append(dt)
        h.append(dt/float(dt_Rfactor)) # E found in reaction part
        dt /= 2  # halve the time step for next simulation
    print 'E:', E
    print 'h:', h
    # Convergence rates for two consecutive experiments
    r = [np.log(E[i]/E[i-1])/np.log(h[i]/h[i-1])
         for i in range(1,num_meshes)]
    return r

def test_convrate_sinexp():
    #n = m = 2
    #L = 1.0
    #u_exact = lambda x, t: np.cos(m*np.pi/L*t)*np.sin(m*np.pi/L*x)

    #sl: arrange so this fu can handle all cases. In particular, f must 
    # either be defined here (no split) or not.
    # Just finish ordinary_splitting first...

    a = 3.5
    b = 1
    L = 1.5
    Nx = 4
    #Nx = 40
    k = np.pi/L
    dx = L/Nx
    dt0 = dx**2/(4*a)   # max dt according to stability req.
    dt_Rfactor=1   # local time step in reaction-part is dt/dt_Rfactor

    def u_exact(x, t):
        '''exact sol. to: du/dt = a*d^2u/dx^2 - b*u'''
        return np.exp(-(a*k**2 + b)*t) * np.sin(k*x)

    # f is defined within reaction part itself when splitting is done.
    # If no splitting, f must be defined here and passed on
    #f = 0      # when splitting
    f = lambda u, t: -b*u   # when FE all the way

    r = convergence_rates(
        u_exact=u_exact,
        I=lambda x: u_exact(x, 0),
        f=f,            
        a=a,
        b=b,
        L=L,
        dx=dx,
        dt_Rfactor=dt_Rfactor,
        dt0=dt0,        
        num_meshes=6,
        T=0.2)
    print 'Computed rates:', [round(r_,2) for r_ in r]
    assert abs(r[-1] - 2) < 0.002


if __name__ == '__main__':
    #if len(sys.argv) < 2:
    #    print """Usage %s function arg1 arg2 arg3 ...""" % sys.argv[0]
    #    sys.exit(0)
    #cmd = '%s(%s)' % (sys.argv[1], ', '.join(sys.argv[2:]))
    #print cmd
    #eval(cmd)
    
    #test_solvers()
    test_convrate_sinexp()
    