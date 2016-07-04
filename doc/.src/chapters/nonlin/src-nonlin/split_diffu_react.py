import sys, time
import scitools.std as plt
import scipy.sparse
import scipy.sparse.linalg
import numpy as np

def diffusion_FE(I, a, f, L, dt, F, t, T, step_no, user_action=None):
    """Diffusion solver, Forward Euler method.
    Note that t always covers the whole global time interval, whether
    splitting is the case or not. T, on the other hand, is 
    the end of the global time interval if there is no split,
    but if splitting, we use T=dt. When splitting, step_no keeps
    track of the time step number (required for lookup in t).
    """

    Nt = int(round(T/float(dt)))
    #t = np.linspace(0, Nt*dt, Nt+1)   # Mesh points in time
    dx = np.sqrt(a*dt/F)
    Nx = int(round(L/dx))
    x = np.linspace(0, L, Nx+1)       # Mesh points in space
    # Make sure dx and dt are compatible with x and t
    dx = x[1] - x[0]
    dt = t[1] - t[0]

    u   = np.zeros(Nx+1)   # solution array
    u_1 = np.zeros(Nx+1)   # solution at t-dt

    # Allow f to be None or 0
    if f is None or f == 0:
        f = lambda x, t: np.zeros((x.size)) \
            if isinstance(x, np.ndarray) else 0

    # Set initial condition  
    if isinstance(I, np.ndarray):   # I is an array
        u_1 = np.copy(I)
    else:                           # I is a function
        for i in range(0, Nx+1):
            u_1[i] = I(x[i])

    if user_action is not None:
        user_action(u_1, x, t, step_no+0)

    for n in range(0, Nt):
        # Update all inner points
        u[1:Nx] = u_1[1:Nx] +  \
                  F*(u_1[0:Nx-1] - 2*u_1[1:Nx] + u_1[2:Nx+1]) +\
                  dt*f(u_1[1:Nx], t[step_no+n])
                  #dt*f(x[1:Nx], t[n])   # if f(x,t), not f(u,t)
        
        # Insert boundary conditions
        u[0] = 0;  u[Nx] = 0

        #sl: ...testing -------------------------
        #print 'time:', t[step_no+n]
        #print 'diff part from diffusion_FE:'
        #print u_1[1:Nx] + F*(u_1[0:Nx-1] - 2*u_1[1:Nx] + u_1[2:Nx+1])
        #print 'react part from diffusion_FE:'
        #print dt*f(u_1[1:Nx], t[step_no+n])
        #print ' '
        #if step_no == 1: sys.exit(0)            
        # ----------------------------------

        if user_action is not None:
            user_action(u, x, t, step_no+(n+1))

        # Switch variables before next step
        u_1, u = u, u_1

    return u_1



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
    u_1 = np.zeros(Nx+1)   # solution at t[n]

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
        b[1:-1] = u_1[1:-1] + \
                  Fr*(u_1[:-2] - 2*u_1[1:-1] + u_1[2:]) + \
                  dt*theta*f(x[1:-1], t[n+1]) + \
                  dt*(1-theta)*f(x[1:-1], t[n])
        b[0] = u_L; b[-1] = u_R  # boundary conditions
        u[:] = scipy.sparse.linalg.spsolve(A, b)

        if user_action is not None:
            user_action(u, x, t, n+1)

        # Update u_1 before next step
        u_1, u = u, u_1

    t1 = time.clock()
    return t1-t0


def reaction_FE(I, f, L, Nx, dt, dt_Rfactor, t, step_no, 
                user_action=None):
    """Reaction solver, Forward Euler method.
    Note that t covers the whole global time interval. 
    dt is the step of the diffustion part, i.e. there
    is a local time interval [0, dt] the reaction_FE
    deals with each time it is called. step_no keeps
    track of the (global) time step number (required 
    for lookup in t).    
    """

    #bypass = True
    #if not bypass:  # original code from sl

    u = np.copy(I)      
    dt_local = dt/float(dt_Rfactor)    
    Nt_local = int(round(dt/float(dt_local)))    
    x = np.linspace(0, L, Nx+1)  
          
    for n in range(Nt_local):
        time = t[step_no] + n*dt_local
        u[1:Nx] = u[1:Nx] + dt_local*f(u[1:Nx], time)  
        
    # BC already inserted in diffusion step, i.e. no action here            
                        
    return u
        
    #else:
    #    return I        
        

def ordinary_splitting(I, a, b, f, L, dt, 
                       dt_Rfactor, F, t, T, 
                       user_action=None):
    '''1st order scheme, i.e. Forward Euler is enough for both
    the diffusion and the reaction part. The time step dt is 
    given for the diffusion step, while the time step for the 
    reaction part is found as dt/dt_Rfactor, where dt_Rfactor >= 1.
    '''
    Nt = int(round(T/float(dt)))
    #t = np.linspace(0, Nt*dt, Nt+1)   # Mesh points, global time
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
        # Note: could avoid the call to diffusion_FE here...
        
        # Diffusion step (one time step dt)
        u_s = diffusion_FE(I=u, a=a, f=0, L=L, dt=dt, F=F, 
                           t=t, T=dt, step_no=n,
                           user_action=None)
        # Reaction step (potentially many smaller steps within dt)
        u = reaction_FE(I=u_s, f=f, L=L, Nx=Nx, 
                        dt=dt, dt_Rfactor=dt_Rfactor, 
                        t=t, step_no=n, 
                        user_action=None)

        if user_action is not None:
            user_action(u, x, t, n+1)


def Strang_splitting_1stOrder(I, a, b, f, L, dt, 
                              dt_Rfactor, F, t, T, 
                              user_action=None):
    '''Strang splitting while still using FE for the diffusion
    step and for the reaction step. Gives 1st order scheme.
    Introduce an extra time mesh t2 for the diffusion part, 
    since it makes twice as many steps (step is dt/2).
    '''
    Nt = int(round(T/float(dt)))
    # diff step makes twice the number of steps ...
    t2 = np.linspace(0, Nt*dt, (Nt+1)+Nt)   # Mesh points in diff    
    dx = np.sqrt(a*dt/F)
    Nx = int(round(L/dx))
    x = np.linspace(0, L, Nx+1)       # Mesh points in space
    u = np.zeros(Nx+1)
        
    # Set initial condition u(x,0) = I(x)
    for i in range(0, Nx+1):
        u[i] = I(x[i])

    for n in range(0, Nt):    
        
        # Diffusion step (1/2 dt: from t_n to t_n+1/2)
        u_s = diffusion_FE(I=u, a=a, f=0, L=L, dt=dt/2.0, F=F/2.0, 
                           t=t2, T=dt/2.0, step_no=2*n,
                           user_action=None)
                                                                             
        # Reaction step (1 dt: from t_n to t_n+1)
        # (potentially many smaller steps within dt)
        u = reaction_FE(I=u_s, f=f, L=L, Nx=Nx, 
                        dt=dt, dt_Rfactor=dt_Rfactor, 
                        t=t, step_no=n, 
                        user_action=None)

        # Diffusion step (1/2 dt: from t_n+1/2 to t_n)
        u_s = diffusion_FE(I=u, a=a, f=0, L=L, dt=dt/2.0, F=F/2.0, 
                           t=t2, T=dt/2.0, step_no=2*n+1,
                           user_action=None)

        # set init. cond. for next step
        u[:] = u_s[:]

        if user_action is not None:
            user_action(u, x, t, n+1)




def convergence_rates(scheme='diffusion'):
    '''Computes empirical conv. rates for the different
    splitting schemes'''
    
    F = 0.5
    T = 1.2
    a = 3.5
    b = 1
    L = 1.5
    k = np.pi/L
    
    def exact(x, t):
        '''exact sol. to: du/dt = a*d^2u/dx^2 - b*u'''
        return np.exp(-(a*k**2 + b)*t) * np.sin(k*x)
        
    def f(u, t):
        return -b*u  

    def I(x):
        return exact(x, 0)
        
    global error    # error computed in the user action function
    error = 0  

    # Convergence study
    def action(u, x, t, n):
        global error            
        if n == 1:      # New simulation, - reset error
            error = 0
        else:
            error = max(error, np.abs(u - exact(x, t[n])).max())        

    E = []
    h = []
    Nx_values = [10, 20, 40, 80, 160]   
    for Nx in Nx_values:        
        dx = L/Nx;  dt = F/a*dx**2
        Nt = int(round(T/float(dt)))
        t = np.linspace(0, Nt*dt, Nt+1)   # Mesh points, global time
        
        if scheme == 'diffusion':
            print 'Running FE on whole eqn...'
            diffusion_FE(I, a, f, L, dt, F, t, T, 
                         step_no=0, user_action=action)
        elif scheme == 'ordinary_splitting':
            print 'Running ordinary splitting...'
            ordinary_splitting(I=I, a=a, b=b, f=f, L=L, dt=dt, 
                               dt_Rfactor=1, F=F, t=t, T=T,
                               user_action=action)        
        elif scheme == 'Strang_splitting_1stOrder':
            print 'Running Strang splitting with 1st order schemes...'
            Strang_splitting_1stOrder(I=I, a=a, b=b, f=f, L=L, dt=dt, 
                                      dt_Rfactor=1, F=F, t=t, T=T,
                                      user_action=action)       
        elif scheme == 'Strang_splitting_2ndOrder':
            print 'Running Strang splitting with 2nd order schemes...'
            print 'Not yet implemented.'        
            sys.exit(0)
        else:
            print 'Unknown scheme requested!'
            sys.exit(0)
                           
                           
        h.append(dt)
        E.append(error)

    print 'E:', E
    print 'h:', h
        
    # Convergence rates 
    r = [np.log(E[i]/E[i-1])/np.log(h[i]/h[i-1])
         for i in range(1,len(Nx_values))]
    print 'Computed rates:', r


if __name__ == '__main__':
    
    schemes = ['diffusion', 
               'ordinary_splitting', 
               'Strang_splitting_1stOrder',
               'Strang_splitting_2ndOrder']
    for scheme in schemes:
        convergence_rates(scheme=scheme)
    
