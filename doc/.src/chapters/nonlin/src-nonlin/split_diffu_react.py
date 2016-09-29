import sys, time
import scitools.std as plt
import scipy.sparse
import scipy.sparse.linalg
import numpy as np

def diffusion_theta(I, a, f, L, dt, F, t, T, step_no, theta=0.5, 
                    u_L=0, u_R=0, user_action=None):
    """
    Full solver for the model problem using the theta-rule
    difference approximation in time (no restriction on F,
    i.e., the time step when theta >= 0.5).     Vectorized 
    implementation and sparse (tridiagonal) coefficient matrix.
    Note that t always covers the whole global time interval, whether
    splitting is the case or not. T, on the other hand, is 
    the end of the global time interval if there is no split,
    but if splitting, we use T=dt. When splitting, step_no 
    keeps track of the time step number (for lookup in t).    
    """

    Nt = int(round(T/float(dt)))
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

    # Time loop
    for n in range(0, Nt):
        b[1:-1] = u_1[1:-1] + \
                  Fr*(u_1[:-2] - 2*u_1[1:-1] + u_1[2:]) + \
                  dt*theta*f(u_1[1:-1], t[step_no+n+1]) + \
                  dt*(1-theta)*f(u_1[1:-1], t[step_no+n])
        b[0] = u_L; b[-1] = u_R  # boundary conditions
        u[:] = scipy.sparse.linalg.spsolve(A, b)

        if user_action is not None:
            user_action(u, x, t, step_no+(n+1))

        # Update u_1 before next step
        u_1, u = u, u_1

    # u is now contained in u_1 (swapping)
    return u_1


def reaction_FE(I, f, L, Nx, dt, dt_Rfactor, t, step_no, 
                user_action=None):
    """Reaction solver, Forward Euler method.
    Note the at t covers the whole global time interval. 
    dt is either one complete,or one half, of the step in the 
    diffusion part, i.e. there is a local time interval 
    [0, dt] or [0, dt/2] that the reaction_FE
    deals with each time it is called. step_no keeps
    track of the (global) time step number (required 
    for lookup in t).    
    """

    u = np.copy(I)      
    dt_local = dt/float(dt_Rfactor)    
    Nt_local = int(round(dt/float(dt_local)))    
    x = np.linspace(0, L, Nx+1)  
          
    for n in range(Nt_local):
        time = t[step_no] + n*dt_local
        u[1:Nx] = u[1:Nx] + dt_local*f(u[1:Nx], time)  
        
    # BC already inserted in diffusion step, i.e. no action here            
    return u

def ordinary_splitting(I, a, b, f, L, dt, 
                       dt_Rfactor, F, t, T, 
                       user_action=None):
    '''1st order scheme, i.e. Forward Euler is enough for both
    the diffusion and the reaction part. The time step dt is 
    given for the diffusion step, while the time step for the 
    reaction part is found as dt/dt_Rfactor, where dt_Rfactor >= 1.
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
    # first for reaction, then for diffusion
    for n in range(0, Nt):    
        # Reaction step (potentially many smaller steps within dt)
        u_s = reaction_FE(I=u, f=f, L=L, Nx=Nx, 
                        dt=dt, dt_Rfactor=dt_Rfactor, 
                        t=t, step_no=n, 
                        user_action=None)        
        
        u = diffusion_theta(I=u_s, a=a, f=0, L=L, dt=dt, F=F,
                              t=t, T=dt, step_no=n, theta=0,
                              u_L=0, u_R=0, user_action=None)                         
                           
        if user_action is not None:
            user_action(u, x, t, n+1)

    return

def Strang_splitting_1stOrder(I, a, b, f, L, dt, dt_Rfactor, 
                              F, t, T, user_action=None):
    '''Strang splitting while still using FE for the reaction
    step and for the diffusion step. Gives 1st order scheme.
    The time step dt is given for the diffusion step, while 
    the time step for the reaction part is found as 
    0.5*dt/dt_Rfactor, where dt_Rfactor >= 1. Introduce an 
    extra time mesh t2 for the reaction part, since it steps dt/2.
    '''
    Nt = int(round(T/float(dt)))
    t2 = np.linspace(0, Nt*dt, (Nt+1)+Nt)   # Mesh points in diff    
    dx = np.sqrt(a*dt/F)
    Nx = int(round(L/dx))
    x = np.linspace(0, L, Nx+1)     
    u = np.zeros(Nx+1)
        
    # Set initial condition u(x,0) = I(x)
    for i in range(0, Nx+1):
        u[i] = I(x[i])

    for n in range(0, Nt):    
        # Reaction step (1/2 dt: from t_n to t_n+1/2)
        # (potentially many smaller steps within dt/2)
        u_s = reaction_FE(I=u, f=f, L=L, Nx=Nx, 
                        dt=dt/2.0, dt_Rfactor=dt_Rfactor, 
                        t=t2, step_no=2*n, 
                        user_action=None)
        # Diffusion step (1 dt: from t_n to t_n+1)
        u_sss = diffusion_theta(I=u_s, a=a, f=0, L=L, dt=dt, F=F,
                              t=t, T=dt, step_no=n, theta=0,
                              u_L=0, u_R=0, user_action=None)                         
        # Reaction step (1/2 dt: from t_n+1/2 to t_n+1)
        # (potentially many smaller steps within dt/2)
        u = reaction_FE(I=u_sss, f=f, L=L, Nx=Nx, 
                        dt=dt/2.0, dt_Rfactor=dt_Rfactor, 
                        t=t2, step_no=2*n+1, 
                        user_action=None)
        
        if user_action is not None:
            user_action(u, x, t, n+1)

    return

def Strang_splitting_2ndOrder(I, a, b, f, L, dt, dt_Rfactor, 
                              F, t, T, user_action=None):
    '''Strang splitting using Crank-Nicolson for the diffusion
    step (theta-rule) and Adams-Bashforth 2 for the reaction step. 
    Gives 2nd order scheme. Introduce an extra time mesh t2 for 
    the reaction part, since it steps dt/2.
    '''
    import odespy
    Nt = int(round(T/float(dt)))
    t2 = np.linspace(0, Nt*dt, (Nt+1)+Nt)   # Mesh points in diff    
    dx = np.sqrt(a*dt/F)
    Nx = int(round(L/dx))
    x = np.linspace(0, L, Nx+1)     
    u = np.zeros(Nx+1)
        
    # Set initial condition u(x,0) = I(x)
    for i in range(0, Nx+1):
        u[i] = I(x[i])

    reaction_solver = odespy.AdamsBashforth2(f)    

    for n in range(0, Nt):    
        # Reaction step (1/2 dt: from t_n to t_n+1/2)
        # (potentially many smaller steps within dt/2)
        reaction_solver.set_initial_condition(u)
        t_points = np.linspace(0, dt/2.0, dt_Rfactor+1)
        u_AB2, t_ = reaction_solver.solve(t_points) # t_ not needed
        u_s = u_AB2[-1,:]  # pick sol at last point in time

        # Diffusion step (1 dt: from t_n to t_n+1)
        u_sss = diffusion_theta(I=u_s, a=a, f=0, L=L, dt=dt, F=F,
                              t=t, T=dt, step_no=n, theta=0.5,
                              u_L=0, u_R=0, user_action=None)                         
        # Reaction step (1/2 dt: from t_n+1/2 to t_n+1)
        # (potentially many smaller steps within dt/2)
        reaction_solver.set_initial_condition(u_sss)
        t_points = np.linspace(0, dt/2.0, dt_Rfactor+1)
        u_AB2, t_ = reaction_solver.solve(t_points) # t_ not needed
        u = u_AB2[-1,:]  # pick sol at last point in time
               
        if user_action is not None:
            user_action(u, x, t, n+1)

    return 

def convergence_rates(scheme='diffusion'):
    
    F = 0.5     # Upper limit for FE (stability). For CN, this 
                # limit does not apply, but for simplicity, we
                # choose F = 0.5 as the initial F value.
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
    Nx_values = [10, 20, 40, 80]   # i.e., dx halved each time
    for Nx in Nx_values:        
        dx = L/Nx           
        if scheme == 'Strang_splitting_2ndOrder':
            print 'Strang splitting with 2nd order schemes...'
            # In this case, E = C*h**r (with r = 2) and since 
            # h = dx = K*dt, the ratio dt/dx must be constant.
            # To fulfill this demand, we must let F change
            # when dx changes. From F = a*dt/dx**2, it follows
            # that halving dx AND doubling F assures dt/dx const.
            # Initially, we simply choose F = 0.5.

            dt = F/a*dx**2
            #print 'dt/dx:', dt/dx            
            Nt = int(round(T/float(dt)))
            t = np.linspace(0, Nt*dt, Nt+1)   # global time            
            Strang_splitting_2ndOrder(I=I, a=a, b=b, f=f, L=L, dt=dt, 
                                      dt_Rfactor=1, F=F, t=t, T=T,
                                      user_action=action)    
            h.append(dx)
            # prepare for next iteration (make F match dx/2)
            F = F*2       # assures dt/dx const. when dx = dx/2                                      
        else:   
            # In these cases, E = C*h**r (with r = 1) and since 
            # h = dx**2 = K*dt, the ratio dt/dx**2 must be constant.
            # This is fulfilled by choosing F = 0.5 (for FE stability)
            # and make sure that F, dx and dt comply to F = a*dt/dx**2.            
            dt = F/a*dx**2
            Nt = int(round(T/float(dt)))
            t = np.linspace(0, Nt*dt, Nt+1)   # global time
            if scheme == 'diffusion':
                print 'FE on whole eqn...'
                diffusion_theta(I, a, f, L, dt, F, t, T,
                                step_no=0, theta=0,
                                u_L=0, u_R=0, user_action=action)                         
                h.append(dx**2)
            elif scheme == 'ordinary_splitting':
                print 'Ordinary splitting...'
                ordinary_splitting(I=I, a=a, b=b, f=f, L=L, dt=dt, 
                                   dt_Rfactor=1, F=F, t=t, T=T,
                                   user_action=action)        
                h.append(dx**2)
            elif scheme == 'Strang_splitting_1stOrder':
                print 'Strang splitting with 1st order schemes...'
                Strang_splitting_1stOrder(I=I, a=a, b=b, f=f, L=L, dt=dt, 
                                          dt_Rfactor=1, F=F, t=t, T=T,
                                          user_action=action)       
                h.append(dx**2)
            else:
                print 'Unknown scheme requested!'
                sys.exit(0)
        
            #print 'dt/dx**2:', dt/dx**2            
        
        E.append(error)
        Nx *= 2         # Nx doubled gives dx/2

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
    
