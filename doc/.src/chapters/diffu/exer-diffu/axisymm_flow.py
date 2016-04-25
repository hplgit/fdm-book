"""
Solve the diffusion equation for axi-symmetric case:

    u_t = 1/r * (r*a(r)*u_r)_r + f(r,t)

on (0,R) with boundary conditions u(0,t)_r = 0 and u(R,t) = 0,
for t in (0,T]. Initial condition: u(r,0) = I(r). 
Pressure gradient f.

The following naming convention of variables are used.

===== ==========================================================
Name  Description
===== ==========================================================
Nx    The total number of mesh cells; mesh points are numbered
      from 0 to Nx.
T     The stop time for the simulation.
I     Initial condition (Python function of x).
a     Variable coefficient (constant).
R     Length of the domain ([0,R]).
r     Mesh points in space.
t     Mesh points in time.
n     Index counter in time.
u     Unknown at current/new time level.
u_1   u at the previous time level.
dr    Constant mesh spacing in r.
dt    Constant mesh spacing in t.
===== ==========================================================

``user_action`` is a function of ``(u, r, t, n)``, ``u[i]`` is the
solution at spatial mesh point ``r[i]`` at time ``t[n]``, where the
calling code can add visualization, error computations, data analysis,
store solutions, etc.
"""

import scipy.sparse
import scipy.sparse.linalg
from numpy import linspace, zeros, random, array, ones, sum, log, sqrt
import time, sys
import sympy as sym    


def solver_theta(I, a, R, Nr, D, T, theta=0.5, u_L=None, u_R=0,
                 user_action=None, f=0):
    """
    The array a has length Nr+1 and holds the values of
    a(x) at the mesh points.

    Method: (implicit) theta-rule in time.

    Nr is the total number of mesh cells; mesh points are numbered
    from 0 to Nr.
    D = dt/dr**2 and implicitly specifies the time step.
    T is the stop time for the simulation.
    I is a function of r.
    u_L = None implies du/dr = 0, i.e. a symmetry condition 
    f(r,t) is pressure gradient with radius.

    user_action is a function of (u, x, t, n) where the calling code
    can add visualization, error computations, data analysis,
    store solutions, etc.
    
    r*alpha is needed midway between spatial mesh points, - use
    arithmetic mean of successive mesh values (i.e. of r_i*alpha_i)
    """
    import time
    t0 = time.clock()

    r = linspace(0, R, Nr+1)   # mesh points in space
    dr = r[1] - r[0]
    dt = D*dr**2   
    Nt = int(round(T/float(dt)))
    t = linspace(0, T, Nt+1)   # mesh points in time

    if isinstance(u_L, (float,int)):
        u_L_ = float(u_L)  # must take copy of u_L number
        u_L = lambda t: u_L_
    if isinstance(u_R, (float,int)):
        u_R_ = float(u_R)  # must take copy of u_R number
        u_R = lambda t: u_R_
    if isinstance(f, (float,int)):
        f_ = float(f)  # must take copy of f number
        f = lambda r, t: f_

    ra = r*a    # help array in scheme

    inv_r = zeros(len(r)-2)    # needed for inner mesh points
    inv_r = 1.0/r[1:-1]

    u   = zeros(Nr+1)   # solution array at t[n+1]
    u_1 = zeros(Nr+1)   # solution at t[n]

    Dl = 0.5*D*theta
    Dr = 0.5*D*(1-theta)

    # Representation of sparse matrix and right-hand side
    diagonal = zeros(Nr+1)
    lower    = zeros(Nr)
    upper    = zeros(Nr)
    b        = zeros(Nr+1)

    # Precompute sparse matrix (scipy format)
    diagonal[1:-1] = 1 + Dl*(ra[2:] + 2*ra[1:-1] + ra[:-2])*inv_r
    lower[:-1] = -Dl*(ra[1:-1] + ra[:-2])*inv_r
    upper[1:]  = -Dl*(ra[2:] + ra[1:-1])*inv_r
    # Insert boundary conditions
    if u_L == None:     # symmetry axis, du/dr = 0
        diagonal[0] = 1 + 8*a[0]*Dl
        upper[0] = -8*a[0]*Dl
    else:
        diagonal[0] = 1
        upper[0] = 0
    diagonal[Nr] = 1
    lower[-1] = 0

    A = scipy.sparse.diags(
        diagonals=[diagonal, lower, upper],
        offsets=[0, -1, 1],
        shape=(Nr+1, Nr+1),
        format='csr')
    #print A.todense()

    # Set initial condition
    for i in range(0,Nr+1):
        u_1[i] = I(r[i])

    if user_action is not None:
        user_action(u_1, r, t, 0)

    # Time loop
    for n in range(0, Nt):
        b[1:-1] = u_1[1:-1] + Dr*(
            (ra[2:] + ra[1:-1])*(u_1[2:] - u_1[1:-1]) -
            (ra[1:-1] + ra[0:-2])*(u_1[1:-1] - u_1[:-2]))*inv_r + \
            dt*theta*f(r[1:-1], t[n+1]) + \
            dt*(1-theta)*f(r[1:-1], t[n])
            
        # Boundary conditions
        if u_L == None:     # symmetry axis, du/dr = 0
            b[0] = u_1[0] + 8*a[0]*Dr*(u_1[1] - u_1[0]) + \
                   dt*theta*f(0, (n+1)*dt) + \
                   dt*(1 - theta)*f(0, n*dt)
        else:               
            b[0]  = u_L(t[n+1])        
        b[-1] = u_R(t[n+1])
        #print b        
        
        # Solve
        u[:] = scipy.sparse.linalg.spsolve(A, b)
        
        if user_action is not None:
            user_action(u, r, t, n+1)

        # Switch variables before next step
        u_1, u = u, u_1

    t1 = time.clock()
    # return u_1, since u and u_1 are switched
    return u_1, t, t1-t0

def compute_rates(h_values, E_values):
    m = len(h_values)
    q = [log(E_values[i+1]/E_values[i])/
         log(h_values[i+1]/h_values[i])
         for i in range(0, m-1, 1)]
    q = [round(q_, 2) for q_ in q]
    return q

def make_a(alpha, r):
    """
    alpha is a func, generally of r, - but may be constant.
    Note: when solution is to be axi-symmetric, alpha
    must be so too.
    """
    a = alpha(r)*ones(len(r))
    return a

def tests_with_alpha_and_u_exact():
    '''
    Test solver performance when alpha is either const or 
    a fu of r, combined with a manufactured sol u_exact 
    that is either a fu of r only, or a fu of both r and t.
    Note: alpha and u_e are defined as symb expr here, since 
    test_solver_symmetric needs to automatically generate 
    the source term f. After that, test_solver_symmetric
    redefines alpha, u_e and f as num functions.
    '''
    R, r, t = sym.symbols('R r t')

    # alpha const ...
    
    # ue = const
    print 'Testing with alpha = 1.5 and u_e = R**2 - r**2...'
    test_solver_symmetric(alpha=1.5, u_exact=R**2 - r**2)
    
    # ue = ue(t)
    print 'Testing with alpha = 1.5 and u_e = 5*t*(R**2 - r**2)...'
    test_solver_symmetric(alpha=1.5, u_exact=5*t*(R**2 - r**2))
    
    # alpha function of r ...
    
    # ue = const 
    print 'Testing with alpha = 1 + r**2 and u_e = R**2 - r**2...'
    test_solver_symmetric(alpha=1+r**2, u_exact=R**2 - r**2)
    
    # ue = ue(t)
    print 'Testing with alpha = 1+r**2 and u_e = 5*t*(R**2 - r**2)...'
    test_solver_symmetric(alpha=1+r**2, u_exact=5*t*(R**2 - r**2))



def test_solver_symmetric(alpha, u_exact):
    '''
    Test solver performance for manufactured solution
    given in the function u_exact. Parameter alpha is 
    either a const or a function of r. In the latter 
    case, an "exact" sol can not be achieved, so then
    testing switches to conv. rates.
    R is tube radius and T is duration of simulation.
    alpha constant:
        Compares the manufactured solution with the 
        solution from the solver at each time step. 
    alpha function of r:
        convergence rates are tested (using the sol
        at the final point in time only).
    '''   
    
    def compare(u, r, t, n):      # user_action function
        """Compare exact and computed solution."""
        u_e = u_exact(r, t[n])
        diff = abs(u_e - u).max()
        #print diff
        tol = 1E-12
        assert diff < tol, 'max diff: %g' % diff

    def pde_source_term(a, u):
        '''Return the terms in the PDE that the source term
        must balance, here du/dt - (1/r) * d/dr(r*a*du/dr).
        a, i.e. alpha, is either const or a fu of r.
        u is a symbolic Python function of r and t.'''
        
        return sym.diff(u, t) - \
               (1.0/r)*sym.diff(r*a*sym.diff(u, r), r)
               
    R, r, t = sym.symbols('R r t')

    # fit source term
    f = sym.simplify(pde_source_term(alpha, u_exact))  

    R = 1.0     # radius of tube
    T = 2.0     # duration of simulation 
   
    if sym.diff(alpha, r) == 0:  
        alpha_is_const = True
    else:
        alpha_is_const = False        

    # make alpha, f and u_exact numerical functions
    alpha = sym.lambdify([r], alpha, modules='numpy')             
    f = sym.lambdify([r, t], f.subs('R', R), modules='numpy')             
    u_exact = sym.lambdify(
        [r, t], u_exact.subs('R', R), modules='numpy')             

    I = lambda r: u_exact(r, 0)

    # some help variables
    FE = 0      # Forward Euler method
    BE = 1      # Backward Euler method
    CN = 0.5    # Crank-Nicolson method

    # test all three schemes 
    for theta in (FE, BE, CN):
        print 'theta: ', theta
        E_values = []
        dt_values = []
        for Nr in (2, 4, 8, 16, 32, 64):
            print 'Nr:', Nr
            r = linspace(0, R, Nr+1)   # mesh points in space
            dr = r[1] - r[0]
            a_values = make_a(alpha, r)   
            if theta == CN:
                dt = dr
            else:   # either FE or BE
                # use most conservative dt as decided by FE
                K = 1.0/(4*a_values.max())              
                dt = K*dr**2                 
            D = dt/dr**2

            if alpha_is_const:  
                u, t, cpu = solver_theta(
                        I, a_values, R, Nr, D, T, 
                        theta, u_L=None, u_R=0,
                        user_action=compare, f=f)   
            else:   # alpha depends on r
                u, t, cpu = solver_theta(
                        I, a_values, R, Nr, D, T, 
                        theta, u_L=None, u_R=0,
                        user_action=None, f=f)   
                        
                # compute L2 error at t = T
                u_e = u_exact(r, t[-1])
                e = u_e - u
                E = sqrt(dr*sum(e**2))
                E_values.append(E)
                dt_values.append(dt)
            
        if alpha_is_const is False:  
            q = compute_rates(dt_values, E_values)        
            print 'theta=%g, q: %s' % (theta, q)
            expected_rate = 2 if theta == CN else 1
            tol = 0.1
            diff = abs(expected_rate - q[-1])
            print 'diff:', diff
            assert diff < tol
    
    
if __name__ == '__main__':
    tests_with_alpha_and_u_exact()        
    print 'This is just a start. More remaining for this Exerc.'