import sympy as sym
import numpy as np

# The code in vib_undamped_verify_mms.py is here generalized
# to treat the model m*u'' + f(u') + c*u = F(t), where the 
# damping term f(u') = 0, b*u' or b*V*abs(V). 

def ode_source_term(u, damping):
    """Return the terms in the ODE that the source term
    must balance, here m*u'' + f(u') + c*u.
    u is a symbolic Python function of t."""
    if damping == 'zero':
        return m*sym.diff(u(t), t, t) + c*u(t)  
    elif damping == 'linear':
        return m*sym.diff(u(t), t, t) + \
               b*sym.diff(u(t), t) + c*u(t)    
    else:  # damping is nonlinear
        return m*sym.diff(u(t), t, t) + \
               b*sym.diff(u(t), t)*abs(sym.diff(u(t), t)) + c*u(t) 

def residual_discrete_eq(u, damping):
    """Return the residual of the discrete eq. with u inserted."""
    if damping == 'zero':
        R = m*DtDt(u, dt) + c*u(t) - F
    elif damping == 'linear':
        R = m*DtDt(u, dt) + b*D2t(u, dt) + c*u(t) - F
    else:   # damping is nonlinear
        R = m*DtDt(u, dt) + b*Dt_p_half(u, dt)*\
            abs(Dt_m_half(u, dt)) + c*u(t) - F
    return sym.simplify(R)

def residual_discrete_eq_step1(u, damping):
    """Return the residual of the discrete eq. at the first
    step with u inserted."""
    half = sym.Rational(1,2)
    if damping == 'zero':
        R = u(t+dt) - u(t) - dt*V - \
            half*dt**2*(F.subs(t, 0)/m) + half*dt**2*(c/m)*I
    elif damping == 'linear':
        R = u(t+dt) - (I + dt*V + \
            half*(dt**2/m)*(-b*V - c*I + F.subs(t, 0)))
    else:   # damping is nonlinear
        R = u(t+dt) - (I + dt*V + \
            half*(dt**2/m)*(-b*V*abs(V) - c*I + F.subs(t, 0)))
    R = R.subs(t, 0)  # t=0 in the rhs of the first step eq.
    return sym.simplify(R)

def DtDt(u, dt):
    """Return 2nd-order finite difference for u_tt.
    u is a symbolic Python function of t.
    """
    return (u(t+dt) - 2*u(t) + u(t-dt))/dt**2

def D2t(u, dt):
    """Return 2nd-order finite difference for u_t.
    u is a symbolic Python function of t.
    """
    return (u(t+dt) - u(t-dt))/(2.0*dt)

def Dt_p_half(u, dt):
    """Return 2nd-order finite difference for u_t, sampled at n+1/2,
    i.e, pluss one half... u is a symbolic Python function of t.
    """
    return (u(t+dt) - u(t))/dt

def Dt_m_half(u, dt):
    """Return 2nd-order finite difference for u_t, sampled at n-1/2,
    i.e, minus one half.... u is a symbolic Python function of t.
    """
    return (u(t) - u(t-dt))/dt

def main(u, damping):
    """
    Given some chosen solution u (as a function of t, implemented
    as a Python function), use the method of manufactured solutions
    to compute the source term f, and check if u also solves
    the discrete equations.
    """
    print '=== Testing exact solution: %s ===' % u(t)
    print "Initial conditions u(0)=%s, u'(0)=%s:" % \
          (u(t).subs(t, 0), sym.diff(u(t), t).subs(t, 0))

    # Method of manufactured solution requires fitting F
    global F  # source term in the ODE
    F = sym.simplify(ode_source_term(u, damping))

    # Residual in discrete equations (should be 0)
    print 'residual step1:', residual_discrete_eq_step1(u, damping)
    print 'residual:', residual_discrete_eq(u, damping)


def linear(damping):
    def u_e(t):
        """Return chosen linear exact solution."""
        # General linear function u_e = c*t + d
        # Initial conditions u(0)=I, u'(0)=V require c=V, d=I
        return V*t + I

    main(u_e, damping)

def quadratic(damping):
    # Extend with quadratic functions
    q = sym.Symbol('q')  # arbitrary constant in quadratic term

    def u_e(t):
        return q*t**2 + V*t + I

    main(u_e, damping)

def cubic(damping):
    r, q = sym.symbols('r q')

    main(lambda t: r*t**3 + q*t**2 + V*t + I, damping)


def solver(I, V, F, c, m, dt, T, damping):
    """
    Solve m*u'' + f(u') + c*u = F for t in (0,T], u(0)=I and u'(0)=V,
    by a central finite difference method with time step dt.
    F(t) is a callable Python function.
    """
    dt = float(dt)
    Nt = int(round(T/dt))
    u = np.zeros(Nt+1)
    t = np.linspace(0, Nt*dt, Nt+1)

    if damping == 'zero':
        u[0] = I
        u[1] = u[0] - 0.5*dt**2*(c/m)*u[0] + 0.5*dt**2*F(t[0]) + dt*V
        for n in range(1, Nt):
            u[n+1] = 2*u[n] - u[n-1] - dt**2*(c/m)*u[n] + dt**2*F(t[n])
    elif damping == 'linear':
        u[0] = I
        u[1] = u[0] + dt*V + 0.5*(dt**2/m)*(-b*V - c*u[0] + F(t[0]))
        for n in range(1, Nt):
            u[n+1] = (2*m*u[n] + (b*dt/2.-m)*u[n-1] + \
                     dt**2*(F(t[n])-c*u[n]))/(m+b*dt/2.)
    else:    # damping is nonlinear
        u[0] = I
        u[1] = u[0] + dt*V + 0.5*(dt**2/m)*(-b*V*abs(V) - c*u[0] + F(t[0]))
        for n in range(1, Nt):
            u[n+1] = 1./(m+b*abs(u[n]-u[n-1])) * \
                     (2*m*u[n] - m*u[n-1] + b*u[n]*\
                     abs(u[n]-u[n-1])+dt**2*(F(t[n])-c*u[n]))
    return u, t

import nose.tools as nt

def test_quadratic_exact_solution(damping):
    # Transform global symbolic variables to functions and numbers
    # for numerical computations

    global p, V, I, b, c, m
    p, V, I, b, c, m = 2.3, 0.9, 1.2, 2.0, 1.5, 1.0
    global F, t
    u_e = lambda t: p*t**2 + V*t + I 	# compute with p, V, I, b, c, m as numbers
    F = ode_source_term(u_e, damping)   # fit source term
    F = sym.lambdify(t, F)            	# turn to numerical Python function

    dt = 2.*m/c
    u, t = solver(I=I, V=V, F=F, c=c, m=m, dt=dt, T=3, damping=damping)
    u_e = u_e(t)
    error = np.abs(u - u_e).max()
    nt.assert_almost_equal(error, 0, delta=1E-12)
    print 'Error in computing a quadratic solution:', error

if __name__ == '__main__':
    damping = ['zero', 'linear', 'nonlinear']
    for e in damping:
        V, t, I, dt, m, b, c = sym.symbols('V t I dt m b c')  # global symbols
        F = None  # global variable for the source term in the ODE
        print '---------------------------------------Damping:', e
        linear(e)  	# linear solution used for MMS
        quadratic(e)   	# quadratic solution for MMS
        cubic(e)       	# ... and cubic
        test_quadratic_exact_solution(e)
