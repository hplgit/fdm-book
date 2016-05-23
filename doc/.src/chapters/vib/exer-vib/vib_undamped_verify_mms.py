import sympy as sym
import numpy as np

V, t, I, w, dt = sym.symbols('V t I w dt')  # global symbols
f = None  # global variable for the source term in the ODE

def ode_source_term(u):
    """Return the terms in the ODE that the source term
    must balance, here u'' + w**2*u.
    u is symbolic Python function of t."""
    return sym.diff(u(t), t, t) + w**2*u(t)

def residual_discrete_eq(u):
    """Return the residual of the discrete eq. with u inserted."""
    R = DtDt(u, dt) + w**2*u(t) - f
    return sym.simplify(R)

def residual_discrete_eq_step1(u):
    """Return the residual of the discrete eq. at the first
    step with u inserted."""
    half = sym.Rational(1,2)
    R = u(t+dt) - I - dt*V - \
        half*dt**2*f.subs(t, 0) + half*dt**2*w**2*I
    R = R.subs(t, 0)  # t=0 in the rhs of the first step eq.
    return sym.simplify(R)

def DtDt(u, dt):
    """Return 2nd-order finite difference for u_tt.
    u is a symbolic Python function of t.
    """
    return (u(t+dt) - 2*u(t) + u(t-dt))/dt**2

def main(u):
    """
    Given some chosen solution u (as a function of t, implemented
    as a Python function), use the method of manufactured solutions
    to compute the source term f, and check if u also solves
    the discrete equations.
    """
    print '=== Testing exact solution: %s ===' % u(t)
    print "Initial conditions u(0)=%s, u'(0)=%s:" % \
          (u(t).subs(t, 0), sym.diff(u(t), t).subs(t, 0))

    # Method of manufactured solution requires fitting f
    global f  # source term in the ODE
    f = sym.simplify(ode_source_term(u))

    # Residual in discrete equations (should be 0)
    print 'residual step1:', residual_discrete_eq_step1(u)
    print 'residual:', residual_discrete_eq(u)


def linear():
    """Test linear function V*t+I: u(0)=I, u'(0)=V."""
    main(lambda t: V*t + I)

def quadratic():
    """Test quadratic function q*t**2 + V*t + I."""
    q = sym.Symbol('q')  # arbitrary constant in t**2 term
    u_e = lambda t: q*t**2 + V*t + I
    main(u_e)

def cubic():
    r, q = sym.symbols('r q')
    main(lambda t: r*t**3 + q*t**2 + V*t + I)

def solver(I, V, f, w, dt, T):
    """
    Solve u'' + w**2*u = f for t in (0,T], u(0)=I and u'(0)=V,
    by a central finite difference method with time step dt.
    f(t) is a callable Python function.
    """
    dt = float(dt)
    Nt = int(round(T/dt))
    u = np.zeros(Nt+1)
    t = np.linspace(0, Nt*dt, Nt+1)

    u[0] = I
    u[1] = u[0] - 0.5*dt**2*w**2*u[0] + 0.5*dt**2*f(t[0]) + dt*V
    for n in range(1, Nt):
        u[n+1] = 2*u[n] - u[n-1] - dt**2*w**2*u[n] + dt**2*f(t[n])
    return u, t

def test_quadratic_exact_solution():
    """Verify solver function via quadratic solution."""
    # Transform global symbolic variables to functions and numbers
    # for numerical computations
    global p, V, I, w
    p, V, I, w = 2.3, 0.9, 1.2, 1.5
    global f, t
    u_e = lambda t: p*t**2 + V*t + I # use p, V, I, w as numbers
    f = ode_source_term(u_e)         # fit source term
    f = sym.lambdify(t, f)           # make function numerical

    dt = 2./w
    u, t = solver(I=I, V=V, f=f, w=w, dt=dt, T=3)
    u_e = u_e(t)
    error = np.abs(u - u_e).max()
    tol = 1E-12
    assert error < tol
    print 'Error in computing a quadratic solution:', error

if __name__ == '__main__':
    linear()
    quadratic()
    cubic()
    test_quadratic_exact_solution()
