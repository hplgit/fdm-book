import sys, os
sys.path.insert(0, os.path.join(os.pardir, 'src'))  # for import vib
import numpy as np
from math import pi

def solver(I, V, m, b, s, F, dt, T, damping='linear'):
    """
    Solve m*u'' + f(u') + s(u) = F(t) for t in (0,T], u(0)=I,
    u'(0)=V by an Euler-Cromer method.
    """
    f = lambda v: b*v if damping == 'linear' else b*abs(v)*v
    dt = float(dt)
    Nt = int(round(T/dt))
    u = np.zeros(Nt+1)
    v = np.zeros(Nt+1)
    t = np.linspace(0, Nt*dt, Nt+1)

    v[0] = V
    u[0] = I
    for n in range(0, Nt):
        v[n+1] = v[n] + dt*(1./m)*(F(t[n]) - s(u[n]) - f(v[n]))
        u[n+1] = u[n] + dt*v[n+1]
        #print 'F=%g, s=%g, f=%g, v_prev=%g' % (F(t[n]), s(u[n]), f(v[n]), v[n])
        #print 'v[%d]=%g u[%d]=%g' % (n+1,v[n+1],n+1,u[n+1])
    return u, v, t

def test_solver():
    """Check 1st order convergence rate."""
    m = 4; b = 0.1
    s = lambda u: 2*u
    f = lambda v: b*v

    import sympy as sym
    def ode(u):
        """Return source F(t) in ODE for given manufactured u."""
        print 'ode:', m*sym.diff(u, t, 2), f(sym.diff(u,t)), s(u)
        return m*sym.diff(u, t, 2) + f(sym.diff(u,t)) + s(u)

    t = sym.symbols('t')
    u = 3*sym.cos(t)
    F = ode(u)
    F = sym.simplify(F)
    print 'F:', F, 'u:', u
    F = sym.lambdify([t], F, modules='numpy')
    u_exact = sym.lambdify([t], u, modules='numpy')
    I = u_exact(0)
    V = sym.diff(u, t).subs(t, 0)
    print 'V:', V, 'I:', I

    # Numerical parameters
    w = np.sqrt(0.5)
    P = 2*pi/w
    dt_values = [P/20, P/40, P/80, P/160, P/320]
    T = 8*P
    error_vs_dt = []
    for n, dt in enumerate(dt_values):
        u, v, t = solver(I, V, m, b, s, F, dt, T, damping='linear')
        error = np.abs(u - u_exact(t)).max()
        if n > 0:
            error_vs_dt.append(error/dt)
    for i in range(len(error_vs_dt)):
        assert abs(error_vs_dt[i]-
                   error_vs_dt[0]) < 0.1

def demo():
    """
    Demonstrate difference between Euler-Cromer and the
    scheme for the corresponding 2nd-order ODE.
    """
    I = 1.2; V = 0.2; m = 4; b = 0.2
    s = lambda u: 2*u
    F = lambda t: 0
    w = np.sqrt(2./4)   # approx freq
    dt = 0.9*2/w  # longest possible time step
    w = 0.5
    P = 2*pi/w
    T = 4*P
    from vib import solver as solver2
    import scitools.std as plt
    for k in range(4):
        u2, t2 = solver2(I, V, m, b, s, F, dt, T, 'quadratic')
        u, v, t = solver(I, V, m, b, s, F, dt, T, 'quadratic')
        plt.figure()
        plt.plot(t, u, 'r-', t2, u2, 'b-')
        plt.legend(['Euler-Cromer', 'centered scheme'])
        plt.title('dt=%.3g' % dt)
        raw_input()
        plt.savefig('tmp_%d' % k + '.png')
        plt.savefig('tmp_%d' % k + '.pdf')
        dt /= 2

if __name__ == '__main__':
    test_solver()
    #demo()
