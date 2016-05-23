import numpy as np
#import matplotlib.pyplot as plt
import scitools.std as plt

def solver_bwdamping(I, V, m, b, s, F, dt, T, damping='linear'):
    """
    Solve m*u'' + f(u') + s(u) = F(t) for t in (0,T],
    u(0)=I and u'(0)=V. All terms except damping is discretized
    by a central finite difference method with time step dt.
    The damping term is discretized by a backward diff. approx.,
    as is the init.cond. u'(0). If damping is 'linear', f(u')=b*u,
    while if damping is 'quadratic', f(u')=b*u'*abs(u').
    F(t) and s(u) are Python functions.
    """
    dt = float(dt); b = float(b); m = float(m) # avoid integer div.
    Nt = int(round(T/dt))
    u = np.zeros(Nt+1)
    t = np.linspace(0, Nt*dt, Nt+1)

    u_original = np.zeros(Nt+1); u_original[0] = I  # for testing

    u[0] = I
    if damping == 'linear':
        u[1] = u[0] + dt*V + dt**2/m*(-b*V - s(u[0]) + F(t[0]))
    elif damping == 'quadratic':
        u[1] = u[0] + dt*V + \
               dt**2/m*(-b*V*abs(V) - s(u[0]) + F(t[0]))
    for n in range(1, Nt):
        if damping == 'linear':
            u[n+1] = (2 - dt*b/m)*u[n] + dt**2/m*(- s(u[n]) + \
                                F(t[n])) + (dt*b/m - 1)*u[n-1]
        elif damping == 'quadratic':
            u[n+1] = 2*u[n] - u[n-1] - b/m*abs(u[n] - \
            u[n-1])*(u[n] - u[n-1]) + dt**2/m*(-s(u[n]) + F(t[n]))
    return u, t


import sympy as sym

def test_constant():
    """Verify a constant solution."""
    u_exact = lambda t: I
    I = 1.2; V = 0; m = 2; b = 0.9
    w = 1.5
    s = lambda u: w**2*u
    F = lambda t: w**2*u_exact(t)
    dt = 0.2
    T = 2
    #u, t = solver(I, V, m, b, s, F, dt, T, 'linear')
    u, t = solver_bwdamping(I, V, m, b, s, F, dt, T, 'linear')
    difference = np.abs(u_exact(t) - u).max()
    print difference
    tol = 1E-13
    assert difference < tol

    #u, t = solver(I, V, m, b, s, F, dt, T, 'quadratic')
    u, t = solver_bwdamping(I, V, m, b, s, F, dt, T, 'quadratic')
    difference = np.abs(u_exact(t) - u).max()
    print difference
    assert difference < tol

def lhs_eq(t, m, b, s, u, damping='linear'):
    """Return lhs of differential equation as sympy expression."""
    v = sym.diff(u, t)
    if damping == 'linear':
        return m*sym.diff(u, t, t) + b*v + s(u)
    else:
        return m*sym.diff(u, t, t) + b*v*sym.Abs(v) + s(u)

def test_quadratic():
    """Verify a quadratic solution."""
    I = 1.2; V = 3; m = 2; b = 0.9
    s = lambda u: 4*u
    t = sym.Symbol('t')
    dt = 0.2
    T = 2

    q = 2  # arbitrary constant
    u_exact = I + V*t + q*t**2
    F = sym.lambdify(t, lhs_eq(t, m, b, s, u_exact, 'linear'))
    u_exact = sym.lambdify(t, u_exact, modules='numpy')
    #u1, t1 = solver(I, V, m, b, s, F, dt, T, 'linear')
    u1, t1 = solver_bwdamping(I, V, m, b, s, F, dt, T, 'linear')
    diff = np.abs(u_exact(t1) - u1).max()
    print diff
    tol = 1E-13
    #assert diff < tol

    # In the quadratic damping case, u_exact must be linear
    # in order to exactly recover this solution
    u_exact = I + V*t
    F = sym.lambdify(t, lhs_eq(t, m, b, s, u_exact, 'quadratic'))
    u_exact = sym.lambdify(t, u_exact, modules='numpy')
    #u2, t2 = solver(I, V, m, b, s, F, dt, T, 'quadratic')
    u2, t2 = solver_bwdamping(I, V, m, b, s, F, dt, T, 'quadratic')
    diff = np.abs(u_exact(t2) - u2).max()
    print diff
    assert diff < tol

def test_sinusoidal():
    """Verify a numerically exact sinusoidal solution when b=F=0."""
    from math import asin

    def u_exact(t):
        w_numerical = 2/dt*np.arcsin(w*dt/2)
        return I*np.cos(w_numerical*t)

    I = 1.2; V = 0; m = 2; b = 0
    w = 1.5  # fix the frequency
    s = lambda u: m*w**2*u
    F = lambda t: 0
    dt = 0.2
    T = 6
    #u, t = solver(I, V, m, b, s, F, dt, T, 'linear')
    u, t = solver_bwdamping(I, V, m, b, s, F, dt, T, 'linear')
    diff = np.abs(u_exact(t) - u).max()
    print diff
    tol = 1E-14
    #assert diff < tol

    #u, t = solver(I, V, m, b, s, F, dt, T, 'quadratic')
    u, t = solver_bwdamping(I, V, m, b, s, F, dt, T, 'quadratic')
    print diff
    diff = np.abs(u_exact(t) - u).max()
    assert diff < tol

def test_mms():
    """Use method of manufactured solutions."""
    m = 4.; b = 1
    w = 1.5
    t = sym.Symbol('t')
    u_exact = 3*sym.exp(-0.2*t)*sym.cos(1.2*t)
    I = u_exact.subs(t, 0).evalf()
    V = sym.diff(u_exact, t).subs(t, 0).evalf()
    u_exact_py = sym.lambdify(t, u_exact, modules='numpy')
    s = lambda u: u**3
    dt = 0.2
    T = 6
    errors_linear = []
    errors_quadratic = []
    # Run grid refinements and compute exact error
    for i in range(5):
        F_formula = lhs_eq(t, m, b, s, u_exact, 'linear')
        F = sym.lambdify(t, F_formula)
        #u1, t1 = solver(I, V, m, b, s, F, dt, T, 'linear')
        u1, t1 = solver_bwdamping(I, V, m, b, s,
                                  F, dt, T, 'linear')
        error = np.sqrt(np.sum((u_exact_py(t1) - u1)**2)*dt)
        errors_linear.append((dt, error))

        F_formula = lhs_eq(t, m, b, s, u_exact, 'quadratic')
        #print sym.latex(F_formula, mode='plain')
        F = sym.lambdify(t, F_formula)
        #u2, t2 = solver(I, V, m, b, s, F, dt, T, 'quadratic')
        u2, t2 = solver_bwdamping(I, V, m, b, s,
                                  F, dt, T, 'quadratic')
        error = np.sqrt(np.sum((u_exact_py(t2) - u2)**2)*dt)
        errors_quadratic.append((dt, error))
        dt /= 2
    # Estimate convergence rates
    tol = 0.05
    for errors in errors_linear, errors_quadratic:
        for i in range(1, len(errors)):
            dt, error = errors[i]
            dt_1, error_1 = errors[i-1]
            r = np.log(error/error_1)/np.log(dt/dt_1)
        # check r for final simulation with (final and) smallest dt
        # note that the method now is 1st order, i.e. r should
        # approach 1.0
        print r
        assert abs(r - 1.0) < tol

import os, sys
sys.path.insert(0, os.path.join(os.pardir, 'src-vib'))
from vib import (plot_empirical_freq_and_amplitude,
                 visualize_front, visualize_front_ascii,
                 minmax, periods, amplitudes,
                 solver as solver2)

def visualize(list_of_curves, legends, title='', filename='tmp'):
    """Plot list of curves: (u, t)."""
    for u, t in list_of_curves:
        plt.plot(t, u)
        plt.hold('on')
    plt.legend(legends)
    plt.xlabel('t')
    plt.ylabel('u')
    plt.title(title)
    plt.savefig(filename + '.png')
    plt.savefig(filename + '.pdf')
    plt.show()

def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--I', type=float, default=1.0)
    parser.add_argument('--V', type=float, default=0.0)
    parser.add_argument('--m', type=float, default=1.0)
    parser.add_argument('--b', type=float, default=0.0)
    parser.add_argument('--s', type=str, default='4*pi**2*u')
    parser.add_argument('--F', type=str, default='0')
    parser.add_argument('--dt', type=float, default=0.05)
    parser.add_argument('--T', type=float, default=20)
    parser.add_argument('--window_width', type=float, default=30.,
                        help='Number of periods in a window')
    parser.add_argument('--damping', type=str, default='linear')
    parser.add_argument('--savefig', action='store_true')
    # Hack to allow --SCITOOLS options
    # (scitools.std reads this argument at import)
    parser.add_argument('--SCITOOLS_easyviz_backend',
                        default='matplotlib')
    a = parser.parse_args()
    from scitools.std import StringFunction
    s = StringFunction(a.s, independent_variable='u')
    F = StringFunction(a.F, independent_variable='t')
    I, V, m, b, dt, T, window_width, savefig, damping = \
       a.I, a.V, a.m, a.b, a.dt, a.T, a.window_width, a.savefig, \
       a.damping

    # compute u by both methods and then visualize the difference
    u, t = solver2(I, V, m, b, s, F, dt, T, damping)
    u_bw, _ = solver_bwdamping(I, V, m, b, s, F, dt, T, damping)
    u_diff = u - u_bw

    num_periods = plot_empirical_freq_and_amplitude(u_diff, t)
    if num_periods <= 40:
        plt.figure()
        legends = ['1st-2nd order method',
                   '2nd order method',
                   '1st order method']
        visualize([(u_diff, t), (u, t), (u_bw, t)], legends)
    else:
        visualize_front(u_diff, t, window_width, savefig)
        #visualize_front_ascii(u_diff, t)
    plt.show()

if __name__ == '__main__':
    main()
    #test_constant()
    #test_sinusoidal()
    #test_mms()
    #test_quadratic()
    raw_input()
