"""Verify exact solution of vib_undamped.solver function."""
from vib_undamped import solver
from numpy import arcsin as asin, pi, cos, abs
from scitools.std import plot

def test_solver_exact_discrete_solution():
    def tilde_w(w, dt):
        return (2./dt)*asin(w*dt/2.)

    def u_numerical_exact(t):
        return I*cos(tilde_w(w, dt)*t)

    w = 2.5
    I = 1.5

    # Estimate period and time step
    P = 2*pi/w
    num_periods = 4
    T = num_periods*P
    N = 5               # time steps per period
    dt = P/N
    u, t = solver(I, w, dt, T)
    u_e = u_numerical_exact(t)
    diff = abs(u_e - u).max()
    # Make a plot in a file, but not on the screen
    plot(t, u, 'bo', t, u_e, 'r-',
         legend=('numerical', 'exact'), show=False, savefig='tmp.png')

    assert diff < 1E-14

test_solver_exact_discrete_solution()
