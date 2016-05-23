import numpy as np
import matplotlib.pyplot as plt
from vib_verify_mms import solver

def u_exact(t, I, V, A, f, c, m):
    """Found by solving mu'' + cu = F in Wolfram alpha."""
    k_1 = I
    k_2 = (V - A*2*np.pi*f/(c - 4*np.pi**2*f**2*m))*\
                                        np.sqrt(m/float(c))
    return A*np.sin(2*np.pi*f*t)/(c - 4*np.pi**2*f**2*m) + \
           k_2*np.sin(np.sqrt(c/float(m))*t) + \
           k_1*np.cos(np.sqrt(c/float(m))*t)

def convergence_rates(N, solver_function, num_periods=8):
    """
    Returns N-1 empirical estimates of the convergence rate
    based on N simulations, where the time step is halved
    for each simulation.
    solver_function(I, V, F, c, m, dt, T, damping) solves
    each problem, where T is based on simulation for
    num_periods periods.
    """

    def F(t):
        """External driving force"""
        return A*np.sin(2*np.pi*f*t)

    b, c, m = 0, 1.6, 1.3  # just some chosen values
    I = 0                  # init. cond. u(0)
    V = 0                  # init. cond. u'(0)
    A = 1.0                # amplitude of driving force
    f = 1.0                # chosen frequency of driving force
    damping = 'zero'

    P = 1/f
    dt = P/30              # 30 time step per period 2*pi/w
    T = P*num_periods

    dt_values = []
    E_values = []
    for i in range(N):
        u, t = solver_function(I, V, F, b, c, m, dt, T, damping)
        u_e = u_exact(t, I, V, A, f, c, m)
        E = np.sqrt(dt*np.sum((u_e-u)**2))
        dt_values.append(dt)
        E_values.append(E)
        dt = dt/2

    #plt.plot(t, u, 'b--', t, u_e, 'r-'); plt.grid(); plt.show()

    r = [np.log(E_values[i-1]/E_values[i])/
         np.log(dt_values[i-1]/dt_values[i])
         for i in range(1, N, 1)]
    print r
    return r

def test_convergence_rates():
    r = convergence_rates(
        N=5,
        solver_function=solver,
        num_periods=8)
    # Accept rate to 1 decimal place
    tol = 0.1
    assert abs(r[-1] - 2.0) < tol

if __name__ == '__main__':
    test_convergence_rates()
