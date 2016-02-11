import numpy as np
import matplotlib.pyplot as plt
from vib_undamped import solver, u_exact, visualize

def convergence_rates(m, solver_function, num_periods=8):
    """
    Return m-1 empirical estimates of the convergence rate
    based on m simulations, where the time step is halved
    for each simulation.
    solver_function(I, w, dt, T) solves each problem, where T
    is based on simulation for num_periods periods.
    """
    from math import pi
    w = 0.35; I = 0.3       # just chosen values
    P = 2*pi/w              # period
    dt = P/30               # 30 time step per period 2*pi/w
    T = P*num_periods
    Energy_const = 0.5*I**2*w**2    # initial energy when V = 0

    dt_values = []
    E_u_values = []         # error in u
    E_Energy_values = []    # error in energy
    for i in range(m):
        u, t = solver_function(I, w, dt, T)
        u_e = u_exact(t, I, w)
        E_u = np.sqrt(dt*np.sum((u_e-u)**2))
        E_u_values.append(E_u)
        Energy = 0.5*((u[2:] - u[:-2])/(2*dt))**2 + 0.5*w**2*u[1:-1]**2
        E_Energy = Energy - Energy_const
        E_Energy_norm = np.abs(E_Energy).max()
        E_Energy_values.append(E_Energy_norm)        
        dt_values.append(dt)
        dt = dt/2

    r_u = [np.log(E_u_values[i-1]/E_u_values[i])/
         np.log(dt_values[i-1]/dt_values[i])
         for i in range(1, m, 1)]
    r_E = [np.log(E_Energy_values[i-1]/E_Energy_values[i])/
         np.log(dt_values[i-1]/dt_values[i])
         for i in range(1, m, 1)]        
    return r_u, r_E

def test_convergence_rates():
    r_u, r_E = convergence_rates(m=5, solver_function=solver, num_periods=8)
    # Accept rate to 1 decimal place
    tol = 0.1
    assert abs(r_u[-1] - 2.0) < tol
    assert abs(r_E[-1] - 2.0) < tol    

def main(solver_function=solver):
    import argparse
    from math import pi
    parser = argparse.ArgumentParser()
    parser.add_argument('--I', type=float, default=1.0)
    parser.add_argument('--w', type=float, default=2*pi)
    parser.add_argument('--dt', type=float, default=0.05)
    parser.add_argument('--num_periods', type=int, default=5)
    parser.add_argument('--savefig', action='store_true')
    # Hack to allow --SCITOOLS options (read when importing scitools.std)
    parser.add_argument('--SCITOOLS_easyviz_backend', default='matplotlib')
    a = parser.parse_args()
    I, w, dt, num_periods, savefig = \
       a.I, a.w, a.dt, a.num_periods, a.savefig
    P = 2*pi/w  # one period
    T = P*num_periods
    u, t = solver_function(I, w, dt, T)
    visualize(u, t, I, w)
    plt.show()


if __name__ == '__main__':
    #main()
    test_convergence_rates()