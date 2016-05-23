import os, sys
sys.path.insert(0, os.path.join(os.pardir, 'src-vib'))
from vib_undamped import solver, u_exact, visualize
import numpy as np

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
    energy_const = 0.5*I**2*w**2    # initial energy when V = 0

    dt_values = []
    E_u_values = []         # error in u
    E_energy_values = []    # error in energy
    for i in range(m):
        u, t = solver_function(I, w, dt, T)
        u_e = u_exact(t, I, w)
        E_u = np.sqrt(dt*np.sum((u_e-u)**2))
        E_u_values.append(E_u)
        energy = 0.5*((u[2:] - u[:-2])/(2*dt))**2 + \
                                    0.5*w**2*u[1:-1]**2
        E_energy = energy - energy_const
        E_energy_norm = np.abs(E_energy).max()
        E_energy_values.append(E_energy_norm)
        dt_values.append(dt)
        dt = dt/2

    r_u = [np.log(E_u_values[i-1]/E_u_values[i])/
         np.log(dt_values[i-1]/dt_values[i])
         for i in range(1, m, 1)]
    r_E = [np.log(E_energy_values[i-1]/E_energy_values[i])/
         np.log(dt_values[i-1]/dt_values[i])
         for i in range(1, m, 1)]
    return r_u, r_E

def test_convergence_rates():
    r_u, r_E = convergence_rates(
        m=5,
        solver_function=solver,
        num_periods=8)
    # Accept rate to 1 decimal place
    tol = 0.1
    assert abs(r_u[-1] - 2.0) < tol
    assert abs(r_E[-1] - 2.0) < tol

if __name__ == '__main__':
    test_convergence_rates()
