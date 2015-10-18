import sys, os
sys.path.insert(0, os.path.join(os.pardir, 'src-vib'))

from vib_undamped import *

# Make new convergence_rates and test_convergence_rates where
# we also compute the rate of the energy measure.

def convergence_rates(m, num_periods=8):
    """
    Return m-1 empirical estimates of the convergence rate
    based on m simulations, where the time step is halved
    for each simulation.
    """
    w = 0.35; I = 0.3
    dt = 2*pi/w/30  # 30 time step per period 2*pi/w
    T = 2*pi/w*num_periods
    dt_values = []
    E_values = []
    energy_error_values = []
    for i in range(m):
        u, t = solver(I, w, dt, T)
        u_e = exact_solution(t, I, w)
        E = sqrt(dt*sum((u_e-u)**2))
        energy = 0.5*((u[2:] - u[:-2])/(2*dt))**2 + 0.5*w**2*u[1:-1]**2
        energy0 = 0.5*w**2*I**2
        energy_error = energy - energy0
        energy_error_norm = np.abs(energy_error).max()

        dt_values.append(dt)
        E_values.append(E)
        energy_error_values.append(energy_error_norm)
        dt = dt/2

    r_u = [
        log(E_values[i-1]/E_values[i])/
        log(dt_values[i-1]/dt_values[i])
        for i in range(1, m, 1)]
    r_energy = [
        log(energy_error_values[i-1]/energy_error_values[i])/
        log(dt_values[i-1]/dt_values[i])
        for i in range(1, m, 1)]
    return r_u, r_energy

def test_convergence_rates():
    r_u, r_energy = convergence_rates(m=5, num_periods=8)
    # Accept rate to 1 decimal place
    nt.assert_almost_equal(r_u[-1], 2.0, places=1)
    nt.assert_almost_equal(r_energy[-1], 2.0, places=1)

if __name__ == '__main__':
    test_convergence_rates()
