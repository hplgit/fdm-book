"""Solve u'' + w*2**u = 0 by the Euler-Cromer method in Odespy."""

import scitools.std as plt
import odespy
import numpy as np

def f(u, t, w=1):
    v, u = u
    return [-w**2*u, v]

def run_solver_and_plot(solver, timesteps_per_period=20,
                        num_periods=1, I=1, w=2*np.pi):
    P = 2*np.pi/w  # duration of one period
    dt = P/timesteps_per_period
    Nt = num_periods*timesteps_per_period
    T = Nt*dt
    t_mesh = np.linspace(0, T, Nt+1)

    solver.set(f_kwargs={'w': w})
    solver.set_initial_condition([0, I])
    u, t = solver.solve(t_mesh)

    from vib_undamped import solver
    u2, t2 = solver(I, w, dt, T)

    plt.plot(t, u[:,1], 'r-', t2, u2, 'b-')
    plt.legend(['Euler-Cromer', '2nd-order ODE'])
    plt.xlabel('t');  plt.ylabel('u')
    plt.savefig('tmp1.png'); plt.savefig('tmp1.pdf')

run_solver_and_plot(odespy.EulerCromer(f),
                    timesteps_per_period=20,
                    num_periods=9)
raw_input()
