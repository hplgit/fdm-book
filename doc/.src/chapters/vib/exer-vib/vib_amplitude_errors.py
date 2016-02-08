import scitools.std as plt
#import matplotlib.pyplot as plt
from vib_empirical_analysis import minmax, amplitudes
import sys
import odespy
import numpy as np

def f(u, t, w=1):
    # v, u numbering for EulerCromer to work well
    v, u = u  # u is array of length 2 holding our [v, u]
    return [-w**2*u, v]

def run_solvers_and_check_amplitudes(solvers, timesteps_per_period=20,
                                     num_periods=1, I=1, w=2*np.pi):
    P = 2*np.pi/w  # duration of one period
    dt = P/timesteps_per_period
    Nt = num_periods*timesteps_per_period
    T = Nt*dt
    t_mesh = np.linspace(0, T, Nt+1)

    file_name = 'Amplitudes'   # initialize filename for plot
    for solver in solvers:
        solver.set(f_kwargs={'w': w})
        solver.set_initial_condition([0, I])
        u, t = solver.solve(t_mesh)

        solver_name = \
               'CrankNicolson' if solver.__class__.__name__ == \
               'MidpointImplicit' else solver.__class__.__name__
        file_name = file_name + '_' + solver_name

        minima, maxima = minmax(t, u[:,0])
        a = amplitudes(minima, maxima)
        plt.plot(range(len(a)), a, '-', label=solver_name)
        plt.hold('on')

    plt.xlabel('Number of periods')
    plt.ylabel('Amplitude (absolute value)')
    plt.legend(loc='upper left')
    plt.savefig(file_name + '.png')
    plt.savefig(file_name + '.pdf')
    plt.show()


# Define different sets of experiments
solvers_CNB2 = [odespy.CrankNicolson(f, nonlinear_solver='Newton'),
                odespy.Backward2Step(f)]
solvers_RK34 = [odespy.RK3(f), 
                odespy.RK4(f)]
solvers_AB = [odespy.AdamsBashforth2(f), 
              odespy.AdamsBashforth3(f)]

if __name__ == '__main__':
    # Default values
    timesteps_per_period = 30
    solver_collection = 'CNB2'
    num_periods = 100
    # Override from command line
    try:
        # Example: python vib_undamped_odespy.py 30 RK34 50
        timesteps_per_period = int(sys.argv[1])
        solver_collection = sys.argv[2]
        num_periods = int(sys.argv[3])
    except IndexError:
        pass # default values are ok
    solvers = eval('solvers_' + solver_collection)  # list of solvers
    run_solvers_and_check_amplitudes(solvers,
                                     timesteps_per_period,
                                     num_periods)

