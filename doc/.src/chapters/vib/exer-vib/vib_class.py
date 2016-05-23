# Reimplementation of vib.py using classes

import numpy as np
import scitools.std as plt
import sympy as sym
from vib import solver as vib_solver
from vib import visualize as vib_visualize
from vib import visualize_front as vib_visualize_front
from vib import visualize_front_ascii as vib_visualize_front_ascii
from vib import plot_empirical_freq_and_amplitude as \
                            vib_plot_empirical_freq_and_amplitude

class Vibration:
    '''
    Problem: m*u'' + f(u') + s(u) = F(t) for t in (0,T],
    u(0)=I and u'(0)=V. The problem is solved
    by a central finite difference method with time step dt.
    If damping is 'linear', f(u')=b*u, while if damping is
    'quadratic', f(u')=b*u'*abs(u'). Zero damping is achieved
    with b=0. F(t) and s(u) are Python functions.
    '''
    def __init__(self, I=1, V=0, m=1, b=0, damping='linear'):
        self.I = I; self.V = V; self.m = m; self.b=b;
        self.damping = damping
    def s(self, u):
        return u
    def F(self, t):
        '''Driving force. Zero implies free oscillations'''
        return 0

class Free_vibrations(Vibration):
    '''F(t) = 0'''
    def __init__(self, s=None, I=1, V=0, m=1, b=0, damping='linear'):
        Vibration.__init__(self, I=I, V=V, m=m, b=b, damping=damping)
        if s != None:
            self.s = s

class Forced_vibrations(Vibration):
    '''F(t)! = 0'''
    def __init__(self, F, s=None, I=1, V=0, m=1, b=0,
                 damping='linear'):
        Vibration.__init__(self, I=I, V=V, m=m, b=b,
                           damping=damping)
        if s != None:
            self.s = s
        self.F = F

class Solver:
    def __init__(self, dt=0.05, T=20):
        self.dt = dt; self.T = T

    def solve(self, problem):
        self.u, self.t = vib_solver(
            problem.I, problem.V,
            problem.m, problem.b,
            problem.s, problem.F,
            self.dt, self.T, problem.damping)

class Visualizer:
    def __init__(self, problem, solver, window_width, savefig):
        self.problem = problem; self.solver = solver
        self.window_width = window_width; self.savefig = savefig
    def visualize(self):
        u = self.solver.u; t = self.solver.t   # short forms
        num_periods = vib_plot_empirical_freq_and_amplitude(u, t)
        if num_periods <= 40:
            plt.figure()
            vib_visualize(u, t)
        else:
            vib_visualize_front(u, t, self.window_width, self.savefig)
            vib_visualize_front_ascii(u, t)
        plt.show()

def main():
    # Note: the reading of parameter values would better be done
    # from each relevant class, i.e. class Problem should read I, V,
    #  etc., while class Solver should read dt and T, and so on.
    # Consult, e.g., Langtangen: "A Primer on Scientific Programming",
    # App E.
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--I', type=float, default=1.0)
    parser.add_argument('--V', type=float, default=0.0)
    parser.add_argument('--m', type=float, default=1.0)
    parser.add_argument('--b', type=float, default=0.0)
    parser.add_argument('--s', type=str, default=None)
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
    if a.s != None:
        s = StringFunction(a.s, independent_variable='u')
    else:
        s = None
    F = StringFunction(a.F, independent_variable='t')

    if a.F == '0':  # free vibrations
        problem = Free_vibrations(s=s, I=a.I, V=a.V, m=a.m, b=a.b,
                                  damping=a.damping)
    else:   # forced vibrations
        problem = Forced_vibrations(lambda t: np.sin(t),
                                    s=s, I=a.I, V=a.V,
                                    m=a.m, b=a.b, damping=a.damping)

    solver = Solver(dt=a.dt, T=a.T)
    solver.solve(problem)

    visualizer = Visualizer(problem, solver,
                            a.window_width, a.savefig)
    visualizer.visualize()

if __name__ == '__main__':
    main()
