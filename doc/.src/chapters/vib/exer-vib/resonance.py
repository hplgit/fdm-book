import sys, os
sys.path.insert(0, os.path.join(os.pardir, 'src-vib'))
from vib import solver, visualize, plt
from math import pi, sin
import numpy as np

beta_values = [0.005, 0.05, 0.2]
beta_values = [0.00005]
gamma_values = [5, 1.5, 1.1, 1]
for i, beta in enumerate(beta_values):
    for gamma in gamma_values:
        u, t = solver(I=1, V=0, m=1, b=2*beta, s=lambda u: u,
                      F=lambda t: sin(gamma*t), dt=2*pi/60,
                      T=2*pi*20, damping='quadratic')
        visualize(u, t, title='gamma=%g' %
                  gamma, filename='tmp_%s' % gamma)
        print gamma, 'max u amplitude:', np.abs(u).max()
    for ext in 'png', 'pdf':
        cmd = 'doconce combine_images '
        cmd += ' '.join(['tmp_%s.' % gamma + ext
                         for gamma in gamma_values])
        cmd += ' resonance%d.' % (i+1) + ext
        os.system(cmd)
raw_input()
