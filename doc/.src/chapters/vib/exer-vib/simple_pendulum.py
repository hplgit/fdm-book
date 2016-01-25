import os, sys
sys.path.insert(0, os.path.join(os.pardir, 'src-vib'))
from vib import solver
import numpy as np
import matplotlib.pyplot as plt

def simulate(Theta, alpha, num_periods=10):
    # Dimensionless model requires the following parameters:
    from math import sin, pi

    I = Theta
    V = 0
    m = 1
    b = alpha
    s = lambda u: sin(u)
    F = lambda t: 0
    damping = 'quadratic'

    # Estimate T and dt from the small angle solution
    P = 2*pi   # One period (theta small, no drag)
    dt = P/40  # 40 intervals per period
    T = num_periods*P

    theta, t =  solver(I, V, m, b, s, F, dt, T, damping)
    omega = np.zeros(theta.size)
    omega[1:-1] = (theta[2:] - theta[:-2])/(2*dt)
    omega[0] = (theta[1] - theta[0])/dt
    omega[-1] = (theta[-1] - theta[-2])/dt

    S = omega**2 + np.cos(theta)
    D = alpha*np.abs(omega)*omega
    return t, theta, S, D

def visualize(t, theta, S, D, filename='tmp'):
    f, (ax1, ax2, ax3) = plt.subplots(3, sharex=True, sharey=False)
    ax1.plot(t, theta)
    ax1.set_title(r'$\theta(t)$')
    ax2.plot(t, S)
    ax2.set_title(r'Dimensonless force in the wire')
    ax3.plot(t, D)
    ax3.set_title(r'Dimensionless drag force')
    plt.savefig('%s.png' % filename)
    plt.savefig('%s.pdf' % filename)

import math
# Rough verification that small theta and no drag gives cos(t)
Theta = 1.0
alpha = 0
t, theta, S, D = simulate(Theta, alpha, num_periods=4)
# Scale theta by Theta (easier to compare with cos(t))
theta /= Theta
visualize(t, theta, S, D, filename='pendulum_verify')

Theta = math.radians(40)
alpha = 0.8
t, theta, S, D = simulate(Theta, alpha)
visualize(t, theta, S, D, filename='pendulum_alpha0.8_Theta40')
plt.show()
