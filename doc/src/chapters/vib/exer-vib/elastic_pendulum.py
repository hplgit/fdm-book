import odespy
import numpy as np
import scitools.std as plt

def simulate(
    beta=0.9,                 # dimensionless parameter
    Theta=30,                 # initial angle in degrees
    epsilon=0,                # initial stretch of wire
    num_periods=6,            # simulate for num_periods
    time_steps_per_period=60, # time step resolution
    plot=True,                # make plots or not
    ):
    from math import sin, cos, pi
    Theta = Theta*np.pi/180  # convert to radians
    # Initial position and velocity
    # (we order the equations such that Euler-Cromer in odespy
    # can be used, i.e., vx, x, vy, y)
    ic = [0,                              # x'=vx
          (1 + epsilon)*sin(Theta),       # x
          0,                              # y'=vy
          1 - (1 + epsilon)*cos(Theta),   # y
          ]

    def f(u, t, beta):
        vx, x, vy, y = u
        L = np.sqrt(x**2 + (y-1)**2)
        h = beta/(1-beta)*(1 - beta/L)  # help factor
        return [-h*x, vx, -h*(y-1) - beta, vy]

    # Non-elastic pendulum (scaled similarly in the limit beta=1)
    # solution Theta*cos(t)
    P = 2*pi
    dt = P/time_steps_per_period
    T = num_periods*P
    omega = 2*pi/P

    time_points = np.linspace(
        0, T, num_periods*time_steps_per_period+1)

    solver = odespy.EulerCromer(f, f_args=(beta,))
    solver.set_initial_condition(ic)
    u, t = solver.solve(time_points)
    x = u[:,1]
    y = u[:,3]
    theta = np.arctan(x/(1-y))

    if plot:
        plt.figure()
        plt.plot(x, y, 'b-', title='Pendulum motion',
                 daspect=[1,1,1], daspectmode='equal',
                 axis=[x.min(), x.max(), 1.3*y.min(), 1])
        plt.savefig('tmp_xy.png')
        plt.savefig('tmp_xy.pdf')
        # Plot theta in degrees
        plt.figure()
        plt.plot(t, theta*180/np.pi, 'b-',
                 title='Angular displacement in degrees')
        plt.savefig('tmp_theta.png')
        plt.savefig('tmp_theta.pdf')
        if abs(Theta) < 10*pi/180:
            # Compare theta and theta_e for small angles (<10 degrees)
            theta_e = Theta*np.cos(omega*t)  # non-elastic scaled sol.
            plt.figure()
            plt.plot(t, theta, t, theta_e,
                     legend=['theta elastic', 'theta non-elastic'],
                     title='Elastic vs non-elastic pendulum, '\
                            'beta=%g' % beta)
            plt.savefig('tmp_compare.png')
            plt.savefig('tmp_compare.pdf')
        # Plot y vs x (the real physical motion)
    return x, y, theta, t

def test_equilibrium():
    """Test that starting from rest makes x=y=theta=0."""
    x, y, theta, t = simulate(
        beta=0.9, Theta=0, epsilon=0,
        num_periods=6, time_steps_per_period=10, plot=False)
    tol = 1E-14
    assert np.abs(x.max()) < tol
    assert np.abs(y.max()) < tol
    assert np.abs(theta.max()) < tol

def test_vertical_motion():
    beta = 0.9
    omega = np.sqrt(beta/(1-beta))
    # Find num_periods. Recall that P=2*pi for scaled pendulum
    # oscillations, while here we don't have gravity driven
    # oscillations, but elastic oscillations with frequency omega.
    period = 2*np.pi/omega
    # We want T = N*period
    N = 5
    # simulate function has T = 2*pi*num_periods
    num_periods = 5/omega
    n = 600
    time_steps_per_period = omega*n

    y_exact = lambda t: -0.1*np.cos(omega*t)
    x, y, theta, t = simulate(
        beta=beta, Theta=0, epsilon=0.1,
        num_periods=num_periods,
        time_steps_per_period=time_steps_per_period,
        plot=False)

    tol = 0.00055 # ok tolerance for the above resolution
    # No motion in x direction is epxected
    assert np.abs(x.max()) < tol
    # Check motion in y direction
    y_e = y_exact(t)
    diff = np.abs(y_e - y).max()
    if diff > tol: # plot
        plt.plot(t, y, t, y_e, legend=['y', 'exact'])
        raw_input('Error in test_vertical_motion; type CR:')
    assert diff < tol, 'diff=%g' % diff

def demo(beta=0.999, Theta=40, num_periods=3):
    x, y, theta, t = simulate(
        beta=beta, Theta=Theta, epsilon=0,
        num_periods=num_periods, time_steps_per_period=600,
        plot=True)

if __name__ == '__main__':
    test_equilibrium()
    test_vertical_motion()
    #demo(0.999, num_periods=1)
    demo(0.93, num_periods=1)
    raw_input('Type CR: ')
