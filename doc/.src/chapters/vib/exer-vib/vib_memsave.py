import numpy as np
import scitools.std as plt

def solve_and_store(filename, I, V, m, b, s, F, dt, T, damping='linear'):
    """
    Solve m*u'' + f(u') + s(u) = F(t) for t in (0,T], u(0)=I and u'(0)=V,
    by a central finite difference method with time step dt.
    If damping is 'linear', f(u')=b*u, while if damping is
    'quadratic', f(u')=b*u'*abs(u'). F(t) and s(u) are Python functions.
    The solution is written to file (filename).
    Naming convention: we use the name u for the new solution to be 
    computed, u_1 for the solution one time step prior to that and u_2
    for the solution two time steps prior to that.
    Returns min and max u values needed for subsequent plotting.
    """
    dt = float(dt); b = float(b); m = float(m) # avoid integer div.
    Nt = int(round(T/dt))
    outfile = open(filename, 'w')
    outfile.write('Time          Position\n')

    u_2 = I
    u_min = u_max = u_2
    outfile.write('%6.3f         %7.5f\n' % (0*dt, u_2))
    if damping == 'linear':
        u_1 = u_2 + dt*V + dt**2/(2*m)*(-b*V - s(u_2) + F(0*dt))
    elif damping == 'quadratic':
        u_1 = u_2 + dt*V + \
               dt**2/(2*m)*(-b*V*abs(V) - s(u_2) + F(0*dt))
    if u_1 < u_2:
        u_min = u_1
    else:  # either equal or u_1 > u_2
        u_max = u_1
    outfile.write('%6.3f         %7.5f\n' % (1*dt, u_1))
 
    for n in range(1, Nt):
        # compute  solution at next time step
        if damping == 'linear':
            u = (2*m*u_1 + (b*dt/2 - m)*u_2 +
                dt**2*(F(n*dt) - s(u_1)))/(m + b*dt/2)
        elif damping == 'quadratic':
            u = (2*m*u_1 - m*u_2 + b*u_1*abs(u_1 - u_2)
                + dt**2*(F(n*dt) - s(u_1)))/\
                (m + b*abs(u_1 - u_2))
        if u < u_min:
            u_min = u
        elif u > u_max:
            u_max = u

        # write solution to file
        outfile.write('%6.3f         %7.5f\n' % ((n+1)*dt, u))
        # switch references before next step
        u_2, u_1, u = u_1, u, u_2

    outfile.close()
    return u_min, u_max

def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--I', type=float, default=1.0)
    parser.add_argument('--V', type=float, default=0.0)
    parser.add_argument('--m', type=float, default=1.0)
    parser.add_argument('--b', type=float, default=0.0)
    parser.add_argument('--s', type=str, default='u')
    parser.add_argument('--F', type=str, default='0')
    parser.add_argument('--dt', type=float, default=0.05)
    parser.add_argument('--T', type=float, default=10)
    parser.add_argument('--window_width', type=float, default=30.,
                        help='Number of periods in a window')
    parser.add_argument('--damping', type=str, default='linear')
    parser.add_argument('--savefig', action='store_true')
    # Hack to allow --SCITOOLS options (scitools.std reads this argument
    # at import)
    parser.add_argument('--SCITOOLS_easyviz_backend', default='matplotlib')
    a = parser.parse_args()
    from scitools.std import StringFunction
    s = StringFunction(a.s, independent_variable='u')
    F = StringFunction(a.F, independent_variable='t')
    I, V, m, b, dt, T, window_width, savefig, damping = \
       a.I, a.V, a.m, a.b, a.dt, a.T, a.window_width, a.savefig, \
       a.damping

    filename = 'vibration_sim.dat'
    u_min, u_max = solve_and_store(filename, I, V, m, b, s, F, dt, T, damping)

    read_and_plot(filename, u_min, u_max)

def read_and_plot(filename, u_min, u_max):
    """
    Read file and plot u vs t line by line in a
    terminal window (only using ascii characters).
    """
    from scitools.avplotter import Plotter
    import time
    umin = 1.2*u_min;  umax = 1.2*u_max
    p = Plotter(ymin=umin, ymax=umax, width=60, symbols='+o')
    fps = 10
    infile = open(filename, 'r')

    # read and treat one line at a time
    infile.readline()   # skip header line
    for line in infile:
        time_and_pos = line.split()  # gives list with 2 elements
        t = float(time_and_pos[0])
        u = float(time_and_pos[1])
        #print 'time: %g   position: %g' % (time, pos)
        print p.plot(t, u), '%.2f' % (t)
        time.sleep(1/float(fps))

if __name__ == '__main__':
    main()
