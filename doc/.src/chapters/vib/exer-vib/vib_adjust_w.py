from numpy import *
from matplotlib.pyplot import *

def solver(I, w, dt, T, adjust_w=True):
    """
    Solve u'' + w**2*u = 0 for t in (0,T], u(0)=I and u'(0)=0,
    by a central finite difference method with time step dt.
    """
    dt = float(dt)
    Nt = int(round(T/dt))
    u = zeros(Nt+1)
    t = linspace(0, Nt*dt, Nt+1)
    if adjust_w:
        w = w*(1 - 1./24*w**2*dt**2)

    u[0] = I
    u[1] = u[0] - 0.5*dt**2*w**2*u[0]
    for n in range(1, Nt):
        u[n+1] = 2*u[n] - u[n-1] - dt**2*w**2*u[n]
    return u, t

def exact_solution(t, I, w):
    return I*cos(w*t)

def visualize(u, t, I, w):
    plot(t, u, 'r--o')
    t_fine = linspace(0, t[-1], 1001)  # very fine mesh for u_e
    u_e = exact_solution(t_fine, I, w)
    hold('on')
    plot(t_fine, u_e, 'b-')
    legend(['numerical', 'exact'], loc='upper left')
    xlabel('t')
    ylabel('u')
    dt = t[1] - t[0]
    title('dt=%g' % dt)
    umin = -1.2*I;  umax = -umin
    axis([t[0], t[-1], umin, umax])
    savefig('tmp1.png'); savefig('tmp1.pdf')
    show()

def convergence_rates(m, num_periods=8, adjust_w=True):
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
    for i in range(m):
        u, t = solver(I, w, dt, T, adjust_w)
        u_e = exact_solution(t, I, w)
        E = sqrt(dt*sum((u_e-u)**2))
        dt_values.append(dt)
        E_values.append(E)
        dt = dt/2

    r = [log(E_values[i-1]/E_values[i])/
         log(dt_values[i-1]/dt_values[i])
         for i in range(1, m, 1)]
    return r

def test_convergence_rates():
    r = convergence_rates(m=5, num_periods=8)
    # Accept rough approximation to rate
    assert abs(r[-1] - 4.0) < 0.1

def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--adjust_w', type=str, default='yes')
    parser.add_argument('--I', type=float, default=1.0)
    parser.add_argument('--w', type=float, default=2*pi)
    parser.add_argument('--dt', type=float, default=0.05)
    parser.add_argument('--num_periods', type=int, default=5)
    # Hack to allow --SCITOOLS options
    # (read when importing scitools.std)
    parser.add_argument('--SCITOOLS_easyviz_backend',
                        default='matplotlib')
    a = parser.parse_args()
    adjust_w, I, w, dt, num_periods = \
            a.adjust_w, a.I, a.w, a.dt, a.num_periods
    adjust_w = True if adjust_w == 'yes' else False

    P = 2*pi/w  # one period
    T = P*num_periods
    u, t = solver(I, w, dt, T, adjust_w)
    if num_periods <= 10:
        visualize(u, t, I, w)
    else:
        visualize_front(u, t, I, w)
        #visualize_front_ascii(u, t, I, w)


def visualize_front(u, t, I, w, savefig=False):
    """
    Visualize u and the exact solution vs t, using a
    moving plot window and continuous drawing of the
    curves as they evolve in time.
    Makes it easy to plot very long time series.
    """
    import scitools.std as st
    from scitools.MovingPlotWindow import MovingPlotWindow

    P = 2*pi/w  # one period
    umin = -1.2*I;  umax = -umin
    plot_manager = MovingPlotWindow(
        window_width=8*P,
        dt=t[1]-t[0],
        yaxis=[umin, umax],
        mode='continuous drawing')
    for n in range(1,len(u)):
        if plot_manager.plot(n):
            s = plot_manager.first_index_in_plot
            st.plot(t[s:n+1], u[s:n+1], 'r-1',
                    t[s:n+1], I*cos(w*t)[s:n+1], 'b-1',
                    title='t=%6.3f' % t[n],
                    axis=plot_manager.axis(),
                    show=not savefig) # drop window if savefig
            if savefig:
                st.savefig('tmp_vib%04d.png' % n)
        plot_manager.update(n)

if __name__ == '__main__':
    main()
    #r = convergence_rates(m=5, num_periods=8, adjust_w=True)
    #print 'convergence rate: %.1f' % r[-1]
