import scitools.std as plt
import sys
import odespy
import numpy as np
import time

def solver_PEFRL(I, V, g, dt, T):
    """
    Solve v' = - g(u,v), u'=v for t in (0,T], u(0)=I and v(0)=V,
    by the PEFRL method.
    """
    dt = float(dt)
    Nt = int(round(T/dt))
    u = np.zeros((Nt+1, len(I)))
    v = np.zeros((Nt+1, len(I)))
    t = np.linspace(0, Nt*dt, Nt+1)

    # these values are from eq (20), ref to paper below
    xi = 0.1786178958448091
    lambda_ = -0.2123418310626054
    chi = -0.06626458266981849

    v[0] = V
    u[0] = I
    # Compare with eq 22 in http://arxiv.org/pdf/cond-mat/0110585.pdf
    for n in range(0, Nt):
        u_ = u[n] + xi*dt*v[n]
        v_ = v[n] + 0.5*(1-2*lambda_)*dt*g(u_, v[n])
        u_ = u_ + chi*dt*v_
        v_ = v_ + lambda_*dt*g(u_, v_)
        u_ = u_ + (1-2*(chi+xi))*dt*v_
        v_ = v_ + lambda_*dt*g(u_, v_)
        u_ = u_ + chi*dt*v_
        v[n+1] = v_ + 0.5*(1-2*lambda_)*dt*g(u_, v_)
        u[n+1] = u_ + xi*dt*v[n+1]
        #print 'v[%d]=%g, u[%d]=%g' % (n+1,v[n+1],n+1,u[n+1])
    return u, v, t

def test_solver_PEFRL():
    """Check 4th order convergence rate, using u'' + u = 0,
    I = 3.0, V = 0, which has the exact solution u_e = 3*cos(t)"""
    def g(u, v):
        return np.array([-u])
    def u_exact(t):
        return np.array([3*np.cos(t)]).transpose()
    I = u_exact(0)
    V = np.array([0])
    print 'V:', V, 'I:', I

    # Numerical parameters
    w = 1
    P = 2*np.pi/w
    dt_values = [P/20, P/40, P/80, P/160, P/320]
    T = 8*P
    error_vs_dt = []
    for n, dt in enumerate(dt_values):
        u, v, t = solver_PEFRL(I, V, g, dt, T)
        error = np.abs(u - u_exact(t)).max()
        print 'error:', error
        if n > 0:
            error_vs_dt.append(error/dt**4)
    for i in range(1, len(error_vs_dt)):
        #print abs(error_vs_dt[i]- error_vs_dt[0])
        assert abs(error_vs_dt[i]-
                   error_vs_dt[0]) < 0.1


def compute_orbit_and_error(f, solver_ID, timesteps_per_period=20,
                            num_1000periods = 100):
    '''
    For one particular solver:
    Calculte the orbit for a multiple of 1000 orbits.
    Returns: time step dt, and, for each 1000'th cycle, the 2D
    position error and cpu time (as lists).
    '''

    w = 1
    P = 2*np.pi/w       # scaled period (1 year becomes 2*pi)
    dt = P/timesteps_per_period
    Nt = num_1000periods*1000*timesteps_per_period
    T = Nt*dt
    t_mesh = np.linspace(0, T, Nt+1)
    E_orbit = []
    CPU_time = []

    # Each iteration: simulate 1000 periods and then check error
    print '        dt:', dt
    T_interval = P*1000
    N = int(round(T_interval/dt))

    # set initial conditions
    if solver_ID == 'RK4':
        A = [1,0,0,1]
    elif solver_ID == 'PEFRL':
        I = np.array([1, 0])
        V = np.array([0, 1])
    else:
        print 'Unknown solver requested!'
        sys.exit(1)

    t1 = time.clock()
    for i in range(num_1000periods):
        if solver_ID == 'RK4':
            time_points = np.linspace(i*T_interval, (i+1)*T_interval, N+1)
            solver = odespy.RK4(f)
            solver.set_initial_condition(A)
            ui, ti = solver.solve(time_points)
            #find deviation after 1000 next orbits (correct pos:  x=1, y=0)
            orbit_error = np.sqrt( (1-ui[-1,0])**2 + (0-ui[-1,2])**2)*100
        else:   # solver_ID == 'PEFRL':
            # Note: every T_inverval is here counted from time 0
            ui, vi, ti = solver_PEFRL(I, V, f, dt, T_interval)
            #find deviation after 1000 next orbits (correct pos:  x=1, y=0)
            orbit_error = np.sqrt( (1-ui[-1,0])**2 + (0-ui[-1,1])**2)*100

        print '        Orbit no. %d,   error (per cent): %g' % \
                                        ((i+1)*1000, orbit_error)
        E_orbit.append(orbit_error)
        t2 = time.clock()
        CPU_time.append((t2 - t1)/(60.0*60.0))    # in hours

        # set init. cond. for next time interval
        if solver_ID == 'RK4':
            A = [ui[-1,0], ui[-1,1], ui[-1,2], ui[-1,3]]
        else:   # solver_ID == 'PEFRL':
            I = [ui[-1,0], ui[-1,1]]
            V = [vi[-1,0], vi[-1,1]]

    #print 'Total CPU time: ', CPU_time
    return dt, E_orbit, CPU_time

def orbit_error_with_dt_10000_years_fixed(f, g, solvers):
    '''
    With each solver in list "solvers":
    Simulate 10 000 orbits with different dt values.
    Collect final 2D position error for each dt and plot all errors.
    '''

    num_1000periods = 10
    num_1000periods = 1
    #timestep_decrease = 50

    for solver_ID in solvers:
        print 'Computing orbit with solver:', solver_ID
        E_values = []
        dt_values = []
        cpu_values = []
        #timesteps_per_period = 500
        #for i in range(10):
        #for timesteps_per_period in 50, 100, 200, 400, 800, 1200, 1600:
        for timesteps_per_period in 400, 800, 1200, 1600:
            print '.......time steps per period: ', timesteps_per_period
            if solver_ID == 'RK4':
                dt, E, cpu_time = compute_orbit_and_error(
                                            f,
                                            solver_ID,
                                            timesteps_per_period,
                                            num_1000periods)
            elif solver_ID == 'PEFRL':
                dt, E, cpu_time = compute_orbit_and_error(
                                            g,
                                            solver_ID,
                                            timesteps_per_period,
                                            num_1000periods)
            else:
                print 'Unknown solver requested!'
                sys.exit(1)

            print 'CPU (in hours):', cpu_time
            dt_values.append(dt)
            E_values.append(E[-1])  # need only after 10 000 cycles
            cpu_values.append(cpu_time[-1])
            #timesteps_per_period -= timestep_decrease
        print 'E_values (10 000 years, changing dt):', E_values
        print 'dt_values (10 000 years):', dt_values
        print 'cpu_values (10 000 years, changing dt):', cpu_values
        plt.figure()
        plt.plot(dt_values, cpu_values, 'b*')
        plt.xlabel('dt')
        plt.ylabel('CPU (in hours) for 10000 years sim')
        plt.title(solver_ID)
        filename = solver_ID + '_CPU_after10000years_changing_dt'
        plt.savefig(filename + '.png')
        plt.savefig(filename + '.pdf')
        plt.show()

        # Now make empirical formula

        def E_of_dt(x, p0, p1, p2, p3, p4):
            return p0*x**4 + p1*x**3 + p2*x**2 + p3*x + p4
        # More general
        def E_of_dt(x, *coeff):
            return sum(coeff[i]*x*i for i in range(len(coeff)))
        E_values =  np.array(E_values)
        dt_values = np.array(dt_values)
        degree = 4
        # BYGG OM TIL r comp
        p = np.polyfit(dt_values, E_values, degree+1)
        p_str = map(str, p)
        formula = ' + '.join([p_str[i] + '*x**' + str(i) for i in range(degree)])
        print 'Empirical formula (E with dt, 10000 years):  ', formula
        plt.figure()
        plt.plot(dt_values,
                 E_values, 'b*',
                 dt_values,
                 E_of_dt(dt_values, *p), 'r--')
        plt.xlabel('dt')
        plt.ylabel('orbit error after 10000 years (per cent)')
        plt.title(solver_ID)
        filename = solver_ID + '_E_after10000years_changing_dt'
        plt.savefig(filename + '.png')
        plt.savefig(filename + '.pdf')
        plt.show()

def orbit_error_with_years_dt_fixed(f, g, solvers):
    '''
    For each solver in the list solvers:
    simulate 100 000 orbits with a fixed dt corresponding to 1000
    steps per cycle (year). Collect 2D position errors for
    each 1000'th run, plot these errors and CPU. Finally, make an
    empirical formula for error and CPU as functions of a number
    of cycles.
    '''
    timesteps_per_period = 1000     # fixed for all runs
    num_1000periods = 100
    num_1000periods = 10

    for solver_ID in solvers:
        print 'Computing orbit with solver:', solver_ID
        if solver_ID == 'RK4':
            dt, E, cpu_time = compute_orbit_and_error(
                f,
                solver_ID,
                timesteps_per_period,
                num_1000periods)
        elif solver_ID == 'PEFRL':
            dt, E, cpu_time = compute_orbit_and_error(
                g,
                solver_ID,
                timesteps_per_period,
                num_1000periods)
        else:
            print 'Unknown solver requested!'
            sys.exit(1)

        # E and cpu_time are lists with data for ev 1000th cycle
        print 'Error after 100 000 years:', E[-1]
        print 'CPU (in hours, 100 000 years):', cpu_time[-1]
        print 'E_values (fixed dt, changing no of years):', E
        print 'cpu_values:', cpu_time
        years = range(0, 100)
        plt.figure()
        plt.plot(cpu_time, 'b*')
        plt.xlabel('Number of 1000 periods')
        plt.ylabel('CPU (in hours, fixed dt)')
        plt.title(solver_ID)
        filename = solver_ID + '_CPU_during100000years_fixed_dt'
        plt.savefig(filename + '.png')
        plt.savefig(filename + '.pdf')
        plt.show()

        # Now make empirical formula

        def E_of_years(x, p0, p1, p2, p3, p4):
            return p0*x**4 + p1*x**3 + p2*x**2 + p3*x + p4
        def E_of_years(x, *coeff):
            return sum(coeff[i]*x*i for i in range(len(coeff)))
        E =  np.array(E)
        years = np.array(years)
        degree = 4
        p = np.polyfit(years, E, degree+1)
        p_str = map(str, p)
        formula = ' + '.join([p_str[i] + '*x**' + str(i) for i in range(degree)])
        #formula = p_str[0] + '*x**4 + ' + p_str[1] + '*x**3 + ' + \
        #      p_str[2] + '*x**2 + ' + p_str[3] + '*x + ' + p_str[4]
        print 'Empirical formula (E with years):  ', formula
        plt.figure()
        plt.plot(years,
                 E, 'b*',
                 years,
                 E_of_years(years, *p), 'r--')
        plt.xlabel('Number of years (in 1000 periods)')
        plt.ylabel('orbit error (in per cent, fixed dt)')
        plt.title(solver_ID)
        filename = solver_ID + '_E_during100000years_fixed_dt'
        plt.savefig(filename + '.png')
        plt.savefig(filename + '.pdf')
        plt.show()

def compute_orbit_error_and_CPU():
    '''Orbit error and CPU times are computed with RK4 and
    PEFRL solvers.'''

    def f(u, t):
        '''Return derivatives for the 1st order system'''
        x, vx, y, vy = u  # u, array, length 4, holding [x, vx, y, vy]
        d = -(x**2 + y**2)**(-3.0/2)
        return [vx, d*x, vy, d*y ]

    def g(u, v):
        d = -(u[0]**2 + u[1]**2)**(-3.0/2)
        return np.array([d*u[0], d*u[1]])

    solvers = ['RK4','PEFRL']

    print 'Find orbit error after 10 000 years as fu. of dt...'
    orbit_error_with_dt_10000_years_fixed(f, g, solvers)

    print 'Compute orbit error as fu. of no of years (fixed dt)...'
    orbit_error_with_years_dt_fixed(f, g, solvers)


if __name__ == '__main__':

    #test_solver_PEFRL()
    compute_orbit_error_and_CPU()
