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


def compute_orbit_and_error(
    f,
    solver_ID,
    timesteps_per_period=20,
    grouped_orbits=10,
    N_grouped_orbits=1000,
    E_measure='max_E'):
    '''
    For one particular solver:
    Calculte the orbits for a multiple of grouped orbits, i.e.
    number of orbits = N_grouped_orbits*grouped_orbits.
    Returns: time step dt, and, for each grouped_orbits cycle,
    the 2D position error and cpu time (as lists).
    '''
    def u_exact(t):
        return np.array([np.cos(t), np.sin(t)])

    w = 1
    P = 2*np.pi/w       # scaled period (1 year becomes 2*pi)
    dt = P/timesteps_per_period
    Nt = N_grouped_orbits*grouped_orbits*timesteps_per_period
    T = Nt*dt
    t_mesh = np.linspace(0, T, Nt+1)
    E_orbit = []
    CPU_time = []

    #print '        dt:', dt
    T_interval = P*grouped_orbits
    N = int(round(T_interval/dt))

    # set initial conditions
    if solver_ID == 'RK4':
        A = [1,0,0,1]
    elif solver_ID == 'EC':
        A = [0,1,1,0]
    elif solver_ID == 'PEFRL':
        I = np.array([1, 0])
        V = np.array([0, 1])
    else:
        print 'Unknown solver requested!'
        sys.exit(1)

    t1 = time.clock()
    for i in range(N_grouped_orbits):
        time_points = np.linspace(i*T_interval, (i+1)*T_interval, N+1)
        u_e = u_exact(time_points).transpose()
        if solver_ID == 'RK4':
            solver = odespy.RK4(f)
            solver.set_initial_condition(A)
            ui, ti = solver.solve(time_points)
            #find error (correct final pos:  x=1, y=0)
            if E_measure == 'final_E':
                orbit_error = np.sqrt((1-ui[-1,0])**2 +
                                      (0-ui[-1,2])**2)*100
            elif E_measure == 'max_E':
                orbit_error = np.sqrt((ui[:,0]-u_e[:,0])**2 +
                                      (ui[:,2]-u_e[:,1])**2).max()
            elif E_measure == 'all_E':
                orbit_error = np.sqrt(
                dt*np.sum((ui[:,0]-u_e[:,0])**2 +
                          (ui[:,2]-u_e[:,1])**2))
            else:
                print 'Unknown error measure requested!'
                sys.exit(1)
        elif solver_ID == 'EC':
            solver = odespy.EulerCromer(f)
            solver.set_initial_condition(A)
            ui, ti = solver.solve(time_points)
            #find error (correct final pos:  x=1, y=0)
            if E_measure == 'final_E':
                orbit_error = np.sqrt((1-ui[-1,1])**2 +
                                      (0-ui[-1,3])**2)*100
            elif E_measure == 'max_E':
                orbit_error = np.sqrt((ui[:,1]-u_e[:,0])**2 +
                                      (ui[:,3]-u_e[:,1])**2).max()
            elif E_measure == 'all_E':
                orbit_error = np.sqrt(
                    dt*np.sum((ui[:,1]-u_e[:,0])**2 +
                              (ui[:,3]-u_e[:,1])**2))
            else:
                print 'Unknown error measure requested!'
                sys.exit(1)
        else:   # solver_ID == 'PEFRL':
            # Note: every T_inverval is here counted from time 0
            ui, vi, ti = solver_PEFRL(I, V, f, dt, T_interval)
            #find error (correct final pos:  x=1, y=0)
            if E_measure == 'final_E':
                orbit_error = np.sqrt((1-ui[-1,0])**2 +
                                      (0-ui[-1,1])**2)*100
            elif E_measure == 'max_E':
                orbit_error = np.sqrt(
                    (ui[:,0]-u_e[:,0])**2 +
                    (ui[:,1]-u_e[:,1])**2).max()
            elif E_measure == 'all_E':
                orbit_error = np.sqrt(
                    dt*np.sum((ui[:,0]-u_e[:,0])**2 +
                              (ui[:,1]-u_e[:,1])**2))
            else:
                print 'Unknown error measure requested!'
                sys.exit(1)

        #print '      Orbit no. %d,   error (per cent): %g' % \
        #                   ((i+1)*grouped_orbits, orbit_error)

        E_orbit.append(orbit_error)
        t2 = time.clock()
        CPU_time.append((t2 - t1)/(60.0*60.0))    # in hours

        # set init. cond. for next time interval
        if solver_ID == 'RK4':
            A = [ui[-1,0], ui[-1,1], ui[-1,2], ui[-1,3]]
        elif solver_ID == 'EC':
            A = [ui[-1,0], ui[-1,1], ui[-1,2], ui[-1,3]]
        else:   # solver_ID == 'PEFRL':
            I = [ui[-1,0], ui[-1,1]]
            V = [vi[-1,0], vi[-1,1]]

    #print 'Total CPU time: ', CPU_time
    return dt, E_orbit, CPU_time

def orbit_error_vs_dt(
    f_EC, f_RK4, g, solvers,
    grouped_orbits=1000,    # orbits run between each error found
    N_grouped_orbits=10,
    E_measure='max_E'):
    '''
    With each solver in list "solvers":
    Simulate 10 000 orbits with different dt values.
    Collect final 2D position error for each dt and plot all errors.
    '''

    for solver_ID in solvers:
        print 'Computing orbit with solver:', solver_ID
        E_values = []
        dt_values = []
        cpu_values = []
        #for timesteps_per_period in 400, 800, 1600, 3200:
        for timesteps_per_period in \
                50, 100, 200, 400, 800, 1600, 3200:
            print '... steps per period: ', timesteps_per_period
            if solver_ID == 'RK4':
                dt, E, cpu_time = compute_orbit_and_error(
                    f_RK4,
                    solver_ID,
                    timesteps_per_period,
                    grouped_orbits,
                    N_grouped_orbits,
                    E_measure)
            elif solver_ID == 'EC':
                dt, E, cpu_time = compute_orbit_and_error(
                    f_EC,
                    solver_ID,
                    timesteps_per_period,
                    grouped_orbits,
                    N_grouped_orbits,
                    E_measure)
            elif solver_ID == 'PEFRL':
                dt, E, cpu_time = compute_orbit_and_error(
                    g,
                    solver_ID,
                    timesteps_per_period,
                    grouped_orbits,
                    N_grouped_orbits,
                    E_measure)
            else:
                print 'Unknown solver requested!'
                sys.exit(1)

            dt_values.append(dt)
            E_values.append(E[-1])  # need only after 10 000 cycles
            cpu_values.append(cpu_time[-1])
            #print 'CPU (in hours):', cpu_time
        #print 'E_values (10 000 years, changing dt):', E_values
        #print 'dt_values (10 000 years):', dt_values
        #print 'cpu_values (10 000 years, changing dt):', cpu_values
        #plt.figure()
        #plt.plot(dt_values, cpu_values, 'b*')
        #plt.xlabel('dt')
        #plt.ylabel('CPU (in hours) for 10000 years sim')
        #plt.title(solver_ID)
        #filename = solver_ID + '_CPU_after10000years_changing_dt'
        #plt.savefig(filename + '.png')
        #plt.savefig(filename + '.pdf')
        #plt.show()

        # Now make empirical formula: C*dt**r
        # i.e., we only estimate the conv. rate r

        E_values =  np.array(E_values)
        dt_values = np.array(dt_values)
        m = len(E_values)

        r = [np.log(E_values[i-1]/E_values[i])/
             np.log(dt_values[i-1]/dt_values[i])
             for i in range(1, m, 1)]
        print 'convergence rates (%s): %s' % (solver_ID, r)

        #tol = 0.1
        #if solver_ID == 'EC':   # should be 1st order
        #    assert abs(r[-1] - 1.0) < tol
        #else:   # RK4 and PEFRL should be 4th order
        #    assert abs(r[-1] - 4.0) < tol


def orbit_error_vs_years(
    f_EC, f_RK4, g, solvers,
    grouped_orbits=1000,
    N_grouped_orbits=100
    ):
    '''
    For each solver in the list solvers:
    simulate 100 000 orbits with a fixed dt corresponding to 1000
    steps per cycle (year). Collect 2D position errors for
    each 1000'th run, plot these errors and CPU. Finally, make an
    empirical formula for error and CPU as functions of a number
    of cycles.
    '''
    timesteps_per_period = 1000     # fixed for all runs

    for solver_ID in solvers:
        print 'Computing orbit with solver:', solver_ID
        if solver_ID == 'RK4':
            dt, E, cpu_time = compute_orbit_and_error(
                f_RK4,
                solver_ID,
                timesteps_per_period,
                grouped_orbits,
                N_grouped_orbits)
        elif solver_ID == 'EC':
            dt, E, cpu_time = compute_orbit_and_error(
                f_EC,
                solver_ID,
                timesteps_per_period,
                grouped_orbits,
                N_grouped_orbits)
        elif solver_ID == 'PEFRL':
            dt, E, cpu_time = compute_orbit_and_error(
                g,
                solver_ID,
                timesteps_per_period,
                grouped_orbits,
                N_grouped_orbits)
        else:
            print 'Unknown solver requested!'
            sys.exit(1)

        # E and cpu_time are for every grouped_orbits cycle
        print 'Error after 100 000 years:', E[-1]
        print 'CPU (in hours, 100 000 years):', cpu_time[-1]
        print 'E_values (fixed dt, changing no of years):', E
        print 'cpu_values:', cpu_time
        years = range(0, N_grouped_orbits)
        plt.plot(years, cpu_time, 'b*')
        plt.xlabel('Number of years')
        plt.ylabel('CPU (in hours, fixed dt)')
        plt.title(solver_ID)
        filename = solver_ID + '_CPU_during100000years_fixed_dt'
        plt.savefig(filename + '.png')
        plt.savefig(filename + '.pdf')
        plt.show()

        # Now make empirical formula

        #def E_of_years(x, p0, p1, p2, p3, p4):
        #    return p0*x**4 + p1*x**3 + p2*x**2 + p3*x + p4
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
        plt.xlabel('Number of years')
        plt.ylabel('orbit error (in per cent, fixed dt)')
        plt.title(solver_ID)
        filename = solver_ID + '_E_during100000years_fixed_dt'
        plt.savefig(filename + '.png')
        plt.savefig(filename + '.pdf')
        plt.show()

def compute_orbit_error_and_CPU(
    grouped_orbits=1000, N_grouped_orbits=10):
    '''
    Orbit error and associated CPU times are computed with
    some solvers: RK4, Euler-Cromer, PEFRL.'''

    def f_EC(u, t):
        '''
        Return derivatives for the 1st order system as
        required by Euler-Cromer.
        '''
        vx, x, vy, y = u  # u, array, holding [vx, x, vy, y]
        d = -(x**2 + y**2)**(-3.0/2)
        return [d*x, vx, d*y, vy ]

    def f_RK4(u, t):
        '''
        Return derivatives for the 1st order system as
        required by RK4.
        '''
        x, vx, y, vy = u  # u, array, holding [x, vx, y, vy]
        d = -(x**2 + y**2)**(-3.0/2)
        return [vx, d*x, vy, d*y ]

    def g(u, v):
        '''
        Return derivatives for the 1st order system as
        required by PEFRL.
        '''
        d = -(u[0]**2 + u[1]**2)**(-3.0/2)
        return np.array([d*u[0], d*u[1]])

    solvers = ['RK4', 'EC', 'PEFRL']
    #error_measure = 'final_E'
    error_measure = 'max_E'
    #error_measure = 'all_E'

    print 'Find orbit error as fu. of dt...'
    orbit_error_vs_dt(
        f_EC, f_RK4, g, solvers,
        grouped_orbits=grouped_orbits,
        N_grouped_orbits=N_grouped_orbits,
        E_measure=error_measure)

    #print 'Compute orbit error as fu. of no of years (fixed dt)...'
    #orbit_error_vs_years(
    #    f_EC, f_RK4, g, solvers,
    #    grouped_orbits=grouped_orbits,
    #    N_grouped_orbits=N_grouped_orbits)

if __name__ == '__main__':

    #test_solver_PEFRL()
    #compute_orbit_error_and_CPU(10000, 10)
    compute_orbit_error_and_CPU(500, 10)
    #compute_orbit_error_and_CPU(10, 1)
    #compute_orbit_error_and_CPU(100, 1)
    #compute_orbit_error_and_CPU(1000, 1)
