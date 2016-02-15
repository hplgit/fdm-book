import numpy as np

def FE_logistic(p, u0, dt, Nt):
    u = np.zeros(Nt+1)
    u[0] = u0
    for n in range(Nt):
        u[n+1] = u[n] + dt*(1 - u[n])**p*u[n]
    return u

def BE_logistic(p, u0, dt, Nt, choice='Picard',
                eps_r=1E-3, omega=1, max_iter=1000):
    # u[n] = u[n-1] + dt*(1-u[n])**p*u[n]
    # -dt*(1-u[n])**p*u[n] + u[n] = u[n-1]
    if choice == 'Picard1':
        choice = 'Picard'
        max_iter = 1

    u = np.zeros(Nt+1)
    iterations = []
    u[0] = u0
    for n in range(1, Nt+1):
        c = -u[n-1]
        if choice == 'Picard':
            def F(u):
                return -dt*(1-u)**p*u + u + c

            u_ = u[n-1]
            k = 0
            while abs(F(u_)) > eps_r and k < max_iter:
                # u*(1-dt*(1-u_)**p) + c = 0
                u_ = omega*(-c/(1-dt*(1-u_)**p)) + (1-omega)*u_
                k += 1
            u[n] = u_
            iterations.append(k)

        elif choice == 'Newton':
            def F(u):
                return -dt*(1-u)**p*u + u + c

            def dF(u):
                return dt*p*(1-u)**(p-1)*u - dt*(1-u)**p + 1

            u_ = u[n-1]
            k = 0
            while abs(F(u_)) > eps_r and k < max_iter:
                u_ = u_ - F(u_)/dF(u_)
                k += 1
            u[n] = u_
            iterations.append(k)
    return u, iterations

def CN_logistic(p, u0, dt, Nt):
    # u[n+1] = u[n] + dt*(1-u[n])**p*u[n+1]
    # (1 - dt*(1-u[n])**p)*u[n+1] = u[n]
    u = np.zeros(Nt+1)
    u[0] = u0
    for n in range(0, Nt):
        u[n+1] = u[n]/(1 - dt*(1 - u[n])**p)
    return u

def test_asymptotic_value():
    T = 100
    dt = 0.1
    Nt = int(round(T/float(dt)))
    u0 = 0.1
    p = 1.8

    u_CN = CN_logistic(p, u0, dt, Nt)
    u_BE_Picard, iter_Picard = BE_logistic(
        p, u0, dt, Nt, choice='Picard',
        eps_r=1E-5, omega=1, max_iter=1000)
    u_BE_Newton, iter_Newton = BE_logistic(
        p, u0, dt, Nt, choice='Newton',
        eps_r=1E-5, omega=1, max_iter=1000)
    u_FE = FE_logistic(p, u0, dt, Nt)

    for arr in u_CN, u_BE_Picard, u_BE_Newton, u_FE:
        expected = 1
        computed = arr[-1]
        tol = 0.01
        msg = 'expected=%s, computed=%s' % (expected, computed)
        print msg
        assert abs(expected - computed) < tol

from scitools.std import *

def demo():
    T = 12
    p = 1.2
    try:
        dt = float(sys.argv[1])
        eps_r = float(sys.argv[2])
        omega = float(sys.argv[3])
    except:
        dt = 0.8
        eps_r = 1E-3
        omega = 1
    N = int(round(T/float(dt)))

    u_FE = FE_logistic(p, 0.1, dt, N)
    u_BE31, iter_BE31 = BE_logistic(p, 0.1, dt, N,
                                    'Picard1', eps_r, omega)
    u_BE3, iter_BE3 = BE_logistic(p, 0.1, dt, N,
                                  'Picard', eps_r, omega)
    u_BE4, iter_BE4 = BE_logistic(p, 0.1, dt, N,
                                  'Newton', eps_r, omega)
    u_CN = CN_logistic(p, 0.1, dt, N)

    print 'Picard mean no of iterations (dt=%g):' % dt, \
          int(round(mean(iter_BE3)))
    print 'Newton mean no of iterations (dt=%g):' % dt, \
          int(round(mean(iter_BE4)))

    t = np.linspace(0, dt*N, N+1)
    plot(t, u_FE, t, u_BE3, t, u_BE31, t, u_BE4, t, u_CN,
         legend=['FE', 'BE Picard', 'BE Picard1', 'BE Newton', 'CN gm'],
         title='dt=%g, eps=%.0E' % (dt, eps_r), xlabel='t', ylabel='u',
         legend_loc='lower right')
    filestem = 'logistic_N%d_eps%03d' % (N, log10(eps_r))
    savefig(filestem + '_u.png')
    savefig(filestem + '_u.pdf')
    figure()
    plot(range(1, len(iter_BE3)+1), iter_BE3, 'r-o',
         range(1, len(iter_BE4)+1), iter_BE4, 'b-o',
         legend=['Picard', 'Newton'],
         title='dt=%g, eps=%.0E' % (dt, eps_r),
         axis=[1, N+1, 0, max(iter_BE3 + iter_BE4)+1],
         xlabel='Time level', ylabel='No of iterations')
    savefig(filestem + '_iter.png')
    savefig(filestem + '_iter.pdf')
    raw_input()

def test_solvers():
    p = 2.5
    T = 5000
    dt = 0.5
    eps_r = 1E-6
    omega_values = [1]
    tol = 0.01
    N = int(round(T/float(dt)))

    for omega in omega_values:
        u_FE = FE_logistic(p, 0.1, dt, N)
        u_BE31, iter_BE31 = BE_logistic(p, 0.1, dt, N,
                                        'Picard1', eps_r, omega)
        u_BE3, iter_BE3 = BE_logistic(p, 0.1, dt, N,
                                      'Picard', eps_r, omega)
        u_BE4, iter_BE4 = BE_logistic(p, 0.1, dt, N,
                                      'Newton', eps_r, omega)
        u_CN = CN_logistic(p, 0.1, dt, N)

        print u_FE[-1], u_BE31[-1], u_BE3[-1], u_CN[-1]
        for u_x in u_FE, u_BE31, u_BE3, u_CN:
            print u_x[-1]
            assert abs(u_x[-1] - 1) < tol, 'u=%.16f' % u_x[-1]

    """
    t = np.linspace(0, dt*N, N+1)
    plot(t, u_FE, t, u_BE3, t, u_BE31, t, u_BE4, t, u_CN,
         legend=['FE', 'BE Picard', 'BE Picard1', 'BE Newton', 'CN gm'],
         title='dt=%g, eps=%.0E' % (dt, eps_r), xlabel='t', ylabel='u',
         legend_loc='lower right')
    filestem = 'tmp_N%d_eps%03d' % (N, log10(eps_r))
    savefig(filestem + '_u.png')
    savefig(filestem + '_u.pdf')
    """

if __name__ == '__main__':
    #demo()
    #test_solvers()
    test_asymptotic_value()
