import numpy as np

def solver(dt, T, f, f_0, f_1):
    """
    Solve u'=f by the Forward Euler method and by ordinary and
    Strang splitting: f(u) = f_0(u) + f_1(u).
    """
    Nt = int(round(T/float(dt)))
    t = np.linspace(0, Nt*dt, Nt+1)
    u_FE = np.zeros(len(t))
    u_split1 = np.zeros(len(t))  # 1st-order splitting
    u_split2 = np.zeros(len(t))  # 2nd-order splitting
    u_split3 = np.zeros(len(t))  # 2nd-order splitting w/exact f_0

    # Set initial values
    u_FE[0] = 0.1
    u_split1[0] = 0.1
    u_split2[0] = 0.1
    u_split3[0] = 0.1

    for n in range(len(t)-1):
        # Forward Euler method
        u_FE[n+1] = u_FE[n] + dt*f(u_FE[n])

        # --- Ordinary splitting ---
        # First step
        u_s_n = u_split1[n]
        u_s = u_s_n + dt*f_0(u_s_n)
        # Second step
        u_ss_n = u_s
        u_ss = u_ss_n + dt*f_1(u_ss_n)
        u_split1[n+1] = u_ss

        # --- Strang splitting ---
        # First step
        u_s_n = u_split2[n]
        u_s = u_s_n + dt/2.*f_0(u_s_n)
        # Second step
        u_sss_n = u_s
        u_sss = u_sss_n + dt*f_1(u_sss_n)
        # Third step
        u_ss_n = u_sss
        u_ss = u_ss_n + dt/2.*f_0(u_ss_n)
        u_split2[n+1] = u_ss

        # --- Strang splitting using exact integrator for u'=f_0 ---
        # First step
        u_s_n = u_split3[n]
        u_s = u_s_n*np.exp(dt/2.)  # exact
        # Second step
        u_sss_n = u_s
        u_sss = u_sss_n + dt*f_1(u_sss_n)
        # Third step
        u_ss_n = u_sss
        u_ss = u_ss_n*np.exp(dt/2.)  # exact
        u_split3[n+1] = u_ss

    return u_FE, u_split1, u_split2, u_split3, t

def solver_compact(dt, T, f, f_0, f_1):
    """
    As solver, but shorter code in the splitting steps.
    """
    Nt = int(round(T/float(dt)))
    t = np.linspace(0, Nt*dt, Nt+1)
    u_FE = np.zeros(len(t))
    u_split1 = np.zeros(len(t))  # 1st-order splitting
    u_split2 = np.zeros(len(t))  # 2nd-order splitting
    u_split3 = np.zeros(len(t))  # 2nd-order splitting w/exact f_0

    # Set initial values
    u_FE[0] = 0.1
    u_split1[0] = 0.1
    u_split2[0] = 0.1
    u_split3[0] = 0.1

    for n in range(len(t)-1):
        # Forward Euler method
        u_FE[n+1] = u_FE[n] + dt*f(u_FE[n])

        # Ordinary splitting
        u_s = u_split1[n] + dt*f_0(u_split1[n])
        u_split1[n+1] = u_s + dt*f_1(u_s)

        # Strang splitting
        u_s = u_split2[n] + dt/2.*f_0(u_split2[n])
        u_sss = u_s + dt*f_1(u_s)
        u_split2[n+1] = u_sss + dt/2.*f_0(u_sss)

        # Strang splitting using exact integrator for u'=f_0
        u_s = u_split3[n]*np.exp(dt/2.)  # exact
        u_ss = u_s + dt*f_1(u_s)
        u_split3[n+1] = u_ss*np.exp(dt/2.)

    return u_FE, u_split1, u_split2, u_split3, t

def demo(dt=0.2):
    u_exact = lambda t: 1./(9*np.exp(-t) + 1)
    u_FE, u_split1, u_split2, u_split3, t = solver(
        dt, 8, f = lambda u: u*(1-u),
        f_0=lambda u: u, f_1=lambda u: -u**2)

    import matplotlib.pyplot as plt
    plt.plot(t, u_FE, 'r-', t, u_split1, 'b-',
             t, u_split2, 'g-', t, u_split3, 'y-',
             t, u_exact(t), 'k--')
    plt.legend(['no split', 'split', 'strang', r'strang w/exact $f_0$', 'exact'],
               loc='lower right')
    plt.xlabel('t'); plt.ylabel('u')
    plt.title('Time step: %g' % dt)
    plt.savefig('tmp1.png'); plt.savefig('tmp1.pdf')
    plt.show()

def test_solver():
    np.set_printoptions(precision=15)
    u_FE_expected = np.array([float(x) for x in list("""
[ 0.1                0.118              0.1388152          0.162724308049792
  0.189973329573694  0.220750022298569  0.255153912289319
  0.293163990955874  0.334607764028414  0.379136845684478
  0.426215265270258  0.47512642785443   0.525002688936174
  0.574877662045366  0.62375632919069   0.67069320338774   0.714865969451186
  0.755632492485546  0.792562898242672  0.825444288357041
  0.854261491392197]"""[2:-1].split())])
    u_split1_expected = np.array([float(x) for x in list("""
[ 0.1                0.11712            0.1365934768128    0.158538732137911
  0.183007734044179  0.209963633605659  0.239259958824966
  0.270625296155645  0.303657796722007  0.33783343550351   0.37253027072271
  0.407068029717088  0.440758773984993  0.472961259290703
  0.503130113545367  0.530851841841463  0.555862750949651
  0.578048082546307  0.597425498363755  0.614118436921094
  0.628325385390187]"""[2:-1].split())])
    u_split2_expected = np.array([float(x) for x in list("""
[ 0.1                0.118338           0.139461146546647  0.1635705540078
  0.190798102531391  0.221174981642529  0.254599657026743
  0.290810238700023  0.329367696455926  0.369656716957107  0.41090943878828
  0.452253464828953  0.492779955548098  0.531621845295344
  0.568028513268957  0.601423369535241  0.631435056657206
  0.657899755122731  0.68083880192866   0.700420209898537
  0.716913803147616]"""[2:-1].split())])
    u_split3_expected = np.array([float(x) for x in list("""
[ 0.1                0.119440558200865  0.142033597399576
  0.168033940732164  0.197614356585594  0.230823935778256
  0.267544980228041  0.307455512661905  0.350006879617614
  0.394426527229859  0.439753524322399  0.484908174583947  0.5287881185737
  0.570374606392027  0.608827962442147  0.643553318053156
  0.674226057210925  0.700777592993228  0.723351459147494
  0.742244162719702  0.757844497686122]"""[2:-1].split())])
    for func in solver, solver_compact:
        u_FE, u_split1, u_split2, u_split3, t = solver(
            dt=0.2, T=4, f=lambda u: u*(1-u),
            f_0=lambda u: u, f_1=lambda u: -u**2)
        for quantity in 'u_FE', 'u_split1', 'u_split2', 'u_split3':
            diff = np.abs(eval(quantity + '_expected') -
                          eval(quantity)).max()
            assert diff < 1E-14

if __name__ == '__main__':
    test_solver()
    demo(0.05)
