import numpy as np

def solver(dt, T, f, f_0, f_1):
    """
    Solve u'=f by the Forward Euler method and by ordinary and
    Strang splitting: f(u) = f_1(u) + f_2(u).
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
        u_0 = u_split1[n] + dt*f_0(u_split1[n])
        u_split1[n+1] = u_0 + dt*f_1(u_0)

        # Strang splitting
        u_0 = u_split2[n] + dt/2.*f_0(u_split2[n])
        u_1 = u_0 + dt*f_1(u_0)
        u_split2[n+1] = u_1 + dt/2.*f_0(u_1)

        # Strang splitting using exact integrator for u'=f_0
        u_0 = u_split3[n]*np.exp(dt/2.)  # exact
        u_1 = u_0 + dt*f_1(u_0)
        u_split3[n+1] = u_1*np.exp(dt/2.)

    return u_FE, u_split1, u_split2, u_split3, t

def demo():
    u_exact = lambda t: 1./(9*np.exp(-t) + 1)
    dt = 0.2
    u_FE, u_split1, u_split2, u_split3, t = solver(
        dt, 8, f = lambda u: u*(1-u),
        f_0=lambda u: u, f_1=lambda u: -u**2)

    import matplotlib.pyplot as plt
    plt.plot(t, u_FE, 'r-', t, u_split1, 'b-',
             t, u_split2, 'g-', t, u_split3, 'y-',
             t, u_exact(t), 'k--')
    plt.legend(['no split', 'split', 'strang', 'exact'],
               loc='lower right')
    plt.xlabel('t'); plt.ylabel('u')
    plt.title('Time step: %g' % dt)
    plt.savefig('tmp1.png'); plt.savefig('tmp1.pdf')
    plt.show()

if __name__ == '__main__':
    demo()
