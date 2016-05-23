import numpy as np

def solver(H, C_R, dt, T, eps_v=0.01, eps_h=0.01):
    """
    Simulate bouncing ball until it comes to rest. Time step dt.
    h(0)=H (initial height). T: maximum simulation time.
    Method: Euler-Cromer.
    """
    dt = float(dt)
    Nt = int(round(T/dt))
    h = np.zeros(Nt+1)
    v = np.zeros(Nt+1)
    t = np.linspace(0, Nt*dt, Nt+1)
    g = 0.81

    v[0] = 0
    h[0] = H
    mode = 'free fall'
    for n in range(Nt):
        v[n+1] = v[n] - dt*g
        h[n+1] = h[n] + dt*v[n+1]

        if h[n+1] < eps_h:
            #if abs(v[n+1]) > eps_v:  # handles large dt, but is wrong
            if v[n+1] < -eps_v:
                # Impact
                v[n+1] = -C_R*v[n+1]
                h[n+1] = 0
                if mode == 'impact':
                    # impact twice
                    return h[:n+2], v[:n+2], t[:n+2]
                mode = 'impact'
            elif abs(v[n+1]) < eps_v:
                mode = 'rest'
                v[n+1] = 0
                h[n+1] = 0
                return h[:n+2], v[:n+2], t[:n+2]
            else:
                mode = 'free fall'
        else:
            mode = 'free fall'
        print '%4d v=%8.5f h=%8.5f %s' % (n, v[n+1], h[n+1], mode)
    raise ValueError('T=%g is too short simulation time' % T)

import matplotlib.pyplot as plt
h, v, t = solver(
    H=1, C_R=0.8, T=100, dt=0.0001, eps_v=0.01, eps_h=0.01)
plt.plot(t, h)
plt.legend('h')
plt.savefig('tmp_h.png'); plt.savefig('tmp_h.pdf')
plt.figure()
plt.plot(t, v)
plt.legend('v')
plt.savefig('tmp_v.png'); plt.savefig('tmp_v.pdf')
plt.show()
