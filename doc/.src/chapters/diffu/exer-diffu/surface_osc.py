import sys, os
sys.path.insert(0, os.path.join(os.pardir, 'src-diffu'))
from diffu1D_vc import viz

sol = []  # store solutions
for Nx, L in [[20, 4], [40, 8]]:
    dt = 0.1
    dx = float(L)/Nx
    D = dt/dx**2
    from math import pi, sin
    T = 2*pi*6
    from numpy import zeros
    a = zeros(Nx+1) + 0.5
    cpu, u_ = viz(
        I=lambda x: 0, a=a, L=L, Nx=Nx, D=D, T=T,
        umin=-1.1, umax=1.1, theta=0.5,
        u_L=lambda t: sin(t),
        u_R=0,
        animate=False, store_u=True)
    sol.append(u_)
    print 'computed solution for Nx=%d in [0,%g]' % (Nx, L)

print sol[0].shape
print sol[1].shape
import scitools.std as plt
counter = 0
for u0, u1 in zip(sol[0][2:], sol[1][2:]):
    x0 = sol[0][0]
    x1 = sol[1][0]
    plt.plot(x0, u0, 'r-', x1, u1, 'b-',
             legend=['short', 'long'],
             savefig='tmp_%04d.png' % counter,
             axis=[x1[0], x1[-1], -1.1, 1.1])
    counter += 1
