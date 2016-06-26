import sys, os
sys.path.insert(0, os.path.join(os.pardir, 'src-vib'))

import numpy as np
import matplotlib.pyplot as plt

def solver_memsave(I, w, dt, T, filename='tmp.dat'):
    """
    As vib_undamped.solver, but store only the last three
    u values in the implementation. The solution is written to
    file `tmp_memsave.dat`.
    Solve u'' + w**2*u = 0 for t in (0,T], u(0)=I and u'(0)=0,
    by a central finite difference method with time step dt.
    """
    dt = float(dt)
    Nt = int(round(T/dt))
    t = np.linspace(0, Nt*dt, Nt+1)
    outfile = open(filename, 'w')

    u_n = I
    outfile.write('%20.12f %20.12f\n' % (0, u_n))
    u = u_n - 0.5*dt**2*w**2*u_n
    outfile.write('%20.12f %20.12f\n' % (dt, u))
    u_nm1 = u_n
    u_n = u
    for n in range(1, Nt):
        u = 2*u_n - u_nm1 - dt**2*w**2*u_n
        outfile.write('%20.12f %20.12f\n' % (t[n], u))
        u_nm1 = u_n
        u_n = u
    return u, t

def test_solver_memsave():
    from vib_undamped import solver
    _, _ = solver_memsave(I=1, dt=0.1, w=1, T=30)
    u_expected, _ = solver        (I=1, dt=0.1, w=1, T=30)
    data = np.loadtxt('tmp.dat')
    u_computed = data[:,1]
    diff = np.abs(u_expected - u_computed).max()
    assert diff < 5E-13, diff


if __name__ == '__main__':
    test_solver_memsave()
    solver_memsave(I=1, w=1, dt=0.1, T=30)
