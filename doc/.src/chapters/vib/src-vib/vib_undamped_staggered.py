import sys, os
sys.path.insert(0, os.path.join(os.pardir, os.pardir, 'vib', 'src-vib'))
from numpy import zeros, linspace

def solver_v1(I, w, dt, T):
    """
    Solve u'=v, v' = - w**2*u for t in (0,T], u(0)=I and v(0)=0,
    by a central finite difference method with time step dt.
    """
    dt = float(dt)
    Nt = int(round(T/dt))
    u = zeros(Nt+1)
    v = zeros(Nt+1)
    t = linspace(0, Nt*dt, Nt+1)  # mesh for u
    t_v = (t + dt/2)[:-1]         # mesh for v

    u[0] = I
    v[0] = 0 - 0.5*dt*w**2*u[0]
    for n in range(1, Nt+1):
        u[n] = u[n-1] + dt*v[n-1]
        v[n] = v[n-1] - dt*w**2*u[n]
    return u, t, v, t_v

class HalfInt:
    """
    Class for allowing to write n+half and mean n,
    while n-half is n-1. Used for nice notation in staggered
    meshes.
    """
    def __radd__(self, other):
        return other

    def __rsub__(self, other):
        return other - 1

half = HalfInt()  # singleton object

def solver(I, w, dt, T):
    """
    Solve u'=v, v' = - w**2*u for t in (0,T], u(0)=I and v(0)=0,
    by a central finite difference method with time step dt on
    a staggered mesh with v as unknown at (i+1/2)*dt time points.
    """
    dt = float(dt)
    Nt = int(round(T/dt))
    u = zeros(Nt+1)
    v = zeros(Nt+1)
    t = linspace(0, Nt*dt, Nt+1)  # mesh for u
    t_v = t + dt/2                # mesh for v

    u[0] = I
    v[0+half] = 0 - 0.5*dt*w**2*u[0]
    for n in range(1, Nt+1):
        u[n] = u[n-1] + dt*v[n-half]
        v[n+half] = v[n-half] - dt*w**2*u[n]
    return u, t, v[:-1], t_v[:-1]

def test_staggered():
    I = 1.2; w = 2.0; T = 5; dt = 2/w
    u, t, v, t_v = solver(I, w, dt, T)
    from vib_undamped import solver as solver2
    u2, t2 = solver2(I, w, dt, T)
    error = abs(u - u2).max()
    tol = 1E-14
    assert error < tol

def test_convergence():
    """Verify 2nd-order convergence."""
    from vib_undamped import convergence_rates
    def wrapped_solver(I, w, dt, T):
        u, t, v, t_v = solver(I, w, dt, T)
        return u, t

    r = convergence_rates(8, wrapped_solver, 8)
    print r
    assert abs(r[-1] - 2) < 1E-5

if __name__ == '__main__':
    test_staggered()
    test_convergence()
