#!/usr/bin/env python
"""
1D wave equation with homogeneous Neumann conditions::

  u, x, t, cpu = solver(I, V, f, c, L, dt, C, T, user_action)

Function solver solves the wave equation

   u_tt = c**2*u_xx + f(x,t) on

(0,L) with du/dn=0 on x=0 and x = L.

dt is the desired time step.
T is the stop time for the simulation.
C is the Courant number (=c*dt/dx).
dx is computed on basis of dt and C.

I and f are functions: I(x), f(x,t).
user_action is a function of (u, x, t, n) where the calling code
can add visualization, error computations, data analysis,
store solutions, etc.

Function viz::

  viz(I, V, f, c, L, dt, C, T, umin, umax, animate=True)

calls solver with a user_action function that can plot the
solution on the screen (as an animation).
"""
import numpy as np

def solver(I, V, f, c, L, dt, C, T, user_action=None):
    """
    Solve u_tt=c^2*u_xx + f on (0,L)x(0,T].
    u(0,t)=U_0(t) or du/dn=0 (U_0=None), u(L,t)=U_L(t) or du/dn=0 (u_L=None).
    """
    Nt = int(round(T/dt))
    t = np.linspace(0, Nt*dt, Nt+1)   # Mesh points in time
    dx = dt*c/float(C)
    Nx = int(round(L/dx))
    x = np.linspace(0, L, Nx+1)       # Mesh points in space
    C2 = C**2; dt2 = dt*dt            # Help variables in the scheme
    # Make sure dx and dt are compatible with x and t
    dx = x[1] - x[0]
    dt = t[1] - t[0]

    # Wrap user-given f, V
    if f is None or f == 0:
        f = (lambda x, t: 0)
    if V is None or V == 0:
        V = (lambda x: 0)

    u     = np.zeros(Nx+1)   # Solution array at new time level
    u_n   = np.zeros(Nx+1)   # Solution at 1 time level back
    u_nm1 = np.zeros(Nx+1)   # Solution at 2 time levels back

    import time;  t0 = time.clock()  # CPU time measurement

    # Load initial condition into u_n
    for i in range(0, Nx+1):
        u_n[i] = I(x[i])

    if user_action is not None:
        user_action(u_n, x, t, 0)

    # Special formula for the first step
    for i in range(0, Nx+1):
        ip1 = i+1 if i < Nx else i-1
        im1 = i-1 if i > 0  else i+1
        u[i] = u_n[i] + dt*V(x[i]) + \
               0.5*C2*(u_n[im1] - 2*u_n[i] + u_n[ip1]) + \
               0.5*dt2*f(x[i], t[0])

    if user_action is not None:
        user_action(u, x, t, 1)

    # Update data structures for next step
    #u_nm1[:] = u_n;  u_n[:] = u  # safe, but slower
    u_nm1, u_n, u = u_n, u, u_nm1

    for n in range(1, Nt):
        for i in range(0, Nx+1):
            ip1 = i+1 if i < Nx else i-1
            im1 = i-1 if i > 0  else i+1
            u[i] = - u_nm1[i] + 2*u_n[i] + \
                   C2*(u_n[im1] - 2*u_n[i] + u_n[ip1]) + \
                   dt2*f(x[i], t[n])

        if user_action is not None:
            if user_action(u, x, t, n+1):
                break

        # Update data structures for next step
        #u_nm1[:] = u_n;  u_n[:] = u  # safe, but slower
        u_nm1, u_n, u = u_n, u, u_nm1

    # Wrong assignment u = u_nm1 must be corrected before return
    u = u_n
    cpu_time = time.clock() - t0
    return u, x, t, cpu_time

from wave1D_u0 import viz

def plug(C=1, Nx=50, animate=True, T=2):
    """Plug profile as initial condition."""
    def I(x):
        if abs(x-L/2.0) > 0.1:
            return 0
        else:
            return 1

    L = 1.
    c = 1
    dt = (L/Nx)/c  # choose the stability limit with given Nx
    cpu = viz(I, None, None, c, L, dt, C, T,
              umin=-1.1, umax=1.1, animate=animate)

def test_plug():
    """
    Check that an initial plug is correct back after one period,
    if C=1.
    """
    L = 1.0
    I = lambda x: 0 if abs(x-L/2.0) > 0.1 else 1

    Nx = 10
    c = 0.5
    C = 1
    dt = C*(L/Nx)/c
    nperiods = 4
    T = L/c*nperiods  # One period: c*T = L
    u, x, t, cpu = solver(
        I=I, V=None, f=None, c=c, L=L,
        dt=dt, C=C, T=T, user_action=None)
    u_0 = np.array([I(x_) for x_ in x])
    diff = np.abs(u - u_0).max()
    tol = 1E-13
    assert diff < tol

if __name__ == '__main__':
    test_plug()
