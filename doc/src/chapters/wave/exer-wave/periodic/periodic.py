#!/usr/bin/env python
"""
1D wave equation with u=0 at the boundary.
The solver function here offers scalar and vectorized versions.
See wave1D_u0_s.py for documentation. The only difference
is that function solver takes an additional argument "version":
version='scalar' implies explicit loops over mesh point,
while version='vectorized' provides a vectorized version.
"""
from numpy import *

def solver(I, V, f, c, L, dt, C, T, user_action=None,
           version='vectorized', bc_left='dudx=0'):
    """Solve u_tt=c^2*u_xx + f on (0,L)x(0,T]."""
    Nt = int(round(T/dt))
    t = linspace(0, Nt*dt, Nt+1)   # Mesh points in time
    dx = dt*c/float(C)
    Nx = int(round(L/dx))
    x = linspace(0, L, Nx+1)       # Mesh points in space
    dx = x[1] - x[0]
    C2 = C**2                      # Help variable in the scheme
    if f is None or f == 0:
        f = (lambda x, t: 0) if version == 'scalar' else \
            lambda x, t: zeros(x.shape)
    if V is None or V == 0:
        V = (lambda x: 0) if version == 'scalar' else \
            lambda x: zeros(x.shape)

    u   = zeros(Nx+1)   # Solution array at new time level
    u_1 = zeros(Nx+1)   # Solution at 1 time level back
    u_2 = zeros(Nx+1)   # Solution at 2 time levels back

    import time;  t0 = time.clock()  # for measuring CPU time

    global periodic
    periodic = False  # True: periodic condition at x=0, False: du/dn=0

    # Load initial condition into u_1
    for i in range(0,Nx+1):
        u_1[i] = I(x[i])

    if user_action is not None:
        user_action(u_1, x, t, 0)

    # Special formula for first time step
    n = 0
    for i in range(1, Nx):
        u[i] = u_1[i] + dt*V(x[i]) + \
               0.5*C2*(u_1[i-1] - 2*u_1[i] + u_1[i+1]) + \
               0.5*dt**2*f(x[i], t[n])
    i = 0
    if bc_left == 'dudx=0':
        u[i] = u_1[i] + dt*V(x[i]) + \
               0.5*C2*(u_1[i+1] - 2*u_1[i] + u_1[i+1]) + \
               0.5*dt**2*f(x[i], t[n])
    else: # open boundary condition
        u[i] = u_1[i] + C*(u_1[i+1] - u_1[i])

    u[Nx] = 0

    if user_action is not None:
        user_action(u, x, t, 1)

    # Switch variables before next step
    u_2[:], u_1[:] = u_1, u

    for n in range(1, Nt):
        # Update all inner points at time t[n+1]

        # Turn to periodic conditions when initial disturbance
        # at x=0 has died out
        if not periodic and u[Nx] > 0.00001:
            periodic = True

        if version == 'scalar':
            for i in range(1, Nx):
                u[i] = - u_2[i] + 2*u_1[i] + \
                       C2*(u_1[i-1] - 2*u_1[i] + u_1[i+1]) + \
                       dt**2*f(x[i], t[n])
        elif version == 'vectorized':   # (1:-1 slice style)
            f_a = f(x, t[n])  # Precompute in array
            u[1:-1] = - u_2[1:-1] + 2*u_1[1:-1] + \
                C2*(u_1[0:-2] - 2*u_1[1:-1] + u_1[2:]) + \
                dt**2*f_a[1:-1]

        # Insert boundary conditions
        u[Nx] = u_1[Nx] - C*(u_1[Nx] - u_1[Nx-1])  # open condition
        if periodic:
            u[0] = u[Nx]
        else:
            i = 0
            if bc_left == 'dudx=0':
                u[i] = - u_2[i] + 2*u_1[i] + \
                       C2*(u_1[i+1] - 2*u_1[i] + u_1[i+1]) + \
                       dt**2*f(x[i], t[n])
            else: # open boundary condition
                u[i] = u_1[i] + C*(u_1[i+1] - u_1[i])

        if user_action is not None:
            if user_action(u, x, t, n+1):
                break

        # Switch variables before next step
        u_2[:], u_1[:] = u_1, u

    cpu_time = t0 - time.clock()
    return u, x, t, cpu_time

def viz(I, V, f, c, L, dt, C, T, umin, umax, animate=True,
        version='vectorized', bc_left='dudx=0'):
    """Run solver and visualize u at each time level."""
    import scitools.std as plt, time, glob, os
    #num_frames = 100 # max no of frames in movie

    def plot_u(u, x, t, n):
        """user_action function for solver."""
        bc0 = 'periodic BC' if periodic else bc_left
        try:
            every = t.size/num_frames
        except NameError:
            every = 1  # plot every frame
        if n % every == 0:
            plt.plot(x, u, 'r-',
                     xlabel='x', ylabel='u',
                     axis=[0, L, umin, umax],
                     title='t=%.3f, x=0: %s, x=L: open BC' % (t[n], bc0))
            # Let the initial condition stay on the screen for 2
            # seconds, else insert a pause of 0.2 s between each plot
            time.sleep(2) if t[n] == 0 else time.sleep(0.2)
            plt.savefig('frame_%04d.png' % n)  # for movie making

    # Clean up old movie frames
    for filename in glob.glob('frame_*.png'):
        os.remove(filename)

    user_action = plot_u if animate else None
    u, x, t, cpu = solver(I, V, f, c, L, dt, C, T,
                          user_action, version, bc_left)
    if not animate:
        return cpu

    # Make movie files
    fps = 6  # Frames per second
    plt.movie('frame_*.png', encoder='html', fps=fps,
              output_file='movie.html')
    # Ex: avconv -r 4 -i frame_%04d.png -vcodec libtheora movie.ogg
    codec2ext = dict(flv='flv', libx64='mp4', libvpx='webm',
                     libtheora='ogg')
    filespec = 'frame_%04d.png'
    movie_program = 'avconv'  # or 'ffmpeg'
    for codec in codec2ext:
        ext = codec2ext[codec]
        cmd = '%(movie_program)s -r %(fps)d -i %(filespec)s '\
              '-vcodec %(codec)s movie.%(ext)s' % vars()
        print cmd
        os.system(cmd)
    return cpu


def plug(C=1, Nx=50, animate=True, T=2, loc=0):
    """Plug profile as initial condition."""
    L = 1.
    c = 1

    I = lambda x: 1 if abs(x-loc) < 0.1 else 0

    bc_left = 'dudx=0' if loc == 0 else 'open'
    dt = (L/Nx)/c  # choose the stability limit with given Nx
    cpu = viz(I, None, None, c, L, dt, C, T,
              umin=-0.3, umax=1.1, animate=animate,
              bc_left=bc_left)

def gaussian(C=1, Nx=50, animate=True, T=2, loc=0):
    """Gaussian bell as initial condition."""
    L = 1.
    c = 1

    def I(x):
        return exp(-0.5*((x-loc)/0.05)**2)

    bc_left = 'dudx=0' if loc == 0 else 'open'
    dt = (L/Nx)/c  # choose the stability limit with given Nx
    cpu = viz(I, None, None, c, L, dt, C, T,
              umin=-0.2, umax=1, animate=animate, bc_left=bc_left)

if __name__ == '__main__':
    import sys
    from scitools.misc import function_UI
    cmd = function_UI([plug, gaussian], sys.argv)
    eval(cmd)
