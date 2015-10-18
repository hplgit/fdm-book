#!/usr/bin/env python
"""
Modification of wave1D_u0_sv.py for simulating waves on a
string with varying density.
"""
from numpy import *

# Change from parameter c to density and tension
# density is a two-list, tension is a constant
# We assume the jump in density is at x=L/2, but
# this could be a parameter.
# Note that the C2 help variable becomes different
# from the original code (C is misleading here since
# it actually has two values through the mesh, it is
# better to just use dt, dx, tension, and rho in
# the scheme!).

def solver(I, V, f, density, tension, L, Nx, C, T, user_action=None,
           version='vectorized'):
    """Solve u_tt=c^2*u_xx + f on (0,L)x(0,T]."""
    x = linspace(0, L, Nx+1)     # Mesh points in space
    dx = x[1] - x[0]
    # Must use largest wave velocity (c=sqrt(tension/density))
    # in stability criterion
    c = sqrt(tension/min(density))
    dt = C*dx/c
    # This means that the given C is applied to the medium with smallest
    # density, while a lower effective C=c*dt/dx appears in the medium
    # with the highest density.
    rho = zeros(Nx+1)
    rho[:len(rho)/2] = density[0]
    rho[len(rho)/2:] = density[1]

    Nt = int(round(T/dt))
    t = linspace(0, Nt*dt, Nt+1) # Mesh points in time
    C2 = tension*dt**2/dx**2     # Help variable in the scheme

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

    # Load initial condition into u_1
    for i in range(0,Nx+1):
        u_1[i] = I(x[i])

    if user_action is not None:
        user_action(u_1, x, t, 0)

    # Special formula for first time step
    n = 0
    for i in range(1, Nx):
        u[i] = u_1[i] + dt*V(x[i]) + \
               1/rho[i]*(
               0.5*C2*(u_1[i-1] - 2*u_1[i] + u_1[i+1]) + \
               0.5*dt**2*f(x[i], t[n]))
    u[0] = 0;  u[Nx] = 0

    if user_action is not None:
        user_action(u, x, t, 1)

    # Switch variables before next step
    u_2[:], u_1[:] = u_1, u

    for n in range(1, Nt):
        # Update all inner points at time t[n+1]

        if version == 'scalar':
            for i in range(1, Nx):
                u[i] = - u_2[i] + 2*u_1[i] + \
                       1/rho[i]*(
                       C2*(u_1[i-1] - 2*u_1[i] + u_1[i+1]) + \
                       dt**2*f(x[i], t[n]))
        elif version == 'vectorized':   # (1:-1 slice style)
            f_a = f(x, t[n])  # Precompute in array
            u[1:-1] = - u_2[1:-1] + 2*u_1[1:-1] + \
                1/rho[1:-1]*(
                C2*(u_1[0:-2] - 2*u_1[1:-1] + u_1[2:]) + \
                dt**2*f_a[1:-1])

        # Insert boundary conditions
        u[0] = 0;  u[Nx] = 0
        if user_action is not None:
            if user_action(u, x, t, n+1):
                break

        # Switch variables before next step
        u_2[:], u_1[:] = u_1, u

    cpu_time = t0 - time.clock()
    return u, x, t, cpu_time

def viz(I, V, f, density, tension, L, Nx, C, T, umin, umax,
        animate=True,
        movie_filename='movie',
        version='vectorized'):
    """Run solver and visualize u at each time level."""
    import scitools.std as plt, time, glob, os
    #num_frames = 100 # max no of frames in movie

    def plot_u(u, x, t, n):
        """user_action function for solver."""
        try:
            every = t.size/num_frames
        except NameError:
            every = 1  # plot every frame
        if n % every == 0:
            plt.plot(x, u, 'r-',
                     xlabel='x', ylabel='u',
                     axis=[0, L, umin, umax],
                     title='t=%f' % t[n])
            # Let the initial condition stay on the screen for 2
            # seconds, else insert a pause of 0.2 s between each plot
            time.sleep(2) if t[n] == 0 else time.sleep(0.2)
            plt.savefig('frame_%04d.png' % n)  # for movie making

    # Clean up old movie frames
    for filename in glob.glob('frame_*.png'):
        os.remove(filename)

    user_action = plot_u if animate else None
    u, x, t, cpu = solver(I, V, f, density, tension, L, Nx, C, T,
                          user_action, version)
    if not animate:
        return cpu

    # Make movie files
    fps = 4  # Frames per second
    plt.movie('frame_*.png', encoder='html', fps=fps,
              output_file='movie.html')
    # Ex: avconv -r 4 -i frame_%04d.png -vcodec libtheora movie.ogg
    #codec2ext = dict(flv='flv', libx64='mp4', libvpx='webm',
    #                 libtheora='ogg')
    codec2ext = dict(libtheora='ogg')
    filespec = 'frame_%04d.png'
    movie_program = 'avconv'  # or 'ffmpeg'
    for codec in codec2ext:
        ext = codec2ext[codec]
        cmd = '%(movie_program)s -r %(fps)d -i %(filespec)s '\
              '-vcodec %(codec)s %(movie_filename)s.%(ext)s' % vars()
        os.system(cmd)
    return cpu

import nose.tools as nt

def test_quadratic():
    """
    Check the scalar and vectorized versions work for
    a quadratic u(x,t)=x(L-x)(1+t/2) that is exactly reproduced.
    """
    # The following function must work for x as array or scalar
    exact_solution = lambda x, t: x*(L - x)*(1 + 0.5*t)
    I = lambda x: exact_solution(x, 0)
    V = lambda x: 0.5*exact_solution(x, 0)
    # f is a scalar (zeros_like(x) works for scalar x too)
    f = lambda x, t: zeros_like(x) + 2*c**2*(1 + 0.5*t)

    L = 2.5
    c = 1.5
    Nx = 3  # Very coarse mesh
    C = 1
    T = 18  # Long time integration

    tension = 1 # just some number
    # density follows from c=sqrt(tension/density)
    density = [tension/c**2, tension/c**2]

    def assert_no_error(u, x, t, n):
        u_e = exact_solution(x, t[n])
        diff = abs(u - u_e).max()
        nt.assert_almost_equal(diff, 0, places=13)

    #solver(I, V, f, density, tension, L, Nx, C, T,
    #       user_action=assert_no_error, version='scalar')
    solver(I, V, f, density, tension, L, Nx, C, T,
           user_action=assert_no_error, version='vectorized')

def guitar():
    """Triangular wave (pulled guitar string)."""
    L = 0.75
    x0 = 0.8*L
    a = 0.005
    freq = 440
    wavelength = 2*L

    def I(x):
        return a*x/x0 if x < x0 else a/(L-x0)*(L-x)

    umin = -1.2*a;  umax = -umin

    # Relevant c from frequency (A tone) and wave length (string length)
    c = freq*wavelength
    # c = sqrt(tension/density)
    # Set tension to 150 Newton (nylon string)
    tension = 150

    for jump in [0.1, 10]:
        # Jump to the right where C=1, while 1/jump is the effective Courant
        # number to the left
        density = array([tension/c**2, jump*tension/c**2])
        # Compute period (in time) for the two pieces
        # (not sure the reasoning for computing T is correct, the
        # largest jump seems to give a shorter time history of the string)
        c_variable = sqrt(tension/density)
        omega = 2*pi*c_variable/wavelength  # omega = 2*pi*freq
        P = 2*pi/omega.min()                # longest period of the two
        num_periods = 1
        T = P*num_periods
        Nx = 50
        C = 1.0  # perfect wave for smallest density

        print '*** Simulating with jump=%g' % jump
        cpu = viz(I, 0, 0, density, tension, L, Nx, C, T,
                  umin, umax, animate=True,
                  movie_filename='wave1D_u0_sv_discont_jump%g' % jump)


if __name__ == '__main__':
    #test_quadratic()  # verify
    guitar()

