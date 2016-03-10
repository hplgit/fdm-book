#!/usr/bin/env python
"""
1D wave equation with u=0 at the boundary.
Simplest possible implementation.

The key function is::

  u, x, t, cpu = solver(I, V, f, c, L, dt, C, T, user_action)

which solves the wave equation u_tt = c**2*u_xx on (0,L) with u=0
on x=0,L, for t in (0,T].  Initial conditions: u=I(x), u_t=V(x).

T is the stop time for the simulation.
dt is the desired time step.
C is the Courant number (=c*dt/dx), which specifies dx.
f(x,t) is a function for the source term (can be 0 or None).
I and V are functions of x.

user_action is a function of (u, x, t, n) where the calling
code can add visualization, error computations, etc.
"""

import numpy as np

def solver(I, V, f, c, L, dt, C, T, user_action=None):
    """Solve u_tt=c^2*u_xx + f on (0,L)x(0,T]."""
    Nt = int(round(T/dt))
    t = np.linspace(0, Nt*dt, Nt+1)   # Mesh points in time
    dx = c*dt/C
    Nx = int(round(L/dx))
    x = np.linspace(0, L, Nx+1)     # Mesh points in space
    C2 = C**2                       # Help variable in the scheme
    # Recompute to make sure dx and dt are compatible with x and t
    dx = x[1] - x[0]
    dt = t[1] - t[0]

    if f is None or f == 0 :
        f = lambda x, t: 0
    if V is None or V == 0:
        V = lambda x: 0

    u   = np.zeros(Nx+1)   # Solution array at new time level
    u_1 = np.zeros(Nx+1)   # Solution at 1 time level back
    u_2 = np.zeros(Nx+1)   # Solution at 2 time levels back

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
               0.5*C2*(u_1[i-1] - 2*u_1[i] + u_1[i+1]) + \
               0.5*dt**2*f(x[i], t[n])
    u[0] = 0;  u[Nx] = 0

    if user_action is not None:
        user_action(u, x, t, 1)

    # Switch variables before next step
    u_2[:] = u_1;  u_1[:] = u

    for n in range(1, Nt):
        # Update all inner points at time t[n+1]
        for i in range(1, Nx):
            u[i] = - u_2[i] + 2*u_1[i] + \
                     C2*(u_1[i-1] - 2*u_1[i] + u_1[i+1]) + \
                     dt**2*f(x[i], t[n])

        # Insert boundary conditions
        u[0] = 0;  u[Nx] = 0
        if user_action is not None:
            if user_action(u, x, t, n+1):
                break

        # Switch variables before next step
        u_2[:] = u_1;  u_1[:] = u

    cpu_time = time.clock() - t0
    return u, x, t, cpu_time

def viz(
    I, V, f, c, L, dt, C, T,  # PDE parameters
    umin, umax,               # Interval for u in plots
    animate=True,             # Simulation with animation?
    tool='matplotlib',        # 'matplotlib' or 'scitools'
    solver_function=solver,   # Function with numerical algorithm
    ):
    """
    Run solver, store and viz. u at each time level with all C values.
    """
    
    class PlotUst:
        def __init__(self):
            self.all_u = []
            self.all_u_for_all_C = []
            self.x_mesh = []   # need each mesh for final plots
        def __call__(self, u, x, t, n):
            """user_action function for solver."""
            self.all_u.append(u.copy())
            if t[n] == T: # i.e., whole time interv. done for this C
                self.x_mesh.append(x)
                self.all_u_for_all_C.append(self.all_u)
                self.all_u = []     # reset to empty list
                
                if len(self.all_u_for_all_C) == len(C):  # all C done
                    print 'Finished all C. Proceed with plots...'
                    # note: n will here be the last index in t[n]
                    for n_ in range(0, n+1):      # for each tn
                        plt.plot(self.x_mesh[0], 
                                 self.all_u_for_all_C[0][n_],
                                 axis=[0, L, umin, umax],
                                 title='Solutions for all \
                                        C at t=%f' % t[n_])
                        plt.hold('on')
                        
                        for j in range(1, len(C)):
                            # build plot at this tn with each 
                            # sol. from the different C values
                            plt.plot(self.x_mesh[j], 
                                     self.all_u_for_all_C[j][n_],
                                     axis=[0, L, umin, umax])
                        plt.xlabel('x'); plt.ylabel('u')
                        plt.hold('off')
                        plt.show()
                        # Let the init. cond. stay on the screen for
                        # 2 sec, else insert a pause of 0.2 s  
                        # between each plot                        
                        time.sleep(2) if t[n_] == 0 else \
                                                    time.sleep(0.2)                        
                        plt.savefig('tmp_%04d.png' % n_)  # for movie
                

    class PlotMatplotlib:
        def __init__(self):
            self.all_u = []
            self.all_u_for_all_C = []
        def __call__(self, u, x, t, n):
            """user_action function for solver."""
            self.all_u.append(u.copy())
            if t[n] == T: # i.e., whole time interv. done for this C
                self.all_u_for_all_C.append(self.all_u)
                self.all_u = []     # reset to empty list
                
                if len(self.all_u_for_all_C) == len(C):  # all C done
                    print 'Finished all C. Proceed with plots...'
                    plt.ion()
                    # note: n will here be the last index in t[n]
                    for n_ in range(0, n+1):      # for each tn
                        plt.plot(x, self.all_u_for_all_C[0][n_])
                        plt.axis([0, L, umin, umax])
                        plt.hold(True)
                        for j in range(1, len(C)):
                            # build plot at this tn with each 
                            # sol. from the different C values
                            plt.plot(x, self.all_u_for_all_C[j][n_])
                        plt.axis([0, L, umin, umax])
                        plt.xlabel('x'); plt.ylabel('u')
                        plt.title('Solutions for all \
                                   C at t=%f' % t[n_])
                        plt.hold(False)
                        plt.draw()
                        # Let the init. cond. stay on the screen for
                        # 2 sec, else insert a pause of 0.2 s  
                        # between each plot                        
                        time.sleep(2) if t[n_] == 0 else \
                                                    time.sleep(0.2)                        
                        plt.savefig('tmp_%04d.png' % n_)  # for movie

    if tool == 'matplotlib':
        import matplotlib.pyplot as plt
        plot_u = PlotMatplotlib()
    elif tool == 'scitools':
        import scitools.std as plt  # scitools.easyviz interface
        plot_u = PlotUst()
    import time, glob, os

    # Clean up old movie frames
    for filename in glob.glob('tmp_*.png'):
        os.remove(filename)

    # Call solver and do the simulaton
    user_action = plot_u if animate else None
    for C_value in C:
        print 'C_value --------------------------------- ', C_value
        u, x, t, cpu = solver_function(
                     I, V, f, c, L, dt, C_value, T, user_action)

    # Make video files
    fps = 4  # frames per second
    codec2ext = dict(flv='flv', libx264='mp4', libvpx='webm',
                     libtheora='ogg')  # video formats
    filespec = 'tmp_%04d.png'
    movie_program = 'ffmpeg'  # or 'avconv'
    for codec in codec2ext:
        ext = codec2ext[codec]
        cmd = '%(movie_program)s -r %(fps)d -i %(filespec)s '\
              '-vcodec %(codec)s movie.%(ext)s' % vars()
        os.system(cmd)

    if tool == 'scitools':
        # Make an HTML play for showing the animation in a browser
        plt.movie('tmp_*.png', encoder='html', fps=fps,
                  output_file='movie.html')
    return cpu

def guitar(C):
    """Triangular wave (pulled guitar string)."""
    L = 0.75
    x0 = 0.8*L
    a = 0.005
    freq = 440
    wavelength = 2*L
    c = freq*wavelength
    from math import pi
    w = 2*pi*freq
    num_periods = 1
    T = 2*pi/w*num_periods
    # Choose dt the same as the stability limit for Nx=50
    dt = L/50./c
    dx = dt*c/float(C)
    # Now dt is considered fixed and a list of C 
    # values is made by reducing increasing the dx value 
    # in steps of 10%. 
    all_C = [C]
    all_C.append(c*dt/(1.1*dx))
    all_C.append(c*dt/(1.2*dx))

    def I(x):
        return a*x/x0 if x < x0 else a/(L-x0)*(L-x)

    umin = -1.2*a;  umax = -umin
    cpu = viz(I, 0, 0, c, L, dt, all_C, T, umin, umax,
                 animate=True, tool='scitools')
    #cpu = viz(I, 0, 0, c, L, dt, all_C, T, umin, umax,
    #             animate=True, tool='matplotlib')
    print 'cpu = ', cpu

if __name__ == '__main__':
    import sys
    try:
        C = float(sys.argv[1])
        print 'C=%g' % C
    except IndexError:
        C = 0.85
    print 'Courant number: %.2f' % C
    # The list of C values will be generated from this C value
    guitar(C)