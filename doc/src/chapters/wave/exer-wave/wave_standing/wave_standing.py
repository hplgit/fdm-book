import sys, os
sys.path.insert(0, os.path.join(os.pardir, os.pardir, 'src-wave', 'wave1D'))

from wave1D_u0 import solver
#from wave1D_u0v import solver  # allows faster vectorized operations
import numpy as np

def viz(I, V, f, c, L, Nx, C, T,
        ymax,                      # y axis: [-ymax, ymax]
        u_exact,                   # u_exact(x, t)
        animate='u and u_exact',   # or 'error'
        movie_filename='movie',
        ):
    """Run solver and visualize u at each time level."""
    import scitools.std as plt
    import time, glob, os

    class Plot:
        def __init__(self, ymax, frame_name='frame'):
            self.max_error = []   # hold max amplitude errors
            self.max_error_t = [] # time points corresponding to max_error
            self.frame_name = frame_name
            self.ymax = ymax

        def __call__(self, u, x, t, n):
            """user_action function for solver."""
            if animate == 'u and u_exact':
                plt.plot(x, u, 'r-',
                         x, u_exact(x, t[n]), 'b--',
                         xlabel='x', ylabel='u',
                         axis=[0, L, -self.ymax, self.ymax],
                         title='t=%f' % t[n], show=True)
            else:
                error = u_exact(x, t[n]) - u
                local_max_error = np.abs(error).max()
                # self.max_error holds the increasing amplitude error
                if self.max_error == [] or \
                       local_max_error > max(self.max_error):
                    self.max_error.append(local_max_error)
                    self.max_error_t.append(t[n])
                # Use user's ymax until the error exceeds that value.
                # This gives a growing max value of the yaxis (but
                # never shrinking)
                self.ymax = max(self.ymax, max(self.max_error))
                plt.plot(x, error, 'r-',
                         xlabel='x', ylabel='error',
                         axis=[0, L, -self.ymax, self.ymax],
                         title='t=%f' % t[n], show=True)
            plt.savefig('%s_%04d.png' % (self.frame_name, n))

    # Clean up old movie frames
    for filename in glob.glob('frame_*.png'):
        os.remove(filename)

    plot = Plot(ymax)
    u, x, t, cpu = solver(I, V, f, c, L, Nx, C, T, plot)

    # Make plot of max error versus time
    plt.figure()
    plt.plot(plot.max_error_t, plot.max_error)
    plt.xlabel('time'); plt.ylabel('max abs(error)')
    plt.savefig('error.png')
    plt.savefig('error.pdf')

    # Make .flv movie file
    fps = 4  # Frames per second
    codec2ext = dict(flv='flv') #, libx64='mp4', libvpx='webm', libtheora='ogg')

    filespec = 'frame_%04d.png'
    movie_program = 'avconv'  # or 'ffmpeg'
    for codec in codec2ext:
        ext = codec2ext[codec]
        cmd = '%(movie_program)s -r %(fps)d -i %(filespec)s '\
              '-vcodec %(codec)s %(movie_filename)s.%(ext)s' % vars()
        os.system(cmd)

def simulations():
    from numpy import sin, cos, pi
    L = 12  # length of domain
    m = 8   # 2L/m is the wave length or period in space (2*pi/k, k=pi*m/L)
    c = 2   # wave velocity
    A = 1   # amplitude
    Nx = 80
    C = 0.8
    P = 2*pi/(pi*m*c/L)  # 1 period in time
    T = 6*P

    def u_exact(x, t):
        return A*sin(pi*m*x/L)*cos(pi*m*c*t/L)

    def I(x):
        return u_exact(x, 0)

    V = 0
    f = 0

    viz(I, V, f, c, L, Nx, C, 10.5*P,
        0.1, u_exact,
        animate='error',
        movie_filename='error')

    # Very long simulation to demonstrate different curves
    viz(I, V, f, c, L, Nx, C, 30*P,
        1.2*A, u_exact,
        animate='u and u_exact',
        movie_filename='solution')

simulations()
