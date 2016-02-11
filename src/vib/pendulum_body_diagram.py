"""
Animate a body diagram for the motion of a pendulum.
The visualization is coupled to Pysketcher.
"""
import sys
try:
    from pysketcher import *
except ImportError:
    print 'Pysketcher must be installed from'
    print 'https://github.com/hplgit/pysketcher'
    sys.exit(1)

# Overall dimensions of sketch
H = 15.
W = 17.

drawing_tool.set_coordinate_system(
    xmin=0, xmax=W, ymin=0, ymax=H,
    axis=False)

def sketch(theta, S, mg, drag, t, time_level):
    """
    Draw pendulum sketch with body forces at a time level
    corresponding to time t. The drag force is in
    drag[time_level], the force in the wire is S[time_level],
    the angle is theta[time_level].
    """
    import math
    a = math.degrees(theta[time_level])  # angle in degrees
    L = 0.4*H         # Length of pendulum
    P = (W/2, 0.8*H)  # Fixed rotation point

    mass_pt = path.geometric_features()['end']
    rod = Line(P, mass_pt)

    mass = Circle(center=mass_pt, radius=L/20.)
    mass.set_filled_curves(color='blue')
    rod_vec = rod.geometric_features()['end'] - \
              rod.geometric_features()['start']
    unit_rod_vec = unit_vec(rod_vec)
    mass_symbol = Text('$m$', mass_pt + L/10*unit_rod_vec)

    rod_start = rod.geometric_features()['start']  # Point P
    vertical = Line(rod_start, rod_start + point(0,-L/3))

    def set_dashed_thin_blackline(*objects):
        """Set linestyle of objects to dashed, black, width=1."""
        for obj in objects:
            obj.set_linestyle('dashed')
            obj.set_linecolor('black')
            obj.set_linewidth(1)

    set_dashed_thin_blackline(vertical)
    set_dashed_thin_blackline(rod)
    angle = Arc_wText(r'$\theta$', rod_start, L/6, -90, a,
                      text_spacing=1/30.)

    magnitude = 1.2*L/2   # length of a unit force in figure
    force = mg[time_level]  # constant (scaled eq: about 1)
    force *= magnitude
    mg_force  = Force(mass_pt, mass_pt + force*point(0,-1),
                      '', text_pos='end')
    force = S[time_level]
    force *= magnitude
    rod_force = Force(mass_pt, mass_pt - force*unit_vec(rod_vec),
                      '', text_pos='end',
                      text_spacing=(0.03, 0.01))
    force = drag[time_level]
    force *= magnitude
    air_force = Force(mass_pt, mass_pt -
                      force*unit_vec((rod_vec[1], -rod_vec[0])),
                      '', text_pos='end',
                      text_spacing=(0.04,0.005))

    body_diagram = Composition(
        {'mg': mg_force, 'S': rod_force, 'air': air_force,
         'rod': rod, 'body': mass
         'vertical': vertical, 'theta': angle,})

    body_diagram.draw(verbose=0)
    drawing_tool.savefig('tmp_%04d.png' % time_level, crop=False)
    # (No cropping: otherwise movies will be very strange!)

def simulate(alpha, Theta, dt, T):
    import odespy

    def f(u, t, alpha):
        omega, theta = u
        return [-alpha*omega*abs(omega) - sin(theta),
                omega]

    import numpy as np
    Nt = int(round(T/float(dt)))
    t = np.linspace(0, Nt*dt, Nt+1)
    solver = odespy.RK4(f, f_args=[alpha])
    solver.set_initial_condition([0, Theta])
    u, t = solver.solve(
        t, terminate=lambda u, t, n: abs(u[n,1]) < 1E-3)
    omega = u[:,0]
    theta = u[:,1]
    S = omega**2 + np.cos(theta)
    drag = -alpha*np.abs(omega)*omega
    return t, theta, omega, S, drag

def animate():
    # Clean up old plot files
    import os, glob
    for filename in glob.glob('tmp_*.png') + glob.glob('movie.*'):
        os.remove(filename)
    # Solve problem
    from math import pi, radians, degrees
    import numpy as np
    alpha = 0.4
    period = 2*pi   # Use small theta approximation
    T = 12*period   # Simulate for 12 periods
    dt = period/40  # 40 time steps per period
    a = 70          # Initial amplitude in degrees
    Theta = radians(a)

    t, theta, omega, S, drag = simulate(alpha, Theta, dt, T)

    # Visualize drag force 5 times as large
    drag *= 5
    mg = np.ones(S.size)  # Gravity force (needed in sketch)

    # Draw animation
    import time
    for time_level, t_ in enumerate(t):
        sketch(theta, S, mg, drag, t_, time_level)
        time.sleep(0.2)  # Pause between each frame on the screen

    # Make videos
    prog = 'ffmpeg'
    filename = 'tmp_%04d.png'
    fps = 6
    codecs = {'flv': 'flv', 'mp4': 'libx264',
              'webm': 'libvpx', 'ogg': 'libtheora'}
    for ext in codecs:
        lib = codecs[ext]
        cmd = '%(prog)s -i %(filename)s -r %(fps)s ' % vars()
        cmd += '-vcodec %(lib)s movie.%(ext)s' % vars()
        print(cmd)
        os.system(cmd)

if __name__ == '__main__':
    animate()
