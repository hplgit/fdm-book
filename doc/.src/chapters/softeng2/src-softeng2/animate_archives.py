"""
Given some archives of .npz files with t,x,u0001,u0002,... data,
make a simultaneous animation of the data in the archives.
Ideal for comparing different simulations.
"""
import glob, os, sys, time
import numpy as np

def animate_multiple_solutions(
    archives, umin, umax, pause=0.2, show=True):
    """
    Animate data in list "archives" of numpy.savez archive files.
    Each archive in archives holds t, x, and u0000, u0001, u0002,
    and so on, where t and x are arrays for the time and space
    meshes, resp., and the u* arrays are solutions u(x) at the
    various time points.
    umin and umax are overall min and max values of the data.
    A pause is inserted between each frame in the screen animation.
    Each frame is stored in a file tmp_%04d.png (for video making).
    No screen animation takes place if show=False.

    The animation applies the time points in the coarsest time
    mesh found in the archives.
    Linear interpolation is used to calculate data at these times
    for the solutions in the other archives.
    """
    import matplotlib.pyplot as plt

    simulations = [np.load(archive) for archive in archives]
    # simulations is list of "dicts", each dict is {'t': array,
    # 'x': array, 'u0': array, 'u1': array, ...}
    # Must base animation on coarsest resolution
    coarsest = np.argmin([len(s['t']) for s in simulations])

    # Create plot with all solutions for very first timestep

    n = 0
    sol = 'u%04d' % n
    # Build plot command
    plot_all_sols_at_tn = 'lines = plt.plot(' + ', '.join(
        ["simulations[%d]['x'], simulations[%d]['%s']" %
         (sim, sim, sol) for sim in range(0, len(simulations))]) \
         + ')'

    if show:
        plt.ion()
    else:
        plt.ioff()

    exec plot_all_sols_at_tn    # run plot command
    plt.savefig('tmp_%04d.png' % n)

    plt.xlabel('x'); plt.ylabel('u')
    # size of domain is the same for all sims
    plt.axis([simulations[0]['x'][0],
              simulations[0]['x'][-1],
              umin, umax])
    plt.legend(['t=%.3f' % simulations[0]['t'][n]])
    # How suppress drawing on the screen?
    plt.draw()
    plt.savefig('tmp_%04d.png' % n)
    time.sleep(1)

    # Find legends
    t = simulations[coarsest]['t']
    x = simulations[coarsest]['x']
    dt = t[1] - t[0]
    dx = x[1] - x[0]
    legends = [r'$\Delta x=%.4f, \Delta t=%.4f$' % (dx, dt)]
    for i in range(len(simulations)):
        if i != coarsest:
            t = simulations[i]['t']
            x = simulations[i]['x']
            dt = t[1] - t[0]
            dx = x[1] - x[0]
            legends.append(r'$\Delta x=%.4f, \Delta t=%.4f$' %
                           (dx, dt))

    # Plot all solutions at each remaining time step

    # At every time step set_ydata has to be executed
    # with the new solutions from all simulations.
    # (note that xdata remains unchanged in time)
    interpolation_details = []  # for testing
    t = simulations[coarsest]['t']
    for n in range(1, len(t)):
        sol = 'u%04d' % (n)
        lines[coarsest].set_ydata(simulations[coarsest][sol])
        interpolation_details.append([])
        for k in range(len(simulations)):
            if k != coarsest:
                # Interpolate simulations at t[n] (in coarsest
                # simulation) among all simulations[k]
                a, i, w = interpolate_arrays(
                    simulations[k], simulations[k]['t'], t[n])
                lines[k].set_ydata(a)
                interpolation_details[-1].append([t[n], i, w])
            else:
                interpolation_details[-1].append('coarsest')

        plt.legend(legends)
        plt.title('t=%.3f' % (simulations[0]['t'][n]))
        plt.draw()
        plt.savefig('tmp_%04d.png' % n)
        time.sleep(pause)
    return interpolation_details


def linear_interpolation(t, tp):
    """
    Given an array of time values, with constant spacing,
    and some time point tp, determine the data for linear
    interpolation: i and w such that
    tp = (1-w)*t[i] + w*t[i+1]. If tp happens to equal t[i]
    for any i, return i and None.
    """
    # Determine time cell
    dt = float(t[1] - t[0])  # assumed constant!
    i = int(tp/dt)
    if abs(tp - t[i]) < 1E-13:
        return i, None
    #tp = t[i] + w*dt
    w = (tp - t[i])/dt
    return i, w

def interpolate_arrays(arrays, t, tp):
    """
    Given a time point tp and a collection of arrays corresponding
    to times in array t, perform linear interpolation among array
    i and i+1 when t[i] < tp < t[i+1].
    arrays can be .npz archive (NpzFile) or list of numpy arrays.
    Return interpolated array, i, w (where w is the interpolation
    weight found in linear_interpolation).
    """
    i, w = linear_interpolation(t, tp)

    if isinstance(arrays, np.lib.npyio.NpzFile):
        # arrays behaves as a dict with keys u0001, u0002, ...
        if w is None:
            return arrays['u%04d' % i], i, None
        else:
            return w*arrays['u%04d' % i] + \
                   (1-w)*arrays['u%04d' % (i+1)], i, w

    elif isinstance(arrays, (tuple,list)) and \
         isinstance(arrays[0], np.ndarray):
        # arrays is list/tuple of arrays
        if w is None:
            return arrays[i], i, None
        else:
            return w*arrays[i] + (1-w)*arrays[i+1], i, w
    else:
        raise TypeError(
            'arrays is %s, must be NpzFile archive or list arrays'
            % type(arrays))


def demo_animate_multiple_solutions():
    '''First run all simulations. Then animate all from archives.'''
    # Must delete all archives so we really recompute them
    # and get their names from the pulse function
    for filename in glob.glob('.*.npz') + glob.glob('tmp_*.png'):
        os.remove(filename)
    archives = []
    umin = umax = 0
    from wave1D_dn_vc import pulse
    for spatial_resolution in [20, 55, 200]:
        archive_name, u_min, u_max = pulse(
            Nx=spatial_resolution, pulse_tp='gaussian')
        archives.append(archive_name)
        if u_min < umin: umin = u_min
        if u_max > umax: umax = u_max

    print archives
    animate_multiple_solutions(archives, umin, umax, show=True)
    cmd = 'ffmpeg -i tmp_%04d.png -r 25 -vcodec libtheora movie.ogg'
    os.system(cmd)

def test_animate_multiple_solutions():
    # Must delete all archives so we really recompute them
    # and get their names from the pulse function
    for filename in glob.glob('.*.npz') + glob.glob('tmp_*.png'):
        os.remove(filename)
    archives = []
    umin = umax = 0
    from wave1D_dn_vc import pulse
    for spatial_resolution in [20, 45, 100]:
        archive_name, u_min, u_max = pulse(
            Nx=spatial_resolution, pulse_tp='gaussian')
        archives.append(archive_name)
        if u_min < umin: umin = u_min
        if u_max > umax: umax = u_max

    print archives
    details = animate_multiple_solutions(
        archives, umin, umax, show=False)
    # Round data:
    for i in range(len(details)):
        for j in range(len(details[i])):
            if details[i][j] == 'coarsest':
                continue
            details[i][j][0] = round(details[i][j][0], 4)
            if isinstance(details[i][j][2], float):
                details[i][j][2] = round(details[i][j][2], 4)
    expected = [
        ['coarsest', [0.05, 2, 0.25], [0.05, 5, None]],
        ['coarsest', [0.1, 4, 0.5], [0.1, 10, None]],
        ['coarsest', [0.15, 6, 0.75], [0.15, 15, None]],
        ['coarsest', [0.2, 9, None], [0.2, 20, None]],
        ['coarsest', [0.25, 11, 0.25], [0.25, 25, None]],
        ['coarsest', [0.3, 13, 0.5], [0.3, 30, None]],
        ['coarsest', [0.35, 15, 0.75], [0.35, 35, None]],
        ['coarsest', [0.4, 18, None], [0.4, 40, None]],
        ['coarsest', [0.45, 20, 0.25], [0.45, 45, None]],
        ['coarsest', [0.5, 22, 0.5], [0.5, 50, None]],
        ['coarsest', [0.55, 24, 0.75], [0.55, 55, None]],
        ['coarsest', [0.6, 27, None], [0.6, 60, None]],
        ['coarsest', [0.65, 29, 0.25], [0.65, 65, None]],
        ['coarsest', [0.7, 31, 0.5], [0.7, 70, None]],
        ['coarsest', [0.75, 33, 0.75], [0.75, 75, None]],
        ['coarsest', [0.8, 36, None], [0.8, 80, None]],
        ['coarsest', [0.85, 38, 0.25], [0.85, 85, None]],
        ['coarsest', [0.9, 40, 0.5], [0.9, 90, None]],
        ['coarsest', [0.95, 42, 0.75], [0.95, 95, None]],
        ['coarsest', [1.0, 45, None], [1.0, 100, None]],
        ['coarsest', [1.05, 47, 0.25], [1.05, 105, None]],
        ['coarsest', [1.1, 49, 0.5], [1.1, 110, None]],
        ['coarsest', [1.15, 51, 0.75], [1.15, 115, None]],
        ['coarsest', [1.2, 54, None], [1.2, 120, None]],
        ['coarsest', [1.25, 56, 0.25], [1.25, 125, None]],
        ['coarsest', [1.3, 58, 0.5], [1.3, 130, None]],
        ['coarsest', [1.35, 60, 0.75], [1.35, 135, None]],
        ['coarsest', [1.4, 63, None], [1.4, 140, None]],
        ['coarsest', [1.45, 65, 0.25], [1.45, 145, None]],
        ['coarsest', [1.5, 67, 0.5], [1.5, 150, None]],
        ['coarsest', [1.55, 69, 0.75], [1.55, 155, None]],
        ['coarsest', [1.6, 72, None], [1.6, 160, None]],
        ['coarsest', [1.65, 74, 0.25], [1.65, 165, None]],
        ['coarsest', [1.7, 76, 0.5], [1.7, 170, None]],
        ['coarsest', [1.75, 78, 0.75], [1.75, 175, None]],
        ['coarsest', [1.8, 81, None], [1.8, 180, None]],
        ['coarsest', [1.85, 83, 0.25], [1.85, 185, None]],
        ['coarsest', [1.9, 85, 0.5], [1.9, 190, None]],
        ['coarsest', [1.95, 87, 0.75], [1.95, 195, None]],
        ['coarsest', [2.0, 90, None], [2.0, 200, None]],
        ]
    assert details == expected

if __name__ == '__main__':
    #test_animate_multiple_solutions()
    umin = float(sys.argv[1])
    umax = float(sys.argv[2])
    archives = sys.argv[3:]
    animate_multiple_solutions(
        archives, umin, umax, pause=0.2, show=True)
