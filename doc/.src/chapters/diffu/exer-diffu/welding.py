import sys, os
sys.path.insert(0, os.path.join(os.pardir, 'src-diffu'))
from diffu1D_vc import solver
import numpy as np

def run(gamma, beta=10, delta=40, scaling=1, animate=False):
    """Run the scaled model for welding."""
    if scaling == 1:
        v = gamma
        a = 1
    elif scaling == 2:
        v = 1
        a = 1.0/gamma

    b = 0.5*beta**2
    L = 1.0
    ymin = 0
    # Need gloal to be able change ymax in closure process_u
    global ymax
    ymax = 1.2

    I = lambda x: 0
    f = lambda x, t: delta*np.exp(-b*(x - v*t)**2)

    import time
    import scitools.std as plt
    plot_arrays = []

    def process_u(u, x, t, n):
        global ymax
        if animate:
            plt.plot(x, u, 'r-',
                     x, f(x, t[n])/delta, 'b-',
                     axis=[0, L, ymin, ymax], title='t=%f' % t[n],
                     xlabel='x', ylabel='u and f/%g' % delta)
        if t[n] == 0:
            time.sleep(1)
            plot_arrays.append(x)
        dt = t[1] - t[0]
        tol = dt/10.0
        if abs(t[n] - 0.2) < tol or abs(t[n] - 0.5) < tol:
            plot_arrays.append((u.copy(), f(x, t[n])/delta))
            if u.max() > ymax:
                ymax = u.max()

    Nx = 100
    D = 10
    T = 0.5
    u_L = u_R = 0
    theta = 1.0
    cpu = solver(
        I, a, f, L, Nx, D, T, theta, u_L, u_R, user_action=process_u)
    x = plot_arrays[0]
    plt.figure()
    for u, f in plot_arrays[1:]:
        plt.plot(x, u, 'r-', x, f, 'b--', axis=[x[0], x[-1], 0, ymax],
                 xlabel='$x$', ylabel=r'$u, \ f/%g$' % delta)
        plt.hold('on')
    plt.legend(['$u,\\ t=0.2$', '$f/%g,\\ t=0.2$' % delta,
                '$u,\\ t=0.5$', '$f/%g,\\ t=0.5$' % delta])
    filename = 'tmp1_gamma%g_s%d' % (gamma, scaling)
    s = 'diffusion' if scaling == 1 else 'source'
    plt.title(r'$\beta = %g,\ \gamma = %g,\ $' % (beta, gamma)
              + 'scaling=%s' % s)
    plt.savefig(filename + '.pdf');  plt.savefig(filename + '.png')
    return cpu

def investigate():
    """Do scienfic experiments with the run function above."""
    # Clean up old files
    import glob
    for filename in glob.glob('tmp1_gamma*') + \
            glob.glob('welding_gamma*'):
        os.remove(filename)

    gamma_values = 1, 40, 5, 0.2, 0.025
    for gamma in gamma_values:
        for scaling in 1, 2:
            run(gamma=gamma, beta=10, delta=20, scaling=scaling)

    # Combine images
    for gamma in gamma_values:
        for ext in 'pdf', 'png':
            cmd = 'doconce combine_images -2 '\
                  'tmp1_gamma%(gamma)g_s1.%(ext)s '\
                  'tmp1_gamma%(gamma)g_s2.%(ext)s '\
                  'welding_gamma%(gamma)g.%(ext)s' % vars()
            os.system(cmd)
            # pdflatex doesn't like 0.2 in filenames...
            if '.' in str(gamma):
                os.rename(
                    'welding_gamma%(gamma)g.%(ext)s' % vars(),
                    ('welding_gamma%(gamma)g' % vars()).replace('.', '_')
                    + '.' + ext)

if __name__ == '__main__':
    #run(gamma=1/40., beta=10, delta=40, scaling=2)
    investigate()
