#!/usr/bin/env python
"""
2D wave equation solved by finite differences::

  dt, cpu_time = solver(I, V, f, c, Lx, Ly, Nx, Ny, dt, T,
                        user_action=None, version='scalar',
                        stability_safety_factor=1)

Solve the 2D wave equation u_tt = u_xx + u_yy + f(x,t) on (0,L) with
u=0 on the boundary and initial condition du/dt=0.

Nx and Ny are the total number of mesh cells in the x and y
directions. The mesh points are numbered as (0,0), (1,0), (2,0),
..., (Nx,0), (0,1), (1,1), ..., (Nx, Ny).

dt is the time step. If dt<=0, an optimal time step is used.
T is the stop time for the simulation.

I, V, f are functions: I(x,y), V(x,y), f(x,y,t). V and f
can be specified as None or 0, resulting in V=0 and f=0.

user_action: function of (u, x, y, t, n) called at each time
level (x and y are one-dimensional coordinate vectors).
This function allows the calling code to plot the solution,
compute errors, etc.
"""
import time, sys, os
from glob import glob
import scitools.std as st
import numpy as np
try:
    import mayavi.mlab as mlab
except:
    # We don't have mayavi
    pass

def solver(I, V, f, c, Lx, Ly, Nx, Ny, dt, T,
           user_action=None, version='scalar'):
    if version == 'vectorized':
        advance = advance_vectorized
    elif version == 'scalar':
        advance = advance_scalar

    x = np.linspace(0, Lx, Nx+1)  # Mesh points in x dir
    y = np.linspace(0, Ly, Ny+1)  # Mesh points in y dir
    # Make sure dx, dy are compatible with x, y
    dx = x[1] - x[0]
    dy = y[1] - y[0]

    xv = x[:,np.newaxis]          # For vectorized function evaluations
    yv = y[np.newaxis,:]

    stability_limit = (1/float(c))*(1/np.sqrt(1/dx**2 + 1/dy**2))
    if dt <= 0:                # max time step?
        safety_factor = -dt    # use negative dt as safety factor
        dt = safety_factor*stability_limit
    elif dt > stability_limit:
        print 'error: dt=%g exceeds the stability limit %g' % \
              (dt, stability_limit)
    Nt = int(round(T/float(dt)))
    t = np.linspace(0, Nt*dt, Nt+1)    # mesh points in time
    dt = t[1] - t[0]                   # ensure compatibility
    Cx2 = (c*dt/dx)**2;  Cy2 = (c*dt/dy)**2    # help variables
    dt2 = dt**2

    # Allow f and V to be None or 0
    if f is None or f == 0:
        f = (lambda x, y, t: 0) if version == 'scalar' else \
            lambda x, y, t: np.zeros((x.shape[0], y.shape[1]))
        # or simpler: x*y*0
    if V is None or V == 0:
        V = (lambda x, y: 0) if version == 'scalar' else \
            lambda x, y: np.zeros((x.shape[0], y.shape[1]))


    order = 'Fortran' if version == 'f77' else 'C'
    u     = np.zeros((Nx+1,Ny+1), order=order)   # Solution array
    u_n   = np.zeros((Nx+1,Ny+1), order=order)   # Solution at t-dt
    u_nm1 = np.zeros((Nx+1,Ny+1), order=order)   # Solution at t-2*dt
    f_a   = np.zeros((Nx+1,Ny+1), order=order)   # For compiled loops

    Ix = range(0, u.shape[0])
    Iy = range(0, u.shape[1])
    It = range(0, t.shape[0])

    import time; t0 = time.clock()  # For measuring CPU time

    # Load initial condition into u_n
    if version == 'scalar':
        for i in Ix:
            for j in Iy:
                u_n[i,j] = I(x[i], y[j])
    else:
        # Use vectorized version (requires I to be vectorized)
        u_n[:,:] = I(xv, yv)

    if user_action is not None:
        user_action(u_n, x, xv, y, yv, t, 0)

    # Special formula for first time step
    n = 0
    # First step requires a special formula, use either the scalar
    # or vectorized version (the impact of more efficient loops than
    # in advance_vectorized is small as this is only one step)
    if version == 'scalar':
        u = advance_scalar(
            u, u_n, u_nm1, f, x, y, t, n,
            Cx2, Cy2, dt2, V, step1=True)

    else:
        f_a[:,:] = f(xv, yv, t[n])  # precompute, size as u
        V_a = V(xv, yv)
        u = advance_vectorized(
            u, u_n, u_nm1, f_a,
            Cx2, Cy2, dt2, V=V_a, step1=True)

    if user_action is not None:
        user_action(u, x, xv, y, yv, t, 1)

    # Update data structures for next step
    #u_nm1[:] = u_n;  u_n[:] = u  # safe, but slow
    u_nm1, u_n, u = u_n, u, u_nm1

    for n in It[1:-1]:
        if version == 'scalar':
            # use f(x,y,t) function
            u = advance(u, u_n, u_nm1, f, x, y, t, n, Cx2, Cy2, dt2)
        else:
            f_a[:,:] = f(xv, yv, t[n])  # precompute, size as u
            u = advance(u, u_n, u_nm1, f_a, Cx2, Cy2, dt2)

        if user_action is not None:
            if user_action(u, x, xv, y, yv, t, n+1):
                break

        # Update data structures for next step
        u_nm1, u_n, u = u_n, u, u_nm1

    # Important to set u = u_n if u is to be returned!
    t1 = time.clock()
    # dt might be computed in this function so return the value
    return dt, t1 - t0



def advance_scalar(u, u_n, u_nm1, f, x, y, t, n, Cx2, Cy2, dt2,
                   V=None, step1=False):
    Ix = range(0, u.shape[0]);  Iy = range(0, u.shape[1])
    if step1:
        dt = sqrt(dt2)  # save
        Cx2 = 0.5*Cx2;  Cy2 = 0.5*Cy2; dt2 = 0.5*dt2  # redefine
        D1 = 1;  D2 = 0
    else:
        D1 = 2;  D2 = 1
    for i in Ix[1:-1]:
        for j in Iy[1:-1]:
            u_xx = u_n[i-1,j] - 2*u_n[i,j] + u_n[i+1,j]
            u_yy = u_n[i,j-1] - 2*u_n[i,j] + u_n[i,j+1]
            u[i,j] = D1*u_n[i,j] - D2*u_nm1[i,j] + \
                     Cx2*u_xx + Cy2*u_yy + dt2*f(x[i], y[j], t[n])
            if step1:
                u[i,j] += dt*V(x[i], y[j])
    # Boundary condition u=0
    j = Iy[0]
    for i in Ix: u[i,j] = 0
    j = Iy[-1]
    for i in Ix: u[i,j] = 0
    i = Ix[0]
    for j in Iy: u[i,j] = 0
    i = Ix[-1]
    for j in Iy: u[i,j] = 0
    return u

def advance_vectorized(u, u_n, u_nm1, f_a, Cx2, Cy2, dt2,
                       V=None, step1=False):
    if step1:
        dt = np.sqrt(dt2)  # save
        Cx2 = 0.5*Cx2;  Cy2 = 0.5*Cy2; dt2 = 0.5*dt2  # redefine
        D1 = 1;  D2 = 0
    else:
        D1 = 2;  D2 = 1
    u_xx = u_n[:-2,1:-1] - 2*u_n[1:-1,1:-1] + u_n[2:,1:-1]
    u_yy = u_n[1:-1,:-2] - 2*u_n[1:-1,1:-1] + u_n[1:-1,2:]
    u[1:-1,1:-1] = D1*u_n[1:-1,1:-1] - D2*u_nm1[1:-1,1:-1] + \
                   Cx2*u_xx + Cy2*u_yy + dt2*f_a[1:-1,1:-1]
    if step1:
        u[1:-1,1:-1] += dt*V[1:-1, 1:-1]
    # Boundary condition u=0
    j = 0
    u[:,j] = 0
    j = u.shape[1]-1
    u[:,j] = 0
    i = 0
    u[i,:] = 0
    i = u.shape[0]-1
    u[i,:] = 0
    return u

def quadratic(Nx, Ny, version):
    """Exact discrete solution of the scheme."""

    def exact_solution(x, y, t):
        return x*(Lx - x)*y*(Ly - y)*(1 + 0.5*t)

    def I(x, y):
        return exact_solution(x, y, 0)

    def V(x, y):
        return 0.5*exact_solution(x, y, 0)

    def f(x, y, t):
        return 2*c**2*(1 + 0.5*t)*(y*(Ly - y) + x*(Lx - x))

    Lx = 5;  Ly = 2
    c = 1.5
    dt = -1 # use longest possible steps
    T = 18

    def assert_no_error(u, x, xv, y, yv, t, n):
        u_e = exact_solution(xv, yv, t[n])
        diff = abs(u - u_e).max()
        tol = 1E-12
        msg = 'diff=%g, step %d, time=%g' % (diff, n, t[n])
        assert diff < tol, msg

    new_dt, cpu = solver(
        I, V, f, c, Lx, Ly, Nx, Ny, dt, T,
        user_action=assert_no_error, version=version)
    return new_dt, cpu


def test_quadratic():
    # Test a series of meshes where Nx > Ny and Nx < Ny
    versions = 'scalar', 'vectorized', 'cython', 'f77', 'c_cy', 'c_f2py'
    for Nx in range(2, 6, 2):
        for Ny in range(2, 6, 2):
            for version in versions:
                print 'testing', version, 'for %dx%d mesh' % (Nx, Ny)
                quadratic(Nx, Ny, version)

def run_efficiency(nrefinements=4):
    from numpy import pi, sin

    def I(x, y):
        return sin(pi*x/Lx)*np.sin(pi*y/Ly)

    Lx = 10;  Ly = 10
    c = 1.5
    T = 100
    versions = ['scalar', 'vectorized', 'cython', 'f77',
               'c_f2py', 'c_cy']
    print ' '*15, ''.join(['%-13s' % v for v in versions])
    for Nx in 15, 30, 60, 120:
        cpu = {}
        for version in versions:
            dt, cpu_ = solver(I, None, None, c, Lx, Ly, Nx, Nx,
                              -1, T, user_action=None,
                              version=version)
            cpu[version] = cpu_
        cpu_min = min(list(cpu.values()))
        if cpu_min < 1E-6:
            print 'Ignored %dx%d grid (too small execution time)' \
                  % (Nx, Nx)
        else:
            cpu = {version: cpu[version]/cpu_min for version in cpu}
            print '%-15s' % '%dx%d' % (Nx, Nx),
            print ''.join(['%13.1f' % cpu[version] for version in versions])

def gaussian(plot_method=2, version='vectorized', save_plot=True):
    """
    Initial Gaussian bell in the middle of the domain.
    plot_method=1 applies mesh function, =2 means surf, =0 means no plot.
    """
    # Clean up plot files
    for name in glob('tmp_*.png'):
        os.remove(name)

    Lx = 10
    Ly = 10
    c = 1.0

    from numpy import exp

    def I(x, y):
        """Gaussian peak at (Lx/2, Ly/2)."""
        return exp(-0.5*(x-Lx/2.0)**2 - 0.5*(y-Ly/2.0)**2)

    if plot_method == 3:
        from mpl_toolkits.mplot3d import axes3d
        import matplotlib.pyplot as plt
        from matplotlib import cm
        plt.ion()
        fig = plt.figure()
        u_surf = None

    def plot_u(u, x, xv, y, yv, t, n):
        """User action function for plotting."""
        if t[n] == 0:
            time.sleep(2)
        if plot_method == 1:
            # Works well with Gnuplot backend, not with Matplotlib
            st.mesh(x, y, u, title='t=%g' % t[n], zlim=[-1,1],
                    caxis=[-1,1])
        elif plot_method == 2:
            # Works well with Gnuplot backend, not with Matplotlib
            st.surfc(xv, yv, u, title='t=%g' % t[n], zlim=[-1, 1],
                  colorbar=True, colormap=st.hot(), caxis=[-1,1],
                  shading='flat')
        elif plot_method == 3:
            print 'Experimental 3D matplotlib...under development...'
            # Probably too slow
            #plt.clf()
            ax = fig.add_subplot(111, projection='3d')
            u_surf = ax.plot_surface(xv, yv, u, alpha=0.3)
            #ax.contourf(xv, yv, u, zdir='z', offset=-100, cmap=cm.coolwarm)
            #ax.set_zlim(-1, 1)
            # Remove old surface before drawing
            if u_surf is not None:
                ax.collections.remove(u_surf)
            plt.draw()
            time.sleep(1)
        elif plot_method == 4:
	    # Mayavi visualization
            mlab.clf()
            extent1 = (0, 20, 0, 20,-2, 2)
            s = mlab.surf(x , y, u,
                          colormap='Blues',
                          warp_scale=5,extent=extent1)
            mlab.axes(s, color=(.7, .7, .7), extent=extent1,
                      ranges=(0, 10, 0, 10, -1, 1),
                      xlabel='', ylabel='', zlabel='',
                      x_axis_visibility=False,
                      z_axis_visibility=False)
            mlab.outline(s, color=(0.7, .7, .7), extent=extent1)
            mlab.text(6, -2.5, '', z=-4, width=0.14)
            mlab.colorbar(object=None, title=None,
                          orientation='horizontal',
                          nb_labels=None, nb_colors=None,
                          label_fmt=None)
            mlab.title('Gaussian t=%g' % t[n])
            mlab.view(142, -72, 50)
            f = mlab.gcf()
            camera = f.scene.camera
            camera.yaw(0)

        if plot_method > 0:
            time.sleep(0) # pause between frames
            if save_plot:
                filename = 'tmp_%04d.png' % n
		if plot_method == 4:
                    mlab.savefig(filename)  # time consuming!
		elif plot_method in (1,2):
                    st.savefig(filename)  # time consuming!

    Nx = 40; Ny = 40; T = 20
    dt, cpu = solver(I, None, None, c, Lx, Ly, Nx, Ny, -1, T,
                     user_action=plot_u, version=version)



if __name__ == '__main__':
    #test_quadratic()
    gaussian(plot_method=2, version='vectorized', save_plot=True)
