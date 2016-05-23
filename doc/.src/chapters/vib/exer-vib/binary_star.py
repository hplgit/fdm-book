#import matplotlib.pyplot as plt
import scitools.std as plt
import odespy
import numpy as np

def solver(alpha, ic, T, dt=0.05):
    def f(u, t):
        x_A, vx_A, y_A, vy_A, x_B, vx_B, y_B, vy_B = u
        distance3 = np.sqrt((x_B-x_A)**2 + (y_B-y_A)**2)**3
        system = [
            vx_A,
            1/(1.0 + alpha)*(x_B - x_A)/distance3,
            vy_A,
            1/(1.0 + alpha)*(y_B - y_A)/distance3,
            vx_B,
            -1/(1.0 + alpha**(-1))*(x_B - x_A)/distance3,
            vy_B,
            -1/(1.0 + alpha**(-1))*(y_B - y_A)/distance3,
            ]
        return system

    Nt = int(round(T/dt))
    t_mesh = np.linspace(0, Nt*dt, Nt+1)

    solver = odespy.RK4(f)
    solver.set_initial_condition(ic)
    u, t = solver.solve(t_mesh)
    x_A = u[:,0]
    x_B = u[:,2]
    y_A = u[:,4]
    y_B = u[:,6]
    return x_A, x_B, y_A, y_B, t

def demo_circular():
    # Mass B is at rest at the origin,
    # mass A is at (1, 0) with vel. (0, 1)
    ic = [1, 0, 0, 1, 0, 0, 0, 0]
    x_A, x_B, y_A, y_B, t = solver(
        alpha=0.001, ic=ic, T=2*np.pi, dt=0.01)
    plt.plot(x_A, x_B, 'r2-', y_A, y_B, 'b2-',
             legend=['A', 'B'],
             daspectmode='equal') # x and y axis have same scaling
    plt.savefig('tmp_circular.png')
    plt.savefig('tmp_circular.pdf')
    plt.show()

def demo_two_stars(animate=True):
    # Initial condition
    ic = [0.6, 0, 0, 1,   # star A: velocity (0,1)
          0, 0, 0, -0.5]  # star B: velocity (0,-0.5)
    # Solve ODEs
    x_A, x_B, y_A, y_B, t = solver(
        alpha=0.5, ic=ic, T=4*np.pi, dt=0.05)
    if animate:
        # Animate motion and draw the objects' paths in time
        for i in range(len(x_A)):
            plt.plot(x_A[:i+1], x_B[:i+1], 'r-',
                     y_A[:i+1], y_B[:i+1], 'b-',
                     [x_A[0], x_A[i]], [x_B[0], x_B[i]], 'r2o',
                     [y_A[0], y_A[i]], [y_B[0], y_B[i]], 'b4o',
                     daspectmode='equal',  # axes aspect
                     legend=['A', 'B', 'A', 'B'],
                     axis=[-1, 1, -1, 1],
                     savefig='tmp_%04d.png' % i,
                     title='t=%.2f' % t[i])
    else:
        # Make a simple static plot of the solution
        plt.plot(x_A, x_B, 'r-', y_A, y_B, 'b-',
                 daspectmode='equal', legend=['A', 'B'],
                 axis=[-1, 1, -1, 1], savefig='tmp_two_stars.png')
    #plt.axes().set_aspect('equal')  # mpl
    plt.show()

if __name__ == '__main__':
    import sys
    if sys.argv[1] == 'circular':
        demo_circular()
    else:
        demo_two_stars(True)
    raw_input()
