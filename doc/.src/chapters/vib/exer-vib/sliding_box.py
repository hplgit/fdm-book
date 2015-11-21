import matplotlib.pyplot as plt
import scitools.std as plt
import numpy as np

def plot_spring():
    alpha_values = [1, 2, 3, 10]
    s = lambda u: 1.0/alpha*np.tanh(alpha*u)
    u = np.linspace(-1, 1, 1001)
    for alpha in alpha_values:
        print alpha, s(u)
        plt.plot(u, s(u))
        plt.hold('on')
    plt.legend([r'$\alpha=%g$' % alpha for alpha in alpha_values])
    plt.xlabel('u');  plt.ylabel('Spring response $s(u)$')
    plt.savefig('tmp_s.png'); plt.savefig('tmp_s.pdf')

def simulate(beta, gamma, delta=0,
             num_periods=8, time_steps_per_period=60):
    # Use oscillations without friction to set dt and T
    P = 2*np.pi
    dt = P/time_steps_per_period
    T = num_periods*P
    t = np.linspace(0, T, time_steps_per_period*num_periods+1)
    import odespy
    def f(u, t, beta, gamma):
        # Note the sequence of unknowns: v, u (v=du/dt)
        v, u = u
        return [-beta*np.sign(v) - 1.0/gamma*np.tanh(gamma*u), v]
        #return [-beta*np.sign(v) - u, v]

    solver = odespy.RK4(f, f_args=(beta, gamma))
    solver.set_initial_condition([delta, 1]) # sequence must match f
    uv, t = solver.solve(t)
    u = uv[:,1]  # recall sequence in f: v, u
    v = uv[:,0]
    return u, t


if __name__ == '__main__':
    plt.figure()
    plot_spring()

    beta_values = [0, 0.05, 0.1]
    gamma_values = [0.1, 1, 5]
    for gamma in gamma_values:
        plt.figure()
        for beta in beta_values:
            u, t = simulate(beta, gamma, 0, 6, 60)
            plt.plot(t, u)
            plt.hold('on')
        plt.legend([r'$\beta=%g$' % beta for beta in beta_values])
        plt.title(r'$\gamma=%g$' % gamma)
        plt.xlabel('$t$');  plt.ylabel('$u$')
        filestem = 'tmp_u_gamma%g' % gamma
        plt.savefig(filestem + '.png')
        plt.savefig(filestem + '.pdf')
    plt.show()
    raw_input()
