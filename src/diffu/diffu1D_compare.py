"""Compare FE, BE, and CN."""
from diffu1D_u0 import solver_theta, solver_FE, solver_BE
F = float(sys.argv[1])
dt = 0.0002
u = {} # solutions for all schemes
for name in 'FE', 'BE', 'theta':
    u[name] = plug(solver='solver_'+name, F=F, dt=dt)
# Note that all schemes employ the same time mesh regardless of F
x = u['FE'][0]
t = u['FE'][1]
for n in range(len(t)):
    plot(x, u['FE'][2+n], 'r-',
         x, u['BE'][2+n], 'b-',
         x, u['theta'][2+n], 'g-',
         legend=['FE', 'BE', 'CN'],
         xlabel='x', ylabel='u',
         title='t=%f' % t[n],
         savefig='tmp_frame%04d.png' % n)
