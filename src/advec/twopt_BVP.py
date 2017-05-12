import numpy as np

def solver(eps, Nx, method='centered'):
    """
    Solver for the two point boundary value problem u'=eps*u'',
    u(0)=0, u(1)=1.
    """
    x = np.linspace(0, 1, Nx+1)       # Mesh points in space
    # Make sure dx and dt are compatible with x and t
    dx = x[1] - x[0]
    u   = np.zeros(Nx+1)

    # Representation of sparse matrix and right-hand side
    diagonal = np.zeros(Nx+1)
    lower    = np.zeros(Nx)
    upper    = np.zeros(Nx)
    b        = np.zeros(Nx+1)

    # Precompute sparse matrix (scipy format)
    if method == 'centered':
        diagonal[:] = 2*eps/dx**2
        lower[:] = -1/dx - eps/dx**2
        upper[:] =  1/dx - eps/dx**2
    elif method == 'upwind':
        diagonal[:] = 1/dx + 2*eps/dx**2
        lower[:] =  1/dx - eps/dx**2
        upper[:] = - eps/dx**2

    # Insert boundary conditions
    upper[0] = 0
    lower[-1] = 0
    diagonal[0] = diagonal[-1] = 1
    b[-1] = 1.0

    # Set up sparse matrix and solve
    diags = [0, -1, 1]
    import scipy.sparse
    import scipy.sparse.linalg
    A = scipy.sparse.diags(
        diagonals=[diagonal, lower, upper],
        offsets=[0, -1, 1], shape=(Nx+1, Nx+1),
        format='csr')
    u[:] = scipy.sparse.linalg.spsolve(A, b)
    return u, x

def u_exact(x, eps):
    return (np.exp(x/eps)-1)/(np.exp(1.0/eps)-1)

def demo(eps = 0.01, method='centered'):
    import matplotlib.pyplot as plt
    x_fine = np.linspace(0, 1, 2001)
    for Nx in (20, 40):
        u, x = solver(eps, Nx, method=method)
        plt.figure()
        plt.plot(x, u, 'o-')
        plt.hold('on')
        plt.plot(x_fine, u_exact(x_fine, eps), 'k--')
        plt.legend(['$N_x=%d$' % Nx, 'exact'], loc='upper left')
        plt.title(method + ' difference scheme, ' + r'$\epsilon=%g$' % eps)
        plt.xlabel('x');  plt.ylabel('u')
        stem = 'tmp1_%s_%d_%s' % (method, Nx, str(eps).replace('.','_'))
        plt.savefig(stem + '.png'); plt.savefig(stem + '.pdf')
    plt.show()

if __name__ == '__main__':
    demo(eps=0.1, method='upwind')
    demo(eps=0.01, method='upwind')
    #demo(eps=0.1, method='centered')
    #demo(eps=0.01, mehtod='centered')
