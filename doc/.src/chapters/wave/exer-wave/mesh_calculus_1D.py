# -*- coding: utf-8 -*-
"""
Calculus with a 1D mesh function.
"""
import numpy as np

class MeshCalculus:
    def __init__(self, vectorized=True):
        self.vectorized = vectorized

    def differentiate(self, f, x):
        '''
        Computes the derivative of f by centered differences, but 
        forw and back difference at the start and end, respectively.
        '''
        dx = x[1] - x[0]
        Nx = len(x) - 1     # number of spatial steps
        num_dfdx = np.zeros(Nx+1)
        # Compute approximate derivatives at end-points first
        num_dfdx[0] = (f(x[1]) - f(x[0]))/dx          # FD approx.
        num_dfdx[Nx] = (f(x[Nx]) - f(x[Nx-1]))/dx     # BD approx.
        # proceed with approximate derivatives for inner mesh points
        if self.vectorized:
            num_dfdx[1:-1] = (f(x[2:]) - f(x[:-2]))/(2*dx)
        else:   # scalar version
            for i in range(1, Nx):
                num_dfdx[i] = (f(x[i+1]) - f(x[i-1]))/(2*dx)
        return num_dfdx
    
    def integrate(self, f, x):
        '''
        Computes the integral of f(x) over the interval 
        covered by x. 
        '''
        dx = x[1] - x[0]
        F = np.zeros(len(x))   
        F[0] = 0    # starting value for iterative scheme
        if self.vectorized:
            all_trapezoids = np.zeros(len(x)-1)    
            all_trapezoids[:] = 0.5*(f(x[:-1]) + f(x[1:]))*dx   
            F[1:] = np.cumsum(all_trapezoids)
        else:   # scalar version
            for i in range(0, len(x)-1):
                F[i+1] = F[i] + 0.5*(f(x[i]) + f(x[i+1]))*dx    
        return F
    
def test_differentiate():
    def f(x):
        return 4*x - 2.5
    def dfdx(x):
        derivatives = np.zeros(len(x))
        derivatives[:] = 4
        return derivatives
        
    a = 0; b = 1; Nx = 10
    x = np.linspace(a, b, Nx+1)    
    exact_dfdx = dfdx(x)        # compute exact derivatives
    # test vectorized version
    calc_v = MeshCalculus(vectorized=True)
    num_dfdx  = calc_v.differentiate(f, x)
    print np.abs(num_dfdx - exact_dfdx)
    diff = np.abs(num_dfdx - exact_dfdx).max()
    tol = 1E-14
    assert diff < tol
    # test scalar version
    calc = MeshCalculus(vectorized=False)
    num_dfdx  = calc.differentiate(f, x)
    print np.abs(num_dfdx - exact_dfdx)
    diff = np.abs(num_dfdx - exact_dfdx).max()
    assert diff < tol
    
def test_integrate():
    def f(x):
        return 4*x - 2.5
    a = 0; b = 1; Nx = 10
#    a = 2.5/4; b = 10; Nx = 2
    x = np.linspace(a, b, Nx+1)    
    # The exact integral amounts to the total area of two triangles
    I_exact = 0.5*abs(2.5/4 - a)*f(a) + 0.5*abs(b - 2.5/4)*f(b)
    # test vectorized version    
    calc_v = MeshCalculus(vectorized=True)    
    F = calc_v.integrate(f, x)
    print F, I_exact
    diff = np.abs(F[-1] - I_exact)
    print diff
    tol = 1E-14
    assert diff < tol    
    # test scalar version
    calc = MeshCalculus(vectorized=False)        
    F = calc.integrate(f, x)
    print F, I_exact
    diff = np.abs(F[-1] - I_exact)
    print diff
    assert diff < tol
    
    
if __name__ == '__main__':
    test_differentiate()
    test_integrate()
