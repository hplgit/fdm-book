# Test effect of vectorization
from diffu1D_u0 import solver_FE_simple, solver_FE
I = lambda x: 1
Nx = 100000
a = 2.0
L = 2.0
dx = L/Nx
dt = dx**2/(2*a)
T = 100*dt
u, x, t, cpu1 = solver_FE_simple(
    I=I, a=a, f=None, L=L, Nx=Nx, F=0.5, T=T)
cpu1
u, x, t, cpu2 = solver_FE(
    I=I, a=a, f=None, L=L, Nx=Nx, F=0.5, T=T, version='scalar')
u, x, t, cpu3 = solver_FE(
    I=I, a=a, f=None, L=L, Nx=Nx, F=0.5, T=T, version='vectorized')
cpu2/cpu3
