import sympy as sym
def pde_residual(u, q):
    f = sym.diff(u, t, t) - sym.diff(q*sym.diff(u, x), x)
    #f = sym.simplify(sym.expand(f))
    f = sym.simplify(f)
    return f

x, t, L = sym.symbols('x t L')

q = 1 + (x-L/2)**4
u = sym.cos(sym.pi*x/L)*sym.cos(sym.pi/L*t)
f = pde_residual(u, q)
print f

q = 1 + sym.cos(sym.pi*x/L)
u = sym.cos(sym.pi*x/L)*sym.cos(sym.pi/L*t)
f = pde_residual(u, q)
print f
