import sys, os
sys.path.insert(0, os.path.join(os.pardir, 'src-wave', 'wave1D'))
from wave1D_n0 import solver
import nose.tools as nt

def test_cubic2():
    import sympy as sym
    x, t, c, L, dx = sym.symbols('x t c L dx')
    T = lambda t: 1 + sym.Rational(1,2)*t  # Temporal term
    # Set u as a 3rd-degree polynomial in space
    X = lambda x: sum(a[i]*x**i for i in range(4))
    a = sym.symbols('a_0 a_1 a_2 a_3')
    u = lambda x, t: X(x)*T(t)
    # Force discrete boundary condition to be zero by adding
    # a correction term the analytical suggestion x*(L-x)*T
    # u_x = x*(L-x)*T(t) - 1/6*u_xxx*dx**2
    R = sym.diff(u(x,t), x) - (
        x*(L-x) - sym.Rational(1,6)*sym.diff(u(x,t), x, x, x)*dx**2)
    # R is a polynomial: force all coefficients to vanish.
    # Turn R to Poly to extract coefficients:
    R = sym.poly(R, x)
    coeff = R.all_coeffs()
    s = sym.solve(coeff, a[1:])  # a[0] is not present in R
    # s is dictionary with a[i] as keys
    # Fix a[0] as 1
    s[a[0]] = 1
    X = lambda x: sym.simplify(sum(s[a[i]]*x**i for i in range(4)))
    u = lambda x, t: X(x)*T(t)
    print 'u:', u(x,t)
    # Find source term
    f = sym.diff(u(x,t), t, t) - c**2*sym.diff(u(x,t), x, x)
    f = sym.simplify(f)
    print 'f:', f

    u_exact = sym.lambdify([x,t,L,dx], u(x,t), 'numpy')
    V_exact = sym.lambdify([x,t,L,dx], sym.diff(u(x,t), t), 'numpy')
    f_exact = sym.lambdify([x,t,L,dx,c], f)

    # Replace symbolic variables by numeric ones
    L = 2.0
    Nx = 3
    C = 0.75
    c = 0.5
    dt = C*(L/Nx)/c
    dx = dt*c/C

    I = lambda x: u_exact(x, 0, L, dx)
    V = lambda x: V_exact(x, 0, L, dx)
    f = lambda x, t: f_exact(x, t, L, dx, c)

    # user_action function asserts correct solution
    def assert_no_error(u, x, t, n):
        u_e = u_exact(x, t[n], L, dx)
        diff = abs(u - u_e).max()
        nt.assert_almost_equal(diff, 0, places=13)

    u, x, t, cpu = solver(
        I=I, V=V, f=f, c=c, L=L,
        dt=dt, C=C, T=4*dt, user_action=assert_no_error)

def test_cubic1():
    import sympy as sym
    x, t, c, L, dx, dt = sym.symbols('x t c L dx dt')
    i, n = sym.symbols('i n', integer=True)

    # Assume discrete solution is a polynomial of degree 3 in x
    T = lambda t: 1 + sym.Rational(1,2)*t  # Temporal term
    a = sym.symbols('a_0 a_1 a_2 a_3')
    X = lambda x: sum(a[q]*x**q for q in range(4))  # Spatial term
    u = lambda x, t: X(x)*T(t)

    DxDx = lambda u, i, n: (u((i-1)*dx, n*dt) - 2*u(i*dx, n*dt) +
                            u((i+1)*dx, n*dt))/dx**2
    DtDt = lambda u, i, n: (u(i*dx, (n-1)*dt) - 2*u(i*dx, n*dt) +
                            u(i*dx, (n+1)*dt))/dt**2
    D2x  = lambda u, i, n: (u((i+1)*dx, n*dt) - u((i-1)*dx, n*dt)) \
                            /(2*dx)
    R_0 = sym.simplify(D2x(u, 0, n))   # residual du/dx, x=0
    Nx = L/dx
    R_L = sym.simplify(D2x(u, Nx, n))  # residual du/dx, x=L
    print R_0
    print R_L
    # We have two equations, let a_0 and a_1 be free parameters,
    # adjust a_2 and a_3 so that R_0=0 and R_L=0.
    # For simplicity in final expressions, set a_0=0, a_1=1.
    R_0 = R_0.subs(a[0], 0).subs(a[1], 1)
    R_L = R_L.subs(a[0], 0).subs(a[1], 1)
    a = list(a)  # enable in-place assignment
    a[0:2] = 0, 1
    s = sym.solve([R_0, R_L], a[2:])
    print s
    a[2:] = s[a[2]], s[a[3]]
    # Calling X(x) will now use new a since a has changed
    print 'u:', u(x, t)
    print 'DxDx(x**3,i,n)', sym.simplify(DxDx(lambda x, t: x**3, i, n))
    f = DtDt(u, i, n) - c**2*DxDx(u, i, n)
    f = sym.expand(f)  # Easier to simplify if expanded first
    f = f.subs(i, x/dx).subs(n, t/dt)
    f = sym.simplify(f)
    print 'f:', f

    u_exact = sym.lambdify([x,t,L,dx], u(x,t), 'numpy')
    V_exact = sym.lambdify([x,t,L,dx], sym.diff(u(x,t), t), 'numpy')
    f_exact = sym.lambdify([x,t,L,dx,dt,c], f, 'numpy')

    # Replace symbolic variables by numeric ones
    L = 2.0
    Nx = 3
    C = 0.75
    c = 0.5
    dt = C*(L/Nx)/c
    dx = dt*c/C

    I = lambda x: u_exact(x, 0, L, dx)
    V = lambda x: V_exact(x, 0, L, dx)
    f = lambda x, t: f_exact(x, t, L, dx, dt, c)

    # user_action function asserts correct solution
    def assert_no_error(u, x, t, n):
        u_e = u_exact(x, t[n], L, dx)
        diff = abs(u - u_e).max()
        nt.assert_almost_equal(diff, 0, places=13)

    u, x, t, cpu = solver(
        I=I, V=V, f=f, c=c, L=L,
        dt=dt, C=C, T=4*dt, user_action=assert_no_error)

test_cubic1()
test_cubic2()
