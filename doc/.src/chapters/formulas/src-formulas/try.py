# Method of manufactured solutions (MMS)

from sympy import *
x, t, rho, dt = symbols('x t rho dt')

def a(u):
    return 1 + u**2

def u_simple(x, t):
    return x**2*(Rational(1,2) - x/3)*t

# Show that u_simple satisfies the BCs
for x_point in 0, 1:
    print 'u_x(%s,t):' % x_point,
    print diff(u_simple(x, t), x).subs(x, x_point).simplify()

print 'Initial condition:', u_simple(x, 0)

# MMS: full nonlinear problem
u = u_simple(x, t)
f = rho*diff(u, t) - diff(a(u)*diff(u, x), x)
print 'f: nonlinear problem:'
print f.simplify()

# MMS: one-Picard-iteration linearized problem
u_1 = u_simple(x, t-dt)
f = rho*diff(u, t) - diff(a(u_1)*diff(u, x), x)
print 'f: 1-Picard-iteration linearized problem:'
print f.simplify()

# MMS: Backward Euler discretization in time +
# one-Picard-iteration linearized problem
def D_t_backward(u):
    return (u(x, t) - u(x, t - dt))/dt

u_1 = u_simple(x, t-dt)
f = rho*D_t_backward(u_simple) - diff(a(u_1)*diff(u, x), x)
print 'f: BE in time + 1-Picard-iteration linearized problem:'
print f.simplify()
