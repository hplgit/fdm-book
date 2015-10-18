from sympy import *
import sys

print '1D analysis:'
x, t, A, B, w, b, c, k, q = symbols('x t A B w b c k q')
u = (A*cos(w*t) + B*sin(w*t))*exp(-c*t)*cos(k*x)

# Constrain B from the initial condition u_t=0
u_t = diff(u, t)
u_t = simplify(u_t.subs(t, 0))
print 'u_t(x,0)=0:', u_t
# B = A*c/w
u = u.subs(B, A*c/w)
print 'u:', u

R = diff(u, t, t) + b*diff(u, t) - q*diff(u, x, x)
R = simplify(R)
terms = collect(R, [cos(w*t), sin(w*t)])
print 'factored terms:', terms
print 'latex terms:', latex(terms, mode='plain')
#  A*((-c**2*w + k**2*q*w - w**3)*cos(t*w) + (-b*c**2 - b*w**2 + c**3 + c*k**2*q + c*w**2)*sin(t*w))*exp(-c*t)*cos(k*x)/w
cos_eq = -c**2*w + k**2*q*w - w**3
sin_eq = -b*c**2 - b*w**2 + c**3 + c*k**2*q + c*w**2
print 'cos_eq has solution', solve(cos_eq, w)
# w = sqrt(k**2*q - c**2), assume c < k*sqrt(q)
sin_eq = sin_eq.subs(w, sqrt(k**2*q - c**2))
sin_eq = simplify(sin_eq)
print 'sin_eq:', sin_eq
print 'sin_eq has solution', solve(sin_eq, c)
# c = b/2

print '2D analysis:'
y, kx, ky = symbols('y kx ky')
u = (A*cos(w*t) + B*sin(w*t))*exp(-c*t)*cos(kx*x)*cos(ky*y)

# Constrain B from the initial condition u_t=0
u_t = diff(u, t)
u_t = simplify(u_t.subs(t, 0))
print 'u_t(x,0)=0:', u_t
# B = A*c/w
u = u.subs(B, A*c/w)
print 'u:', u

R = diff(u, t, t) + b*diff(u, t) - q*diff(u, x, x) - q*diff(u, y, y)
R = simplify(R)
terms = collect(R, [cos(w*t), sin(w*t)])
print 'factored terms 2D:', terms
print 'latex terms 2D:', latex(terms, mode='plain')
#  A*((-c**2*w + kx**2*q*w + ky**2*q*w - w**3)*cos(t*w) + (-b*c**2 - b*w**2 + c**3 + c*kx**2*q + c*ky**2*q + c*w**2)*sin(t*w))*exp(-c*t)*cos(kx*x)*cos(ky*y)/w
cos_eq = -c**2*w + kx**2*q*w + ky**2*q*w - w**3
sin_eq = -b*c**2 - b*w**2 + c**3 + c*kx**2*q + c*ky**2*q + c*w**2
print 'cos_eq has solution', solve(cos_eq, w)
# w = sqrt(k**2*q - c**2), assume c < k*sqrt(q)
sin_eq = sin_eq.subs(w, sqrt(kx**2*q + ky**2*q - c**2))
sin_eq = simplify(sin_eq)
print 'sin_eq:', sin_eq
print 'sin_eq, b and c terms:', collect(sin_eq, [b, c])
print 'sin_eq has solution', solve(sin_eq, c)

