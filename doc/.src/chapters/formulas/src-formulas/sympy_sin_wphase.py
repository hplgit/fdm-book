from sympy import *


t, a, Q, w = symbols('t a A w', positive=True, real=True)
u = Q*exp(-a*t)*sin(w*t)
dudt = diff(u, t)
dudt
factor(dudt)
simplify(dudt)
# Alternative, manually derived expression
phi = atan(-a/w)
A = Q*sqrt(a**2 + w**2)
dudt2 = exp(-a*t)*A*cos(w*t - phi)

simplify(expand_trig(dudt2))
simplify(expand_trig(dudt2 - dudt))  # are they equal?
s = solve(dudt2, t)
s
