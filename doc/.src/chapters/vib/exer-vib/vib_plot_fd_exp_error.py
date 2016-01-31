import matplotlib.pyplot as plt
import numpy as np
import sympy as sym

def E_fraction(p):
    return (2./p)**2*(np.sin(p/2.))**2

a = 0; b = np.pi
p = np.linspace(a, b, 100)
E_values = np.zeros(len(p))

# create 4th degree Taylor polynomial (also plotted)
p_ = sym.symbols('p_')
E = (2./p_)**2*(sym.sin(p_/2.))**2
E_series = E.series(p_, 0, 4).removeO()
print E_series
E_pyfunc = sym.lambdify([p_], E_series, modules='numpy')

# To avoid division by zero when p is 0, we rather take the limit
E_values[0] = sym.limit(E, p_, 0, dir='+')  # ...when p --> 0, E --> 1
E_values[1:] = E_fraction(p[1:])

plt.plot(p, E_values, 'k-', p, E_pyfunc(p), 'k--')
plt.xlabel('p'); plt.ylabel('Error fraction')
plt.legend(['E', 'E Taylor'])
plt.savefig('tmp_error_fraction.png')
plt.savefig('tmp_error_fraction.pdf')
plt.show()
