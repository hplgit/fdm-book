"""
Explore algebraic forms arising from the integral f(u)*v in the finite
element method.

   |              phi_im1            phi_i            phi_ip1
   +1                /\               /\                /\
   |                /  \             /  \              /  \
   |               /    \           /    \            /    \
   |              /      \         /      \          /      \
   |             /        \       /        \        /        \
   |            /          \     /          \      /          \
   |           /            \   /            \    /            \
   |          /              \ /              \  /              \
   |         /                /                \/                \
   |        /                / \               /\                 \
   |       /                /   \             /  \                 \
   |      /                /     \           /    \                 \
   |     /                /       \         /      \                 \
   |    /                /         \       /        \                 \
   |   /                /           \     /          \                 \
   |  /                /             \   /            \                 \
   | /                /               \ /              \                 \
--------------------------------------------------------------------------
                     i-1               i               i+1

                          cell L             cell R
"""
from sympy import *
import sys

x, u_im1, u_i, u_ip1, u, h, x_i = symbols('x u_im1 u_i u_ip1 u h x_i')

# Left cell:  [x_im1, x_i]
# Right cell: [x_i, x_ip1]
x_im1 = x_i - h
x_ip1 = x_i + h

phi = {'L':  # Left cell
       {'im1': 1 - (x - x_im1)/h,
        'i': (x - x_im1)/h},
       'R':  # Right cell
       {'i': 1 - (x-x_i)/h,
        'ip1': (x - x_i)/h}}

u = {'L': u_im1*phi['L']['im1'] + u_i  *phi['L']['i'],
     'R': u_i  *phi['R']['i']   + u_ip1*phi['R']['ip1']}

f = lambda u: eval(sys.argv[1])

integral_L = integrate(f(u['L'])*phi['L']['i'], (x, x_im1, x_i))
integral_R = integrate(f(u['R'])*phi['R']['i'], (x, x_i, x_ip1))
expr_i = simplify(expand(integral_L + integral_R))
print expr_i
latex_code = latex(expr_i, mode='plain')
# Replace u_im1 sympy symbol name by latex symbol u_{i-1}
latex_code = latex_code.replace('im1', '{i-1}')
# Replace u_ip1 sympy symbol name by latex symbol u_{i+1}
latex_code = latex_code.replace('ip1', '{i+1}')
print latex_code
# Escape (quote) latex_code so it can be sent as HTML text
import cgi
html_code = cgi.escape(latex_code)
print html_code
# Make a file with HTML code for displaying the LaTeX formula
f = open('tmp.html', 'w')
# Include an image that can be clicked on to yield a new
# page with an interactive editor and display area where the
# formula can be further edited
text = """
<a href="http://www.codecogs.com/eqnedit.php?latex=%(html_code)s"
 target="_blank">
<img src="http://latex.codecogs.com/gif.latex?%(html_code)s"
 title="%(latex_code)s"/>
</a>
 """ % vars()
f.write(text)
f.close()
# load tmp.html into a browser


