"""
Illustrate a 2D spatial mesh with mesh points.
"""

from pysketcher import *

xaxis = 2
xmin = ymin = 0; xmax = ymax = 10
drawing_tool.set_coordinate_system(xmin-1, xmax+1.5, ymin-1, ymax+1, axis=False)
drawing_tool.set_linecolor('black')
drawing_tool.set_linewidth(1)

Nx = 3
Ny = 4
Nx = 3
Ny = 2
dx = (xmax - xmin)/float(Nx)
dy = (ymax - ymin)/float(Ny)

# Draw mesh lines
mesh = Composition({})
for i in range(Nx+1):
    mesh['line_x%0d2' % i] = Line((xmin + i*dx, ymin), (xmin + i*dx, ymax))
for j in range(Ny+1):
    mesh['line_y%0d2' % j] = Line((xmin, ymin + j*dy), (xmax, ymin + j*dy))

# Draw (i,j) and global unknown number
N = (Nx+1)*(Ny+1)
numbering = Composition({})
A = []  # sparse matrix representation
index = lambda i,j: j*(Nx+1) + i
for j in range(Ny+1):
    for i in range(Nx+1):
        x = xmin + i*dx;  y = ymin + j*dy
        ij = index(i,j)
        numbering['pt(%d,%d)' % (i,j)] = \
              Text('(%d,%d): %d' % (i,j, ij), (x+0.05, y+0.05), 'left', fontsize=9)
    for i in range(Nx+1):
        ij = index(i,j)
        A.append([])
        if A[-1]:
            print 'A[%d]' % ij, 'is not empty!'
        if i != 0 and i != Nx and j != 0 and j != Ny:
            if j > 0:
                A[-1].append(((i,j-1), index(i,j-1)))
            if i > 0:
                A[-1].append(((i-1,j), index(i-1,j)))
        A[-1].append(((i,j), ij))
        if i != 0 and i != Nx and j != 0 and j != Ny:
            if i < Nx:
                A[-1].append(((i+1,j), index(i+1,j)))
            if j < Ny:
                A[-1].append(((i,j+1), index(i,j+1)))

# Show matrix as latex code
code = code_x = r"""
\left(\begin{array}{%s}
""" % ('c'*N, )
for ij in range(len(A)):
    line1 = []
    line2 = []
    nz = [A[ij][i][1] for i in range(len(A[ij]))]  # nonzeros
    for ji in range(len(A)):
        if ji in nz:
            line1.append('(%d,%d)' % (ij,ji))
            line2.append('\\bullet')
        else:
            line1.append('%s 0%s ' % (' '*(len(str(ij))), ' '*(len(str(ji)))))
            line2.append('\\cdot')
    code += ' & '.join(line1)
    code_x += ' & '.join(line2)
    if i < len(A)-1:
        code += ' \\\\' + '\n'
        code_x += ' \\\\' + '\n'
code += r'\end{array}\right)' + '\n'
code_x += r'\end{array}\right)' + '\n'
with open('tmp.tex', 'w') as f:
    f.write(code)
    f.write('\n\n')
    f.write(code_x)

fig = Composition(dict(mesh=mesh, numbering=numbering))
drawing_tool.erase()
fig.draw()
drawing_tool.display()
drawing_tool.savefig('tmp1')

raw_input()
