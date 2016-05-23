"""Illustration of staggered mesh in time. Two mesh functions."""
from pysketcher import *

u = SketchyFunc3()
v = SketchyFunc1()
Nt = 5
t_mesh = linspace(0, 6, Nt+1)
dt = t_mesh[1] - t_mesh[0]

# Add 20% space to the left and 30% to the right of the coordinate system
t_axis_extent = t_mesh[-1] - t_mesh[0]
t_min = t_mesh[0] - 0.2*t_axis_extent
t_max = t_mesh[-1] + 0.3*t_axis_extent
#u_max = 1.3*max([u(t) for t in t_mesh])
u_max = 1.1*max([u(t) for t in t_mesh])
u_min = -0.75*u_max

drawing_tool.set_coordinate_system(t_min, t_max, u_min, 1.2*u_max,
                                   axis=True, xkcd=False)
drawing_tool.set_linecolor('black')

r = 0.005*(t_max-t_min)     # radius of circles placed at mesh points
#import random; random.seed(12)
perturbations = [0, 0.1, -0.4, 0.2, -0.4, -0.1]
u_points = {}
u_values = []
v_points = {}
v_values = []
for i, t in enumerate(t_mesh):
    u_value = u(t) + perturbations[i]
    u_values.append(u_value)
    u_points[i] = Composition(dict(
        circle=Circle(point(t, u_value), r).set_filled_curves('black'),
        u_point=Text('$u^%d$' % i,
                     point(t, u_value) + (point(0,3*r)
                                          if i > 0 else point(-3*r,0)))))
    if 0 <= i < len(t_mesh)-1:
        v_value = v(t+dt/2) - perturbations[i]
        v_values.append(v_value)
        v_points[i] = Composition(dict(
            circle=Circle(point(t+dt/2, v_value), r).set_linecolor('black'),
            #v_point=Text(r'$v^%d+\frac{1}{2}}$' % i,
            v_point=Text(r'$v^{%d/2}$' % (2*i+1),
                         point(t+dt/2, v_value) + (point(0,3*r)
                                              if i > 0 else point(-3*r,0)))))
u_discrete = Composition(u_points)
v_discrete = Composition(v_points)

b = 0  # Spacing below original x axis (spacing needed for code)

axes = Composition(dict(
    x=Axis(point(0,-b), t_mesh[-1] + 0.2*t_axis_extent, '$t$',
           label_spacing=(1/45.,-1/30.)),
    y=Axis(point(0,-b), 1.1*u_max, '',
           rotation_angle=90)))

h = 0.03*u_max  # tickmarks height
u_nodes = Composition({i: Composition(dict(
    node=Line(point(t,h-b), point(t,-b-h))
    #name=Text('$t_%d$' % i, point(t,-b-3.5*h)),
    )) for i, t in enumerate(t_mesh)})
h /= 2
v_nodes = Composition({i: Composition(dict(
    node=Line(point(t+dt/2,h-b), point(t+dt/2,-b-h))
    #name=Text('$t_%d$' % i, point(t,-b-3.5*h)),
    )) for i, t in enumerate(t_mesh)})

illustration = Composition(dict(u=u_discrete,
                                v=v_discrete,
                                mesh_u=u_nodes,
                                mesh_v=v_nodes,
                                axes=axes)).set_name('staggered')
drawing_tool.erase()
# Draw t_mesh with discrete u points
illustration.draw()

# Add exact u line (u is a Spline Shape that applies 500 intervals by default
# for drawing the curve)
exact_u = u.set_linestyle('dashed').set_linewidth(1)
exact_v = v.set_linestyle('dashed').set_linewidth(1).set_linecolor('blue')
exact = Composition(dict(u=exact_u, v=exact_v))
#exact.draw()
#drawing_tool.display()

drawing_tool.display()
drawing_tool.savefig('tmp1')
input()
