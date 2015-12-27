"""As oscillator_sketch1.py, but without wheels."""
from pysketcher import *

L = 12.
H = L/6
W = L/6

xmax = L
drawing_tool.set_coordinate_system(xmin=-L, xmax=xmax,
                                   ymin=-1, ymax=L+H,
                                   axis=False,
                                   instruction_file='tmp_mpl.py')
x = 0
drawing_tool.set_linecolor('black')

def make_dashpot(x):
    d_start = (-L,2*H)
    d = Dashpot(start=d_start, total_length=L+x, width=W,
                bar_length=3*H/2, dashpot_length=L/2, piston_pos=H+x)
    d.rotate(-90, d_start)
    return d

def make_spring(x):
    s_start = (-L,4*H)
    s = Spring(start=s_start, length=L+x, bar_length=3*H/2, teeth=True)
    s.rotate(-90, s_start)
    return s

d = make_dashpot(0)
s = make_spring(0)

M = Rectangle((0,H), 4*H, 4*H).set_linewidth(4)
left_wall = Rectangle((-L,H),H/10,L-H).set_filled_curves(pattern='/')
ground = Wall(x=[-L/2,L], y=[H,H], thickness=-H/10)

fontsize = 18
text_m = Text('$m$', (2*H, H+2*H), fontsize=fontsize)
text_kx = Text('$s(u)$', (-L/2, H+4*H), fontsize=fontsize)
text_bv = Text('$f(u)$', (-1.25*L/2, H), fontsize=fontsize)
x_axis = Axis((2*H, L), H, '$u(t)$', fontsize=fontsize,
              label_spacing=(0.04, -0.01))
x_axis_start = Line((2*H, L-H/4), (2*H, L+H/4)).set_linewidth(4)

fig = Composition({
    'spring': s, 'mass': M, 'left wall': left_wall,
    'ground': ground,
    'text_m': text_m, 'text_kx': text_kx,
    'x_axis': x_axis, 'x_axis_start': x_axis_start})

fig.draw()
drawing_tool.display()
drawing_tool.savefig('tmp_oscillator_sliding')

drawing_tool.erase()

fig['dashpot'] = d
fig['text_bv'] = text_bv

# or fig = Composition(dict(fig=fig, dashpot=d, text_bv=text_bv))
fig.draw()

drawing_tool.display()
drawing_tool.savefig('tmp_oscillator_sliding_dashpot')

drawing_tool.erase()

F_force = Force((4*H, H+2*H), (4*H+H, H+2*H), '$F(t)$',
                text_spacing=(0.057, -0.007), text_alignment='left', fontsize=fontsize)
fig['F_force'] = F_force
fig.draw()
drawing_tool.savefig('tmp_oscillator_sliding_dashpot_force')

input()
