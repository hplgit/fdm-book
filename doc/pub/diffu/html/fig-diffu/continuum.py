"""
Mapping of arbitrary continnum volume as often needed when
deriving PDEs in continuum mechanics.
"""
from pysketcher import *

H = 7.
W = 5.

drawing_tool.set_coordinate_system(xmin=0, xmax=W,
                                   ymin=2, ymax=H,
                                   axis=False)
#drawing_tool.set_grid(True)
drawing_tool.set_linecolor('black')

volume1 = ArbitraryVolume((2*W/10, H/2,), vector_field_symbol='q',
                          width=W/4)

fig = Composition({'volume1': volume1})
fig.draw()
drawing_tool.display()
drawing_tool.savefig('tmp1')

input()
