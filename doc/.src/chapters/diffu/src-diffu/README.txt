Sparse matrices: Need 1D *and* 2D example. Find 2D examples around
with sparse mat before implementing the real 5620 solver.

Try scipy.
Try pysparse (python-pysparse, see examples, pysparse.sourceforge.net).
Try PyAMG (see examples, ask Luke).
Try petsc4py: petsc4py-tutorial has nice examples (vec, mat, distributed
arrays, stencils!) See also examples kspsolve/petsc-mat.py in petsc4py
source tree.

Might be smart to create a little common interface after the three
examples above, to allow switching and seeing what the key operations are.
Nice class, good for advanced scientific programming (Sverker/Uppsala).

Also:
see SfePy for a finite element code in Python.
http://docs.sfepy.org/doc-devel/index.html
http://labs.freehackers.org/projects/pythonmeshviewer
