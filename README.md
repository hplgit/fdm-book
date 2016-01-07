## fdm-book

Resources for the book *The Craft of Finite Difference Computing with Partial Differential Equations* by Hans Petter Langtangen and Svein Linge.

> This easy-to-read book introduces the basics of solving partial differential
> equations by finite difference methods. The emphasis is on constructing
> finite difference schemes, formulating algorithms, implementing
> algorithms, verifying implementations, analyzing the physical behavior
> of the numerical solutions, and applying the methods and software
> to solve problems from physics and biology.



## Directory structure

 * `src`: source code for book examples
 * `src/X`: source code from chapter `X`
 * `doc`: documents for the book
 * `doc/pub`: (preliminary) [published documents](http://hplgit/github.io/fdm-book/doc/web/index.html)
 * `doc/pub/book`: complete (preliminary) published book
 * `doc/pub/X`: (preliminary) published chapter `X`
 * `doc/.src`: DocOnce source code and computer code
 * `doc/.src/book`: DocOnce source code for assembling the complete book
 * `doc/.src/chapters/X`: DocOnce source code for chapter `X`

Source files are compiled in `doc/.src/chapters/X` (a specific chapter)
or `/doc/.src/book` (the entire book) and copied to some subdirectory
under `doc/pub` for publication in the gh-pages branch.

Development takes place in the master branch, and the gh-pages branch
is always merge with master.

