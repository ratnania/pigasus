pigasus
=======

Pigasus is a generic Python package for solving (system of) Partial Differential Equations. Its core is written in Fortran. The aim of Pigasus is to discretize spatial differential operators, which makes it easier to write different time schemes thanks to the oriented object aspect of Python.

The basic distribution of Pigasus offers some examples (see the gallery module). Other solvers classes can be adressed as Plugins. These classes will be shared only on request.

For the beta version, MPI and OpenCL were deactivated, as their validation is not yet finished.
