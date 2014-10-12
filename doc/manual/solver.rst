.. role:: envvar(literal)
.. role:: command(literal)
.. role:: file(literal)
.. _solver:


Solver
======

In this section, we show how to use the predefined linear solvers, *scipy* and *petsc*. The use of *Pastix* and *Murge* will be
explained in another section.


Predefined solvers
^^^^^^^^^^^^^^^^^^

Pigasus offers predefined solvers like *Conjugate-Gradient*, *Gauss-Seidel*, *etc*. in order to import the *solver* class, the user can simply call::

  from pigasus.solver.solver import *

  
.. todo:: Give examples for CG/GS/.../PASTIX


Using Scipy
^^^^^^^^^^^

The following example lists the use of *scipy.sparse.linalg* solvers.

.. literalinclude:: include/demo/test_solver_scipy.py
    :linenos:
    :language: python

Which gives the following results

.. literalinclude:: include/demo/test_solver_scipy.txt

Using PyAMG
^^^^^^^^^^^


.. literalinclude:: include/demo/test_solver_pyamg.py
    :linenos:
    :language: python

Which gives the following results

.. literalinclude:: include/demo/test_solver_pyamg.txt

Using Petsc 
^^^^^^^^^^^


.. literalinclude:: include/demo/test_solver_petsc4py.py
    :linenos:
    :language: python

Which gives the following results

.. literalinclude:: include/demo/test_solver_petsc4py.txt

Geometric MultiGrid
^^^^^^^^^^^^^^^^^^^

**Pigasus** includes a Geometric MultiGrid. More details can be found in see :ref:`multigrid`.


.. Local Variables:
.. mode: rst
.. End:
