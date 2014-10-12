.. role:: envvar(literal)
.. role:: command(literal)
.. role:: file(literal)
.. _boundcond:


Boundary Conditions
*******************

In the current version of **Pigasus**, we only handle *Dirichlet, Neumann* and *periodic* boundary conditions. In this section, we
show how to implement them.

The treatment of Dirichlet boundary conditions must be done carefuly. There is a difference between Homogeneous Dirichlet boundary conditions and non homogeneous ones. The first one can be treated as::

  # Homogeneous Dirichlet boundary condition can be done in two ways
  # using the flag AllDirichlet
  testcase['AllDirichlet'] = True
  # by giving all external faces, using the geometry object 
  testcase['Dirichlet'] = geo.external_faces
  # by giving specific faces, for each patch
  testcase['Dirichlet'] = [[0,3], [1,2]]

Non Homogeneous Dirichlet doundary conditions can be given as a dictionary with a key specified by the *patch-id* and the *face-id*::

  bc_dirichlet={}
  bc_dirichlet [0,0] = g 
  ...
  # save the boundary condition in the testcase dictionary
  testcase['bc_dirichlet'] = bc_dirichlet

.. note::

  The latter dictionary takes a vectorial function object even if the value at the boundary is only a fixed scalar. 

A *1D* example

.. literalinclude:: include/demo/test_dirichlet_1d.py
    :linenos:
    :language: python
    :lines: 14-74

The complete script can be found here

.. plot:: include/demo/test_dirichlet_1d.py

Which gives the numerical solution

.. image:: include/demo/test_dirichlet_1d.png

A *2D* example

.. literalinclude:: include/demo/test_dirichlet_2d.py
    :linenos:
    :language: python
    :lines: 14-88

The complete script can be found here

.. plot:: include/demo/test_dirichlet_2d.py

Which gives the numerical solution

.. image:: include/demo/test_dirichlet_2d.png


In the case of Neumann boundary conditions, the user can do it like for non homogeneous Dirichlet case::

  bc_neumann = {}
  bc_neumann [0,0] = g 
  ...
  # save the boundary condition in the testcase dictionary
  testcase['bc_neumann'] = bc_neumann

.. note::

  The Neumann boundary condition must be given as a vectorial function.

We give here some examples depending on the geometry

.. literalinclude:: include/demo/test_neumann_2d.py
    :linenos:
    :language: python
    :lines: 84-181

The complete script can be found here

.. plot:: include/demo/test_neumann_2d.py

Which gives the numerical solution

.. image:: include/demo/test_neumann_2d.png

On a circle domain

.. literalinclude:: include/demo/test_neumann_circle.py
    :linenos:
    :language: python
    :lines: 53-122

The complete script can be found here

.. plot:: include/demo/test_neumann_circle.py

Which gives the numerical solution

.. image:: include/demo/test_neumann_circle.png

On a quart-circle domain

.. literalinclude:: include/demo/test_neumann_quartcircle.py
    :linenos:
    :language: python
    :lines: 72-131

The complete script can be found here

.. plot:: include/demo/test_neumann_quartcircle.py

Which gives the numerical solution

.. image:: include/demo/test_neumann_quartcircle.png

.. todo:: Give an example for Periodic BC and the connectivity function

.. Local Variables:
.. mode: rst
.. End:
