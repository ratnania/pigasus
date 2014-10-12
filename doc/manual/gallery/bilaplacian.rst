.. role:: envvar(literal)
.. role:: command(literal)
.. role:: file(literal)
.. _gallery.bilaplacian:


Bilaplacian
***********

Pigasus offers a pre-implemented *Bilaplacian* solver.

Let us consider the Biharmonic equation

*For a given function* :math:`f`, *find* :math:`u` *such that*

.. math::

  - \nabla^4 u &= f , ~~~~ \Omega
     \\
     u &= 0  , ~~~~ \partial \Omega

.. literalinclude:: include/demo/test_bilaplacian_2d.py
    :linenos:
    :language: python
    :lines: 54-65

See :download:`the complete script <include/demo/test_bilaplacian_2d.py>`.

The numerical solution is then

.. image:: include/demo/test_bilaplacian_2d.png
   :width: 9cm
   :height: 9cm

.. Local Variables:
.. mode: rst
.. End:
