.. role:: envvar(literal)
.. role:: command(literal)
.. role:: file(literal)
.. _gallery.poisson:

Poisson
*******

Let us consider the Poisson equation

*For given functions* :math:`f,g,k`, *find* :math:`u` *such that*

.. math::

  - \nabla^2 u &= f , ~~~~ \Omega
     \\
     u &= g  , ~~~~ \Gamma_D 
     \\
     \nabla u \cdot \mathbf{n} &= k, ~~~~ \Gamma_N    

You can use the *basicPDE* class to write your own Poisson solver, or use the pre-implemented *poisson* class solver

Example of the *1D* solver

.. literalinclude:: include/demo/test_poisson_1d.py
    :linenos:
    :language: python
    :lines: 75-89

See :download:`the complete script <include/demo/test_poisson_1d.py>`.

the resulting numerical solution is

.. image:: include/demo/test_poisson_1d.png
   :width: 9cm
   :height: 9cm

Example of the *2D* solver

.. literalinclude:: include/demo/test_poisson_2d.py
    :linenos:
    :language: python
    :lines: 96-109

See :download:`the complete script <include/demo/test_poisson_2d.py>`.

the resulting numerical solution is

.. image:: include/demo/test_poisson_2d.png
   :width: 9cm
   :height: 9cm

An example of solving the Poisson equation on a circle domain is given here :download:`script <include/demo/test_poisson_circle.py>`.


the resulting numerical solution is

.. image:: include/demo/test_poisson_circle.png
   :width: 9cm
   :height: 9cm

The user can define an external *Metric* to be associated to the discrete functional space. An example is given here :download:`script <include/demo/test_poisson_circle_metric.py>`.

.. Local Variables:
.. mode: rst
.. End:
