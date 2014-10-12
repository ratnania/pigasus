.. role:: envvar(literal)
.. role:: command(literal)
.. role:: file(literal)
.. _gallery.parabolic:


Parabolic Partial Differential Equations
****************************************

Pigasus offers a pre-defined *parabolic* solver. For the moment, only a one-step method has been implemented.

Let us consider the parabolic equation

*For a given function* :math:`f`, *find* :math:`u` *such that*

.. math::

     \partial_t u &=  \mathcal{L} (u), ~~~~ \Omega
     \\
     u &= g  , ~~~~ \Gamma_D 
     \\
     \nabla u \cdot \mathbf{n} &= k, ~~~~ \Gamma_N 

Where :math:`\mathcal{L}` is the Differential Operator described by the *basicPDE* class.

This parabolic equation can be descritized in time using a one-step method, and can be written in the form 

.. math:: 

  \mathcal{I} u^{n+1} = \mathcal{E} u^{n}

The following example shows how to solve the Anisotropic Diffusion problem.

.. include:: gallery/anidiff.rst
   :width: 9cm
   :height: 9cm


.. Local Variables:
.. mode: rst
.. End:
