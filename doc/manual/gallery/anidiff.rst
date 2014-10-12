.. role:: envvar(literal)
.. role:: command(literal)
.. role:: file(literal)
.. _gallery.anidiff:


Anisotropic Diffusion
^^^^^^^^^^^^^^^^^^^^^

In this example, we solve the 2D anisotropic Diffusion problem, with Homogeneous Dirichlet boundary conditions, using Pigasus.

.. math::

  \partial_t u &=  \nabla \cdot \left( A_{\epsilon} \nabla u \right), ~~~~ \Omega
  \\
  u &= 0, ~~~~  \partial \Omega

Which can be descritized in time as

.. math:: 

  u^{n+1} - (1-\alpha) dt \nabla \cdot \left( A_{\epsilon} \nabla u^{n+1}\right) = u^{n} + \alpha dt \nabla \cdot \left( A_{\epsilon} \nabla u^{n}\right)  

By introducing the Operators

.. math::

  I &: u \rightarrow u - (1-\alpha) dt \nabla \cdot \left( A_{\epsilon} \nabla u\right)
  \\
  E &: u \rightarrow u + \alpha dt \nabla \cdot \left( A_{\epsilon} \nabla u\right)


The latter equation becomes

.. math::

  Iu^{n+1} = Eu^{n}
  

Using Pigasus, we'll need to define two *basicPDE* objects, one for each Operator::

  E = basicPDE(geometry=geo, testcase=tc_E)
  I = basicPDE(geometry=geo, testcase=tc_I)

* Assembly the Explicit and Implicit operators **E** and **I**::

     E.assembly() 
     I.assembly() 

* Define the *fields* :math:`u^{n}` and :math:`u^{n+1}` as **un** and **unew**::

     un   = E.unknown
     unew = I.unknown

* Get the values of the field at initialization::

     u0 = E.rhs
     un.set(u0)

* Hence, the algorithm for solving the heat equation becomes ::

     for i in range(0,niter):
         rhs = E.dot(un)
         I.solve(rhs) # updates unew automatically
         un.set(unew)

The user can also use the one-step method from the *parabolic* class as the following

.. literalinclude:: include/demo/test_anisotropicDiffusion.py
    :linenos:
    :language: python

Which gives the numerical solution

.. image:: include/demo/test_anisotropicDiffusion.png
   :width: 9cm
   :height: 9cm


.. Local Variables:
.. mode: rst
.. End:
