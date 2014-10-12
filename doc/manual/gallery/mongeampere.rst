.. role:: envvar(literal)
.. role:: command(literal)
.. role:: file(literal)
.. _gallery.mongeampere:


MongeAmpere
***********

Here, we solve the Monge-Ampere equation using the *Benamou-Froese-Oberman* method (**BFO**)

The Monge-Ampere equation (MA) writes:

.. math::
   :label: MA-eq

   \rho_1(\nabla \phi) \det H(\phi) = \rho_0(\boldsymbol\xi),~ ~ ~ ~ \boldsymbol\xi \in \Omega_c

where :math:`H(\phi)` is the Hessian matrix of the potential :math:`\phi`.

Let us define the operator :math:`T: H^2(\Omega) \rightarrow H^2(\Omega)` by

.. math::

   T[u] = \left( \nabla^{2}\right)^{-1} \sqrt{(\nabla^2 u )^2 + 2 (f - \det{H(u)})}

which in **2D** writes

.. math::

   T[u] = \left( \nabla^{2}\right)^{-1} \sqrt{u_{xx}^2 + u_{yy}^2  + 2 u_{xy}^2 + 2 f }

It is proved in :cite:`Benamou2010a` that when :math:`u \in H^2(\Omega)` is a solution of the Monge-Amp\`ere equation, then it is a fixed point of the operator :math:`T`. 

.. math::
  :label: MA-eq-BFO

   u = T[u]


Picard Algorithm
^^^^^^^^^^^^^^^^

In :cite:`Benamou2010a`, Picard's method is used to solve :ref:`MA-eq-BFO`.

The algorithm is

* Given an initial value :math:`u^0`,
* Compute :math:`u^{n+1}` as the solution of 

.. math::
   :label: MA-eq-Benamou

   - \nabla^2 u^{n+1} = - \sqrt{(\nabla^2 u^n )^2 + 2 (f - \det{H(u^n)})} 

.. literalinclude:: include/demo/test_monge_ampere_picard.py
    :linenos:
    :language: python
    :lines: 34-113

See :download:`the complete script <include/demo/test_monge_ampere_picard.py>`.

.. image:: include/demo/test_monge_ampere_picard.png
   :width: 15cm
   :height: 9cm

.. Local Variables:
.. mode: rst
.. End:
