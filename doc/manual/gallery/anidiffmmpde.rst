.. role:: envvar(literal)
.. role:: command(literal)
.. role:: file(literal)
.. _gallery.anidiffmmpde:


Moving FEM for Anistropic Diffusion
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In this section, we consider the discretization of a general parabolic problem using Moving Finite Element Method. The mesh will
be generated using the Equidistribution principle, by solving the Monge-Amp\`ere problem, as described previously.

We consider again the equation

.. math::

  \partial_t u - \nabla \cdot \left( A_{\epsilon} \nabla u \right) = f, ~~~~ \Omega
  \\
  u &= 0, ~~~~  \partial \Omega

where we have added a time dependent term source :math:`f:=f(\mathbf{x}, t)`.

because of the time dependence of the mesh, the time derivative must e changed

.. math::

  \partial_t u = \dot{u} + \dot{x} \partial_x u + \dot{y} \partial_y u 

Which can be descritized in time as

.. math:: 

  u^{n+1} - (1-\alpha) dt \nabla \cdot \left( A_{\epsilon} \nabla u^{n+1}\right) + \frac{1-\beta}{2} \dot{x} \partial_x u^{n+1} + \frac{1-\beta}{2} \dot{y} \partial_y u^{n+1} 
  
  = u^{n} + \alpha dt \nabla \cdot \left( A_{\epsilon} \nabla u^{n}\right) - \frac{\beta}{2} \dot{x} \partial_x u^{n} - \frac{\beta}{2} \dot{y} \partial_y u^{n}    


By introducing the Operators

.. math::

  I &: u \rightarrow u - (1-\alpha) dt \nabla \cdot \left( A_{\epsilon} \nabla u\right) + \frac{1-\beta}{2} \dot{x} \partial_x u + \frac{1-\beta}{2} \dot{y} \partial_y u
  \\
  E &: u \rightarrow u + \alpha dt \nabla \cdot \left( A_{\epsilon} \nabla u\right) - \frac{\beta}{2} \dot{x} \partial_x u - \frac{\beta}{2} \dot{y} \partial_y u


The latter equation becomes

.. math::

  Iu^{n+1} = Eu^{n}

Both :math:`\dot{x}` and :math:`\dot{y}` will be computed as

.. math::

  \dot{\mathbf{x}} = \frac{\mathbf{x}^{n+1} - \mathbf{x}^n}{dt}

This leads to the final algorithm ::

  while n < niter:
        updateMesh()
        rhs = E.dot(un)
        I.assembly()
        I.solve(rhs)
        unew.set(un)

        n += 1 ; t += dt

The routine **updateMesh** computes the new mapping :math:`\mathbf{x}^{n+1} = F^{n+1}(\boldsymbol\xi)`

.. note::

  We can use different time steps for discretizing the problem :ref:`MMPDE_anidiff` and the mesh equation. A general form for the
  routine **updateMesh** is the following::

    # mesh generation is done each nfreq step
    if n % nfreq > 0:
          return

    solveMongeAmpere()


.. Local Variables:
.. mode: rst
.. End:
