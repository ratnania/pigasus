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

Anisotropic Diffusion eigenvalue problem
****************************************

In this example, we solve the 2D anisotropic Diffusion steady state problem.

.. math::

  - \nabla \cdot \left( \mathbf{b} \otimes \mathbf{b} \nabla u \right) &= \omega^2 u, ~~~~ \Omega
  \\
  u &= 0, ~~~~  \partial \Omega

Naive formulation
^^^^^^^^^^^^^^^^^

The weak formulation writes

.. math::

   \sum_j [u]^j \int ( \mathbf{b} \cdot \nabla \phi_i ) ( \mathbf{b} \cdot \nabla \phi_j ) = \omega^2 \sum_j [u]^j \int \phi_i \phi_j 

which leads to the linear eigenvalue problem

.. math::

  \mathcal{S} [u] = \omega^2 \mathcal{M} [u]
  
Mixed formulation
^^^^^^^^^^^^^^^^^

Let us introduce the auxiliary vector function variable :math:`\mathbf{p}` such that

.. math::

  p_{\parallel} := \mathbf{b} \cdot \mathbf{p} = \mathbf{b} \cdot \nabla u

In this case, our annistropic diffusion eigenvalue problem writes

.. math::

  \mathbf{b} \cdot \mathbf{p} = \mathbf{b} \cdot \nabla u
  \\
  - \nabla \cdot ( \mathbf{b} p_{\parallel} ) = \omega^2 u

We assume that :math:`u \in V` and :math:`\mathbf{p} \in W` where :math:`V` and :math:`W` are two spaces that we will introduce later. 

.. math::

  \sum_j [\mathbf{p}]^j \int (\mathbf{b} \cdot \Psi_j) (\mathbf{b} \cdot \Psi_i) = \sum_j [u]^j \int (\mathbf{b} \cdot \nabla \phi_j) (\mathbf{b} \cdot \Psi_i)
  \\
  \sum_j [\mathbf{p}]^j( \mathbf{b} \cdot \Psi_j ) ( \mathbf{b} \cdot \nabla \phi_i) = \omega^2 \sum_j [u]^j \int \phi_i \phi_j   

By introducing the matrices

.. math::

  \mathcal{M}_W = (\int ( \mathbf{b} \cdot \Psi_j) ( \mathbf{b} \cdot \Psi_i) )_{1 \leq i,j \leq n_W}
  \\
  \mathcal{M}_V = (\int \phi_j \phi_i )_{1 \leq i,j \leq n_V}
  \\
  \mathcal{K}_{WV} = (\int (\mathbf{b} \cdot \nabla \phi_j) ( \mathbf{b} \cdot \Psi_i) )_{1 \leq i \leq n_W, 1 \leq j \leq n_V}
  \\ 
  \mathcal{K}_{VW} = \mathcal{K}_{WV}^T  

The linear system writes

.. math:: 

  \mathcal{M}_W [\mathbf{p}] = \mathcal{K}_{WV} [u]
  \\
  \mathcal{K}_{VW} [\mathbf{p}] = \omega^2 \mathcal{M}_{V} [u]

which can be written as

.. math::

  \mathcal{A} \begin{pmatrix}
  [\mathbf{p}] \\
  [u] 
  \end{pmatrix}
  = \omega^2
  \mathcal{M} \begin{pmatrix}
  [\mathbf{p}] \\
  [u] 
  \end{pmatrix}

with

.. math::

  \mathcal{A} = \begin{pmatrix}
   -\mathcal{M}_W  &  \mathcal{K}_{WV} \\
  \mathcal{K}_{VW}   & 0 
  \end{pmatrix}

and  

.. math::

  \mathcal{M} = \begin{pmatrix}
   0 & 0 \\
   0 & \mathcal{M}_V 
  \end{pmatrix}

Standard discretization
_______________________

In this case we consider that :math:`\mathbf{p}` coordinates lives in :math:`V`.

Compatible discretization
_________________________

In this case


.. Local Variables:
.. mode: rst
.. End:
