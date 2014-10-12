.. role:: envvar(literal)
.. role:: command(literal)
.. role:: file(literal)
.. _gallery.mhdequilibrium:


MHD Equilibrium
***************

In MHD, the helical symmetry implies that any physical quantity will depend only on :math:`r` and :math:`u=l \phi + k z` ( more details can be found in :cite:`biskamp_book`). The general solution of Gauss's law for magnetism :math:`\divs \cdot \mathbf{B} = 0` in helical symmetry can be written in the form:

.. math::

  \mathbf{B} = \mathbf{h} \times \nabla \psi(r,u) + \mathbf{h} f(r,u).

where:

* :math:`\mathbf{h}` is defined by:

.. math::  

  \mathbf{h} = \frac{r}{l^2+k^2r^2} \nabla r \times \nabla u,

is tangent to the helix :math:`r=\mbox{const}`, :math:`u=\mbox{const}`.

* :math:`\psi` is the helical flux function. :math:`\psi` is essantially the component of the vector potentional in the direction of the ignorable coordinate :cite:`biskamp_book`,
  
.. math::

  \psi = - \frac{\mathbf{A} \cdot \mathbf{h} }{\| \mathbf{h} \|^2} = - (l A_z - krA_{\phi})

* :math:`f` is the helical magnetic field. We have:

.. math::

  f = \frac{\mathbf{B} \cdot \mathbf{h} }{\| \mathbf{h} \|^2} = (l B_z - krB_{\phi})


The starting point for the magnetostatic equilibrium (using primitive variables) is the following balance equations (\textit{c.f} MHD's equations in :cite:`jardin_book`):

.. math::

  \rho ( \partial_t \mathbf{v} + \mathbf{v} \cdot \nabla \mathbf{v} ) = - \nabla p + \mathbf{J} \times \mathbf{B} - \rho \nabla \phi_g

where :math:`p` is the pressure and :math:`\phi_g` is the gravitational potential. Usually, in laboratory, this term is negligible. However, it is very important in astrophysical systems.

In the case of a steady state, :math:`\partial_t \cdot = 0`, we end up with,

.. math::

  \rho \mathbf{v} \cdot \nabla \mathbf{v}  = - \nabla p + \mathbf{J} \times \mathbf{B} 

which traduces the balance of forces inside the electrical fluid.

Equilibrium in the absence of toroidal flux
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When we can neglect the term :math:`\rho \mathbf{v} \cdot \nabla \mathbf{v}`, the equilibrium writes simply,

.. math::

  \nabla p = \mathbf{J} \times \mathbf{B}

The general form of the MHD equilibrium, in this case, is 

.. math::

  \mathcal{L}_{k,l} \psi = \frac{2kl}{l^2+k^2r^2}f - f f^{\prime} - (l^2+k^2r^2)p^{\prime}

where the differential operator :math:`\mathcal{L}_{k,l}` is :

.. math::

  \mathcal{L}_{k,l} \psi = \frac{l^2+k^2r^2}{r}\{ \partial_r \frac{r}{l^2+k^2r^2}\partial_r \psi + \frac{1}{r} \partial_{uu} \psi \}


Plane case
^^^^^^^^^^


In this case, we have :math:`l=1` and :math:`k=0`, so that :math:`\mathcal{L}_{0,1}=\Delta`. This leads to the equation:

.. math::

  \nabla^{2} \psi = r^{-1} \partial_r (r\partial_r \psi) + r^{-2}\partial_{\phi \phi} \psi = -ff^{\prime} - p^{\prime}


Axisymmetric case
_________________

In this case, we have :math:`l=0` and :math:`k=1`, so that :math:`\mathcal{L}_{1,0}=\Delta^{\star}`. It is one of the most important case in tokamaks. This leads to the Grad-Shafranov equation.

.. math::

  \Delta^{\star} \psi = r \partial_r (r^{-1} \partial_r \psi) + \partial_{zz} \psi = -ff^{\prime} - r^2p^{\prime}


Soloviev equilibrium
____________________

Soloviev equilibrium is a special case of the Grad-Shafranov equation. This happens when the right hand side is independent of :math:`\psi`. Thus:

.. math::

  p = - | p^{\prime} | \psi + p_0,~~~~\mbox{and}~~~~f^2 = f_0^2 + \frac{2 \gamma | p^{\prime}|}{1+\alpha^2}\psi

The quantities :math:`p_0` and :math:`\frac{f_0}{R_0}` are the pressure and the toroidal field, on the magnetic axis :math:`r=R_0,~z=0,~ \psi=0`.

In this configuration, the Grad-Shafranov equation writes:

.. math::

  \Delta^{\star} \psi = | p^{\prime} |  (r^2-\frac{\gamma}{1+\eta^2})

a solution of such equation is :

.. math::

  \psi = \frac{| p^{\prime}|}{2(1+ \eta ^2)}((r^2-\gamma)z^2+\frac{\eta ^2}{4}(r^2-R_{0}^2)^2)

In the sequel, we shall use the cartesian coordinates :math:`(x,y)`:

.. math::

  r = R_0 + a x = R_0 (1+\epsilon x),~~~~\mbox{and}~~~~z = ay

where :math:`\epsilon = \frac{a}{R_0}`.

The Soloviev equilibrium writes :

.. math::

  - \nabla \cdot \left(  \frac{\nabla \psi}{1+\epsilon x} \right) = a^2 \alpha R_0^2 (1+\epsilon x) + \frac{a^2 \beta}{1+\epsilon x}

this can be solved under homogeneous Dirichlet boundary condition, as there exists a level surface where the solution vanishes. We have introduced the quantities:

.. math::

  \beta = -\frac{\lambda}{b^2 \epsilon},~~~~\alpha = \frac{4(a^2+b^2)\epsilon+a^2(2\lambda-\epsilon^3)}{2R_0^2\epsilon a^2 b^2}

Notice that the points :math:`(\pm1,0)` and :math:`(0,\pm \frac{b}{a})` are on the level surface :math:`\psi=0`.


Equilibrium with toroidal flux
______________________________

In this case we keep the term :math:`\rho \mathbf{v} \cdot \nabla \mathbf{v}`, we have,

.. math::

  \rho \mathbf{v} \cdot \nabla \mathbf{v} + \nabla p = \mathbf{J} \times \mathbf{B}

we get a generalized Grad-Shafranov equation:

.. math::

  \Delta^{\star} \psi = r \partial_r (r^{-1} \partial_r \psi) + \partial_{zz} \psi = -ff^{\prime} - r^2 \partial_{\psi} p

:math:`f` is still a function of the only variable :math:`\psi`, but now the pressure :math:`p` is a function of :math:`r` and :math:`\psi`.


Nonlinear equilibrium
_____________________

In this case, we look for the solution of the Grad-Shafranov nonlinear equation :

.. math::

  \Delta^{\star} \psi = F(r,\psi)

where, :math:`F` is a nonlinear function of :math:`\psi`:

.. math::

  F(r,\psi) := -(r^2 f_1(\psi) + f_2(\psi)) 

This non linear partial differential equation can be solved using Picard or Newton methods.


Grad-Shafranov as a nonlinear Poisson equation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Let us consider the solution of the Grad-Shafranov equation 

.. math::

	\Delta^{\star} \Psi = \Psi_{RR} - \frac{1}{R} \Psi_{R} + \Psi_{ZZ} = F(\Psi, R, Z)
	\label{eq:grad_shafranov}

as the solution of the following nonlinear Poisson problem

.. math::

	\nabla^{2} \Psi = \Psi_{RR} + \Psi_{ZZ} = F(\Psi, R, Z) + \frac{1}{R} \Psi_{R}
	\label{eq:grad_shafranov2}

In order to avoid the variations due to the :math:`\Psi_R` term, we use the following change of variables

.. math::

	U(R,Z) = \frac{\Psi(R,Z)}{\sqrt{R}}
	\label{eq:GS_change_var}

Now the Grad-Shafranov writes

.. math::

	\nabla^2 U = \mathcal{F}(U, R,Z),  \Omega
	\\
	U = \mbox{cte},  \partial \Omega
	\label{eq:grad_shafranov3}

where

.. math::

	\mathcal{F}(U,R,Z) = \frac{1}{\sqrt{R}} F(\sqrt{R} U, R, Z) + \frac{3}{4R^2} U

The following plot shows the analytical solution for **ITER** relevant parameters

.. image:: include/mhd_equilibrium/soloviev_iter.png
   :width: 9cm
   :height: 10cm

In the sequel we show how to construct the corresponding polar mesh.

.. image:: include/mhd_equilibrium/ex5.png
   :width: 12cm
   :height: 14cm

The complete script can be found here :download:`script <include/mhd_equilibrium/ex5.py>`

We Have enforced an interpolation constraint the extremeties of the curve. This is done by calling the function *MakeConstraint* as the following::

   # ... Define the first point A on the face = 0
   A = [0.,0.] ; face = 0
   constraint = MakeConstraint("C0", face, A)
   constraints.append(constraint)

.. note:: As we can expect, the approximation is not good where the number of input data is not suffisant.

Using a closed curve, we can generate the corresponding *2D* description, which leads to the following mesh.

.. image:: include/mhd_equilibrium/ex6.png
   :width: 14cm
   :height: 14cm

The *cad_geometry* module contains a function for this purpose::

  geo = geo_f.polarExtrude()

As a first validation, we test our solver on the following problem

.. math::

	\Delta^{\star} \Psi = R^2

which we solve both numerically and analytically. The analytical solution is given by

.. math::

   Psi(R,Z) = R**4/8 + d1 + d2 * R**2 + d3 * (R**4 - 4 * (R*Z)**2)

with

.. math::

   d1, d2, d3 = [ 0.07538503, -0.20629496, -0.03143371]

The numerical solution is given in the following figure

.. image:: include/mhd_equilibrium/soloviev_iter.png
   :width: 18cm
   :height: 14cm
  


.. Local Variables:
.. mode: rst
.. End:
