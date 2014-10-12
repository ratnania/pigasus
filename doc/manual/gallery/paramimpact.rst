.. role:: envvar(literal)
.. role:: command(literal)
.. role:: file(literal)
.. _gallery.paramimpact:



Parametrization Impact
**********************

In this section, we will show the impact of the Mesh parametrization. There are mainly 2 kinds of meshes that are of interest (fig :ref:`MESH1` and :ref:`MESH2`).

Example of the **MESH1** type on the unit disk

.. image:: include/paramimpact/meshes_1.png
   :width: 8cm
   :height: 8cm
..    :label: MESH1

Example of the **MESH2** type on the unit disk

.. image:: include/paramimpact/meshes_2.png
   :width: 8cm
   :height: 8cm
..    :label: MESH2

Let us consider the following problem:

.. math::
   :label: pb_A
   
   \left\{
   \begin{aligned}
   \partial_t u &= [u,\psi], &~~ \forall (t,\mathbf{x}) \in [0,T] \times \Omega
   \\
   u (t,\cdot) &= 0, &~~ \forall (t,\mathbf{x}) \in [0,T] \times \partial \Omega
   \\
   u (t=0,\cdot) &= \psi(\cdot), &~~ \forall \mathbf{x} \in \Omega
   \end{aligned} 
   \right.

where :math:`\psi (x,y) = \sin ( 2 \pi (x^2 + y^2) )`.

Because of the boundary condition, we have 

.. math::

   \int_{\Omega} [u,\psi] \phi = - \int_{\Omega} [\phi,\psi] u

Therefor, we can proove the conservation of the following quantities:

Some basic properties
^^^^^^^^^^^^^^^^^^^^^

Energy Conservation
___________________

.. math::

  \frac{d}{dt} \int_{\Omega} u^2 = 0
  \\
  \frac{d}{dt} \int_{\Omega} u \psi = 0

TODO

which leads to the following variational formulation of (Eq. :ref:`pb_A`),

.. math::
   :label: eq_u

    \frac{d}{dt} \int_{\Omega} u \phi = \frac{1}{2} \int_{\Omega} \nabla \times \psi (\phi \nabla u - u \nabla \phi) 

therefor, we get the linear system,

.. math::

   M \dot{[u]} = \frac{1}{2} (A-A^T) [u]

where :math:`M_{b,b^{\prime}} = \int_{\Omega} \phi_b \phi_{b^{\prime}}` and 
:math:`A_{b,b^{\prime}} = \int_{\Omega} (\nabla \times \psi \cdot \nabla \phi_b ) \phi_{b^{\prime}}` are the Mass and Advection matrices.

It is easy to see (from Eq. :ref:`eq_u`), that the :math:`L^2` norm is conserved. We use Crank-Nicolson time scheme and get,

.. math::

   \left( M - \frac{dt}{2} A \right)[u]^{n+1} = \left( M + \frac{dt}{2} A \right) [u]^{n}


Discrete Energy Conservation
____________________________

We have

* :math:`\int_{\Omega} |u^{n+1}|^2 = \int_{\Omega} |u^{n}|^2`
* :math:`\int_{\Omega} u^{n+1} \psi = \int_{\Omega} u^{n} \psi`

Let's go back to the equation 

.. math::
  
   u^{n+1} = u^{n} + \frac{1}{2}[u^{n+1},\psi] + \frac{1}{2}[u^{n},\psi]

we multiply it by :math:`u^{n+1}` and again by :math:`u^{n}` and then integrate to get

.. math::
   :label: eq_1

   \int_{\Omega}|u^{n+1}|^2 = \int_{\Omega} u^{n} u^{n+1} + \frac{1}{2}\int_{\Omega} [u^{n+1},\psi]u^{n+1} + \frac{1}{2}\int_{\Omega} [u^{n},\psi]u^{n+1}

.. math::
   :label: eq_2

   \int_{\Omega} u^{n} u^{n+1} = \int_{\Omega}|u^{n}|^2 + \frac{1}{2}\int_{\Omega} [u^{n+1},\psi]u^{n} + \frac{1}{2}\int_{\Omega} [u^{n},\psi]u^{n} 

but, we have 

.. math::

   \int_{\Omega} [u^{n+1},\psi]u^{n+1} = \int_{\Omega} [u^{n},\psi]u^{n} = 0

and 

.. math::

   \int_{\Omega} [u^{n+1},\psi]u^{n+1} =- \int_{\Omega} [u^{n+1},\psi]u^{n}

so that by summing (Eq :ref:`eq_1` and :ref:`eq_2`), we get the desired relation.


Numerical Results
_________________

Simulations were done on a :math:`16 \times 16` grid, with quadratic *NURBS*
and :math:`dt = 0.5`.

Next we plot the evolution (Long times) of the :math:`L^2` norm for **MESH1**.

.. image:: include/paramimpact/MESH1/energy.png
   :width: 9cm
   :height: 9cm

Next we plot the evolution (Long times) of the :math:`L^2` norm for **MESH2**.

.. image:: include/paramimpact/MESH2/energy.png
   :width: 9cm
   :height: 9cm

Next we plot the evolution (Long times) of the helicity norm for **MESH1**.

.. image:: include/paramimpact/MESH1/energyH.png
   :width: 9cm
   :height: 9cm

Next we plot the evolution (Long times) of the helicity norm for **MESH2**.

.. image:: include/paramimpact/MESH2/energyH.png
   :width: 9cm
   :height: 9cm

Notice that we get a linear dependece on time, in the case of **MESH2**. A linear regression gives the slops :math:`6.01 10^{-15}` and 
:math:`3.0 10^{-15}` for the :math:`L^2` norm and the helicity quantity.


Impact of the Jacobian singularities
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In the next figures, we show the impact of the singularities of the Jacobian on the evolution of a pulse. As expected, the solution will blow up when the pulse pass through the *polar-like* singularity. 

Impact of the jacobian singularities for **MESH1**


.. image:: include/paramimpact/MESH1/pulse/u_0.png
   :width: 8cm
   :height: 8cm
..    :label: impact_sing_mesh1a

.. image:: include/paramimpact/MESH1/pulse/u_1.png
   :width: 8cm
   :height: 8cm
..    :label: impact_sing_mesh1b

.. image:: include/paramimpact/MESH1/pulse/u_2.png
   :width: 8cm
   :height: 8cm
..    :label: impact_sing_mesh1c

.. image:: include/paramimpact/MESH1/pulse/u_3.png
   :width: 8cm
   :height: 8cm
..    :label: impact_sing_mesh1d

.. image:: include/paramimpact/MESH1/pulse/u_4.png
   :width: 8cm
   :height: 8cm
..    :label: impact_sing_mesh1e

.. image:: include/paramimpact/MESH1/pulse/u_5.png
   :width: 8cm
   :height: 8cm
..    :label: impact_sing_mesh1f


Impact of the jacobian singularities for **MESH2**


.. image:: include/paramimpact/MESH2/pulse/u_0.png
   :width: 8cm
   :height: 8cm
..    :label: impact_sing_mesh2a

.. image:: include/paramimpact/MESH2/pulse/u_1.png
   :width: 8cm
   :height: 8cm
..    :label: impact_sing_mesh2b

.. image:: include/paramimpact/MESH2/pulse/u_2.png
   :width: 8cm
   :height: 8cm
..    :label: impact_sing_mesh2c

.. image:: include/paramimpact/MESH2/pulse/u_3.png
   :width: 8cm
   :height: 8cm
..    :label: impact_sing_mesh2d

.. image:: include/paramimpact/MESH2/pulse/u_4.png
   :width: 8cm
   :height: 8cm
..    :label: impact_sing_mesh2e

.. image:: include/paramimpact/MESH2/pulse/u_5.png
   :width: 8cm
   :height: 8cm
..    :label: impact_sing_mesh2f


.. Local Variables:
.. mode: rst
.. End:
