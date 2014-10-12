.. role:: envvar(literal)
.. role:: command(literal)
.. role:: file(literal)
.. _gallery.rrefinement:


R-refinement
************

In this section, we show an application of the *r-refinement*. 

The *r-refinement* strategy is not new and has been used by many authors :cite:`Cao2001a`, :cite:`Cao2001b`, :cite:`Cao1999a`.

In order to do the r-refinement, an *a-posteriori estimator* can be used. 

In this example, we solve the 2D anisotropic Diffusion problem, with Homogeneous Dirichlet boundary conditions, using Pigasus.

.. math::

  - \nabla \cdot \left( A_{\epsilon} \nabla u \right) = f, ~~~~ \Omega
  \\
  u &= 0, ~~~~  \partial \Omega


The idea will be to use the Monge-Kantorovich optimization in order to construct an adaptive mesh. Then, the resolution of the previous equation will be done on the new mesh.

As described earlier, the resolution of the Monge-Kantorovich optimization problem leads to a gradiant-map. We start our example by showing how to construct such maps as a *cad_geometry* object and how to generate the corresponding metric.


Gradiant-Maps using **pigasus**
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A gradiant map is a mapping :math:`\mathbf{F}: \boldsymbol\xi \rightarrow \mathbf{x}` such that :math:`\mathbf{F}= \nabla \boldsymbol \Psi` for  a smooth function :math:`\boldsymbol\Psi`. 

In the context of IsoGeometric Analysis, the function :math:`\boldsymbol\Psi` is described by a *B-spline* or a *NURBS* curve, surface or volume. Moreover, in order to treate general geometries, we also need to add a metric, so that the gradiant definition is consistent. Using **pigasus** this can be achieved using the *cad_op_nurbs* class which relies on the *opNURBS* object. In the sequel, we show how to define and manipulate such objects. 

.. code-block:: python

   from igakit.nurbs import NURBS
   from igakit.op_nurbs import grad

   # ... creation of an arc of circle using a rational Bezier curve
   C = [[0, 1], [1, 1], [1, 0]] # 3x2 grid of 2D control points
   w = [1, np.sqrt(2)/2, 1]     # rational weigths
   U = [0,0,0, 1,1,1]           # knot vector
   crv = NURBS([U], C, weights=w)

   # ... creation of the gradiant map
   gcrv = grad(crv)

The *opNURBS* object implements similar functions to those we can find in the *NURBS* object. For instance, we can refine, elevate the spline degree of these mappings, or evaluate them at giving parametric sites in the same way as for the *NURBS* object.

.. code-block:: python

   from igakit.nurbs import NURBS
   from igakit.op_nurbs import grad
   import numpy as np

   # ... Create a random surface
   C = np.random.rand(3,3,3)
   U = [0,0,0,1,1,1]
   V = [0,0,0.5,1,1]
   nrb = NURBS([U,V], C)

   # ... create the gradiant map
   s1 = grad(nrb)

   # ... clone and elevate the spline degree
   s2 = s1.clone().elevate(0, 1).elevate(1, 1)

   # ... evaluate the two mappings
   u = v = np.linspace(0,1,100)
   xyz1 = s1.evaluate(u, v)
   xyz2 = s2.evaluate(u, v)
   print np.allclose(xyz1, xyz2, rtol=0, atol=1e-5)

A basic example is the use of gradiant map to construct the identity mapping. If we consider the following scalar function

.. math::

   f (x,y) = \frac{1}{2}\left( x^2 + y^2 \right)

we have 

.. math::

   \partial_x f (x,y) = x
   \\
   \partial_y f (x,y) = y

The construction of the gradiant map needs to define a metric, which in our case is nothing but the identity transformation, given by the **square** object from **cad_geometry**. The following script shows how to construct an interpolation of :math:`f` over the patch and then create the gradiant map

.. code-block:: python

   # ... construct a gradiant map that leads to the identity mapping
   f   = lambda x,y : 0.5 * (x**2 + y**2)
   
   from igakit.cad_geometry import square as domain
   from pigasus.interpolate.interpolation import surfint
   geo = domain(n=[31,31], p=[2,2])
   
   interpolator = surfint(geo)
   y   = interpolator.interpolate(f)

   # ... define the metric
   met = geo[0]
   C = np.zeros_like(met.points)
   C[:,:,:2] = met.points[:,:,:2]
   C[:,:,2]  = y
   # ... construct the gradiant map
   nrb  = NURBS(met.knots, C, weights=met.weights)
   c  = grad(nrb)
   
   u = np.linspace(0.,1.,20)
   v = np.linspace(0.,1.,30)

   Q = c.evaluate(u, v)
   x   = Q[:,:,0]
   y   = Q[:,:,1]

   X,Y = np.meshgrid(u,v)
   
   print np.allclose(x, X.transpose(), rtol=0, atol=1e-12)
   print np.allclose(y, Y.transpose(), rtol=0, atol=1e-12)
   
   # ... derivatives evaluation. expected values: xdu = 1, xdv = 0, ydu = 0, ydv = 1
   Q = c.evaluate_deriv(u, v, nderiv=1)
   x   = Q[0,:,:,0]
   y   = Q[0,:,:,1]
   xdu = Q[1,:,:,0]
   ydu = Q[1,:,:,1]
   xdv = Q[2,:,:,0]
   ydv = Q[2,:,:,1]
   
   print np.allclose(xdu, np.ones_like(xdu), rtol=0, atol=1e-10)
   print np.allclose(xdv, np.zeros_like(xdv), rtol=0, atol=1e-10)
   print np.allclose(ydu, np.zeros_like(ydu), rtol=0, atol=1e-10)
   print np.allclose(ydv, np.ones_like(ydv), rtol=0, atol=1e-10)


.. For the moment, **pigasus** can not handle these kind of mappings automatically. By this, we mean that we can not create a *basicPDE* operator while giving a *opNURBS* object as a geometry. For the moment, we need to use the *metric* object. 
.. 
.. .. code-block:: python
.. 
..    from igakit.nurbs import NURBS
..    from igakit.op_nurbs import grad
..    import numpy as np
.. 
..    # ... Create a random surface
..    C = np.random.rand(3,3,3)
..    U = [0,0,0,1,1,1]
..    V = [0,0,0.5,1,1]
..    nrb = NURBS([U,V], C)
.. 
..    # ... create the gradiant map
..    s1 = grad(nrb)
.. 
..    # ... create a cad_op_nurbs object for additional information about the geometry 
..    from igakit.cad_geometry import cad_op_nurbs
..    cad_s1 = cad_op_nurbs(s1)
.. 
..    # ... create a cad_geometry object that will contain the gradiant map
..    from igakit.cad_geometry import cad_geometry
..    geo_m = cad_geometry() 
..    geo_m.append(cad_s1)
.. 
..    # ... create a Metric object
..    from pigasus.fem.metric import metric
..    Metric = metric(geometry=geo_m)


The **adaptiveMesh** module 
^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Pigasus** offers an implemented module and class that generates adaptive meshes. For the moment, the mesh generation process is based on the resolution of the Monge-Kantorovitch optimization problem. In the future, other methods will be implemented.

In the next example, we show how to generate a *cad_geometry* object using the **adaptiveMeshMA** class, which can be found in the **plugin** directory.

.. code-block:: python

   adapt = adaptiveMeshMA(geo_h, geo_H=geo_H, verbose=verbose, bc_neumann=bc_neumann)
   geo_f = adapt.construct(rho0, rho1)
   geo_f.save("adaptive_mesh.xml")

In order to construct an initial guess for the Picard solver, we may use a two grids method. In this case, you need to specify a **coarse** geometry for the argument *geo_H*.   

.. note:: The interior knots of *geo_H* must be included in *geo_h*.

The reader can check that the generated geometry inherits all data information from **geo_h**. Which means that both **geo_h** and **geo_f** have the same internal/external faces and connectivity. Moreover, the generated patch is stored in the **XML** file with the **gradiant** attribut.


A-posteriori estimates for Elliptic PDE
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We define the bilinear form 

.. math::

   B(\psi, \phi) = \int_{\Omega}  A \nabla \psi \cdot \nabla \phi,~~~\forall \phi \in H^1_{0}(\Omega)

Let :math:`e_h = \psi - \psi_h` be the error of the IsoGeometric approximation :math:`\psi_h`. We have

.. math::

   B(e_h,\phi) = B(\psi, \phi) - B(\psi_h, \phi) = (f,\phi) - B(\psi_h, \phi) = \sum_{Q \in \mathcal{Q}_h} \left(\int_{Q} f \phi - \int_{Q} A \nabla \psi_h \cdot \nabla \phi \right)

Moreover, we get, after integrating by parts, 

.. math::

   B(e_h,\phi) = \sum_{Q \in \mathcal{Q}_h} \left(\int_{Q} (f + \nabla^2 \psi_h )\phi - \int_{\partial Q} ( A \nabla \psi_h \cdot \mathbf{n}) \phi \right)

by rearranging terms, we get

.. math::

   B(e_h,\phi) = \sum_{Q \in \mathcal{Q}_h} \int_{Q} R_Q \phi + \sum_{\gamma \in \mathcal{I}_h} \int_{\gamma} J_{\gamma} \phi 

where :math:`\mathcal{I}_h` is the set of interior *curved* edges, 

.. math::

   J_{\gamma} = \mathbf{n}^- \cdot \nabla \psi_h^- + \mathbf{n}^+ \cdot \nabla \psi_h^+

for any interior *curved* edge :math:`\gamma`, while

.. math::
  :label: aposteriori_est
   
   R_Q = f + \nabla^2 \psi_h,~~~\mbox{in} ~Q

Due to the continuity of :math:`\nabla \psi_h`, we have :math:`J_{\gamma}=0`. Therefor,

.. math::

   B(e_h,\phi) = \sum_{Q \in \mathcal{Q}_h} \int_{Q} R_Q \phi

Finaly, we get the classical estimation

.. math::
   
   \| e_h \|^2 \lesssim \sum_{Q \in \mathcal{Q}_h} h_Q^2 \| R_Q \|_{L^2(Q)}^2


.. Local Variables:
.. mode: rst
.. End:
