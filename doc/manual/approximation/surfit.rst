.. role:: envvar(literal)
.. role:: command(literal)
.. role:: file(literal)
.. _approximation.surfit:

Surface fitting
***************


The following example shows the approximation of a function on the unit square

.. code-block:: python

   import matplotlib.pyplot    as plt
   import numpy                as np
   from igakit.cad_geometry import square as patch
   from pigasus.fit.surfit import surfit

   sin = np.sin ; pi = np.pi

   #-----------------------------------
   # ...
   nx = 15 ; ny = 15
   px = 3 ; py = 3
   # ...

   # ...
   u = np.linspace(0.,1.,100)
   v = np.linspace(0.,1.,200)
   n = len(u)*len(v)

   U,V = np.meshgrid(u,v)

   u = U.reshape(n)
   v = V.reshape(n)

   f = sin(2*pi*u) * sin(3*pi*v)
   # ...
   #-----------------------------------

   # ...
   geo = patch(n=[nx,ny], p=[px,py])
   fit = surfit(geometry=geo, constraints=[], alpha=1., rational=0)
   list_nrb = fit.construct([f], uvk=[[u], [v]], exportGeometry=False)
   # ...

   # ...
   nrb = list_nrb[0]

   u = np.linspace(0.,1.,200)
   v = np.linspace(0.,1.,200)
   U,V = np.meshgrid(u,v)
   U = U.transpose()
   V = V.transpose()

   P = nrb.evaluate(u,v)
   f = P[:,:,0]

   plt.contourf(U,V,f)
   plt.colorbar(); plt.show()
   # ...

which leads to the following plot

.. image:: include/approximation/surfit_ex1.png
   :width: 12cm
   :height: 9cm


Mesh generation using **surfit**
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Another application of **surfit** is the generation of aligned meshes. In the following example, we construct an aligned adaptive mesh for the anisotropic diffusion problem

.. math::
   :label: anidiff_eq
   
   - \nabla \cdot (A_{\epsilon} \nabla u_{\epsilon}) = f

on the unit square, with Dirichlet, Neumann or periodic boundary conditions. 

The matrix :math:`A_{\epsilon}` is of the form:

.. math::

    A_{\epsilon} = \epsilon \mathbb{I} + (1-\epsilon) \mathbf{b} \otimes \mathbf{b}

where :math:`\mathbf{b} = \frac{\mathbf{B}}{\| \mathbf{B} \|}` and :math:`\epsilon << 1` is the characteristic parameter for the anisotrpy problem.

We consider here the case where the solution :math:`u` is of the form:

.. math::

   u(x,y) = \epsilon \operatorname{sin}\left(\pi n y\right) \operatorname{cos}\left(\pi m x\right) + \operatorname{sin}\left(\alpha \left(y^{2} - y\right) \operatorname{cos}\left(\pi m x\right) + \pi n y\right)

.. image:: include/approximation/ex1_mapping2D_potential.png
  :width: 9cm
  :height: 9cm

which corresponds to the following magnetic field 

.. math::   

   \mathbf{B}(x,y) = \begin{pmatrix}
   - \pi \epsilon n \operatorname{cos}\left(\pi m x\right) \operatorname{cos}\left(\pi n y\right) - \left(\alpha \left(2 y -1\right) \operatorname{cos}\left(\pi m x\right) + \pi n\right) \operatorname{cos}\left(\alpha \left(y^{2} - y\right) \operatorname{cos}\left(\pi m x\right) + \pi n y\right)
   \\
   - \pi \alpha m \left(y^{2} - y\right) \operatorname{sin}\left(\pi m x\right) \operatorname{cos}\left(\alpha \left(y^{2} - y\right) \operatorname{cos}\left(\pi m x\right) + \pi n y\right) - \pi \epsilon m \operatorname{sin}\left(\pi m x\right) \operatorname{sin}\left(\pi n y\right)
   \end{pmatrix}

The Neumann boundary is :math:`\Gamma_N = \{ x=0 \} \cup \{ x=1 \}`, while the Dirichlet's one is :math:`\Gamma_D = \{ y=0 \} \cup \{ y=1 \}`.

The resulting mesh is given in the next plot

.. image:: include/approximation/ex1_mapping2D_mesh.png
  :width: 9cm
  :height: 9cm

The following script describes the different steps that lead to this construction

.. code-block:: python

   import matplotlib.pyplot    as plt
   import numpy                as np
   from igakit.cad_geometry import square as patch
   from pigasus.fit.surfit import surfit, compute_uk

   sin = np.sin ; cos = np.cos ; pi = np.pi

   #-----------------------------------
   # ...
   #method = "uniform"
   method = "chord"
   #method = "centripetal"
   # ...

   # ...
   nx = 31 ; ny = 31
   px = 3 ; py = 3
   # ...

   # ...
   eps = 1.e-3
   alpha = 2.
   m = 1
   n = 1
   psi   = lambda x,y: eps*sin(n*pi*y)*cos(m*pi*x) \
   + sin(alpha*(y**2 - y)*cos(m*pi*x) + n*pi*y)
   psidx = lambda x,y:  -pi*alpha*m*(y**2 - y)*sin(pi*m*x) \
   * cos(alpha*(y**2 - y)*cos(pi*m*x) + pi*n*y) \
   - pi*eps*m*sin(pi*m*x)*sin(pi*n*y)
   psidy = lambda x,y: pi*eps*n*cos(pi*m*x)*cos(pi*n*y) \
   + (alpha*(2*y - 1)*cos(pi*m*x) + pi*n) \
   * cos(alpha*(y**2 - y)*cos(pi*m*x) + pi*n*y)
   # ...

   # ...
   from pigasus.utils.impeqpy import impeqpy
   imp=impeqpy()
   xk = [] ; yk = []
   uk = [] ; vk = []
   xgrid  = list(np.linspace(0.,1.,100))
   ygrid = np.zeros_like(xgrid)

   list_level = np.linspace(0,0.99,50)
   for level in list_level:
       y0 = ygrid
       imp.solve2Dx(psi,psidy,level,xgrid,ygrid,y0=y0)

       levela = level / 2
       xk += list(xgrid)
       yk += list(ygrid)
       list_Q = zip(xgrid, ygrid)
       _uk = compute_uk(list_Q, method=method)
       uk      += list(_uk)
       vk      += list(levela * np.ones_like(_uk))
   #    plt.plot(xgrid, ygrid, '-b')

       levelb = (2-level)/2
       zgrid = 1.-ygrid[::-1]
       xk += list(xgrid)
       yk += list(zgrid)
       list_Q = zip(xgrid, zgrid)
       _uk = compute_uk(list_Q, method=method)
       uk      += list(_uk)
       vk      += list(levelb * np.ones_like(_uk))
   #    plt.plot(xgrid, zgrid, '-r')

   #plt.show()
   #-----------------------------------

   # ...
   geo = patch(n=[nx,ny], p=[px,py])
   fit = surfit(geometry=geo, constraints=[], alpha=1., rational=0)
   # ...

   # ...
   geo_f = fit.construct([xk, yk], uvk=[[uk], [vk]], exportGeometry=True)
   # ...

   # ...
   srf = geo_f[0]
   D = srf.evalMesh(10)
   for d in D:
       plt.plot(d[:,0], d[:,1], '-k')
   plt.show()
   # ...

.. note:: The numerical solutions of the following tests are given :download:`here <include/approximation/anidiff_tests_potential.py>`   

Test 101 
________

We consider here the case where the solution :math:`u` is of the form:

.. math::
   
   u(x,y) = \epsilon \operatorname{sin}\left(\pi n y\right) \operatorname{cos}\left(\pi m x\right) + \operatorname{sin}\left(\pi n y\right)

which corresponds to the following magnetic field 

.. math::   

   \mathbf{B}(x,y) = \begin{pmatrix}
   1
   \\
   0
   \end{pmatrix}

.. todo:: rajouter le script et les figures

The analytical solution for the test case **101**: (left) :math:`\epsilon=10^{-1}`, (middle) :math:`\epsilon=10^{-4}` and (right) the limit case.

.. image:: include/approximation/test101/u_eps_1e-1.png
  :width: 9cm
  :height: 9cm

.. image:: include/approximation/test101/u_eps_1e-4.png
  :width: 9cm
  :height: 9cm

.. image:: include/approximation/test101/u_limit.png
  :width: 9cm
  :height: 9cm

Test 102 
________

We consider here the case where the solution :math:`u` is of the form:

.. math::

   u(x,y) = \epsilon \operatorname{sin}\left(\pi n y\right) \operatorname{cos}\left(\pi m x\right) + \operatorname{sin}\left(4 \pi n y \left(- y + 1\right)\right)

which corresponds to the following magnetic field 

.. math::   

   \mathbf{B}(x,y) = \begin{pmatrix}
   1
   \\
   0
   \end{pmatrix}

.. todo:: rajouter le script et les figures

The analytical solution for the test case **102**: (left) :math:`\epsilon=10^{-1}`, (middle) :math:`\epsilon=10^{-4}` and (right) the limit case.

.. image:: include/approximation/test102/u.png
  :width: 9cm
  :height: 9cm

.. image:: include/approximation/test102/u_eps_1e-4.png
  :width: 9cm
  :height: 9cm

.. image:: include/approximation/test102/u_limit.png
  :width: 9cm
  :height: 9cm

Test 110 
________

We consider here the case where the solution :math:`u` is of the form:

.. math::

   u(x,y) = \epsilon \operatorname{sin}\left(\pi n y\right) \operatorname{cos}\left(\pi m x\right) + \operatorname{sin}\left(\alpha \left(y^{2} - y\right) \operatorname{cos}\left(\pi m x\right) + \pi n y\right)

which corresponds to the following magnetic field 

.. math::   

   \mathbf{B}(x,y) = \begin{pmatrix}
   - \pi \epsilon n \operatorname{cos}\left(\pi m x\right) \operatorname{cos}\left(\pi n y\right) - \left(\alpha \left(2 y -1\right) \operatorname{cos}\left(\pi m x\right) + \pi n\right) \operatorname{cos}\left(\alpha \left(y^{2} - y\right) \operatorname{cos}\left(\pi m x\right) + \pi n y\right)
   \\
   - \pi \alpha m \left(y^{2} - y\right) \operatorname{sin}\left(\pi m x\right) \operatorname{cos}\left(\alpha \left(y^{2} - y\right) \operatorname{cos}\left(\pi m x\right) + \pi n y\right) - \pi \epsilon m \operatorname{sin}\left(\pi m x\right) \operatorname{sin}\left(\pi n y\right)
   \end{pmatrix}

 
.. todo:: rajouter le script et les figures

The analytical solution for the test case **110**: (left) :math:`\epsilon=10^{-1}`, (middle) :math:`\epsilon=10^{-4}` and (right) the limit case.

.. image:: include/approximation/test110/u_eps_1e-1.png
  :width: 9cm
  :height: 9cm

.. image:: include/approximation/test110/u_eps_1e-4.png
  :width: 9cm
  :height: 9cm

.. image:: include/approximation/test110/u_limit.png
  :width: 9cm
  :height: 9cm

.. Local Variables:
.. mode: rst
.. End:
