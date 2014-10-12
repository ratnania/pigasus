.. role:: envvar(literal)
.. role:: command(literal)
.. role:: file(literal)
.. _basicPDE:


=================
Dive into Pigasus
=================

The *basicPDE* class
********************

The *basicPDE* class describes the following kind of partial differential equations

.. math::

   - \nabla \cdot (A \nabla u) + ( \mathbf{v} - \mathbf{w} )  \cdot \nabla u + b u &= f, ~~~~ \Omega  
     \\
     u &= g, ~~~~ \Gamma_D 
     \\
     A \nabla u \cdot \mathbf{n} &= k, ~~~~ \Gamma_N

which can be written in a Weak Formulation, for any test function :math:`\phi` that vanishes on :math:`\Gamma_D`,

.. math::
  
  \int_{\Omega} A \nabla u \cdot \nabla \phi ~d\Omega + \int_{\Omega} \left( \mathbf{v} \cdot \nabla u \right) \phi ~d\Omega + \int_{\Omega}  u \mathbf{w} \cdot \nabla \phi~d\Omega + \int_{\Omega} ( b + \nabla \cdot \mathbf{w} ) u \phi  ~d\Omega = 
  \\
  \int_{\Omega} f \phi ~d\Omega + \int_{\Gamma_N}  \left( A \nabla u + u \mathbf{w} \right) \cdot \mathbf{n}  \phi ~d\partial\Omega 

It is assumed that the domain :math:`\Omega` has a smooth boundary :math:`\partial \Omega` such that

.. math::

  \Gamma_D \cup \Gamma_N &=\partial \Omega
  \\
  \Gamma_D \cap \Gamma_N &= \emptyset

All you need to provide is 

* the geometry
* the boundary conditions
* the functions operators:
   1. The matrix function :math:`A`
   2. The vector functions :math:`\mathbf{v}` and :math:`\mathbf{w}`
   3. The scalar functions :math:`b` and :math:`f` for the right hand side [#f1]_

What you should know about *basicPDE*
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following facts are very important

* The user can access to the unknown *field* by calling::

    u = PDE.unknown

* Accessing the *right hand side* of the *basicPDE* object can be done like this::

    f = PDE.rhs  

* The *basicPDE* class has a private object *_system*. This is the linear system that descritize the PDE. The user can call the *property system* to get the *Matrix* object::
  
    M = PDE.system

This will be very useful when trying to solve more complicated problems.

A first example
^^^^^^^^^^^^^^^

   Now that you know the form of a *basicPDE* model, let's solve the following 2D example

*For given functions* :math:`A,f`, *find* :math:`u` *such that*

.. math::

  -\nabla \cdot ( A \nabla u ) = f  &  ,\Omega \\
  u = 0 &  ,\partial \Omega 

In this example we will take 

.. math::

  A = A_{\epsilon} = \begin{pmatrix}
  \epsilon \cos(\phi)^2 + \sin(\phi)^2 & (\epsilon - 1)\sin(\phi)\cos(\phi)  \\
  (\epsilon - 1)\sin(\phi)\cos(\phi)  & \epsilon \sin(\phi)^2 + \cos(\phi)^2
 \end{pmatrix}

Construction of the pde
^^^^^^^^^^^^^^^^^^^^^^^

All information needed by Pigasus can be passed through a Python dictionary namely *testcase*::

  testcase = {}
  testcase['AllDirichlet'] = True
  eps = 1.0e-6
  phi = 2 * (1./3) * pi
  c = cos(phi)
  s = sin(phi)  
  testcase['A'] = lambda x,y : [  eps * c**2 + s**2  \
                               , ( eps - 1 ) * c * s \
                               , ( eps - 1 ) * c * s \
                               , eps * s**2 + c**2]  
  testcase['f'] = lambda x,y : [1.]

.. For the example we are interested in  

.. .. literalinclude:: include/introduction/testcase.py
..     :linenos:
..     :language: python 

We are ready to create the PDE::

  PDE = basicPDE(geometry=geo, testcase=testcase)


Assembling and solving the pde
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We can assemble the PDE by calling::

  PDE.assembly()

Then solve it::

  PDE.solve()

The visualization can be done using *matplotlib* or *pylab*, by calling::

  PDE.plot() ; pl.show()
  

Atomic operators
^^^^^^^^^^^^^^^^

Each term of the last equation can be viewed as an atomic operator: *MASS, ADVECTION, transpose(ADVECTION), STIFFNESS and SECOND_DERIV*.    

* *MASS* operator implements the following contribution 

.. math:: 

  \int_{\Omega} b([u],\mathbf{x}) ~\varphi_I ~\varphi_{J} ~d\Omega


* *STIFFNESS* operator implements the following contribution 

.. math::

  \int_{\Omega} ( A([u],\mathbf{x}) ~ \nabla \varphi_I ) ~ \cdot \nabla \varphi_{J} ~d\Omega

* *ADVECTION* operator implements the following contribution 
  
.. math::

  \int_{\Omega} \varphi_{J} ~\mathbf{v}([u],\mathbf{x}) \cdot ~\nabla \varphi_{I} ~d\Omega

.. note:: The *transpose(ADVECTION)* operator is given by switching the indices of the basis functions when computing the contribution term. This writes :math:`\int_{\Omega} \varphi_{I} ~\mathbf{v}([u],\mathbf{x}) \cdot ~\nabla \varphi_{J} ~d\Omega`


* *SECOND_DERIV* operator implements the following contribution  
 
.. math::
 
   \int_{\Omega} A D^2 \varphi_I \cdot D^2 \varphi_J ~d\Omega

where in **2D**, the matrix function **A** is **3x3** and the operator D2 denotes

.. math:: 

  D^2 u = (u_{xx}, u_{xy}, u_{yy})  

.. note:: Remark that the *SECOND_DERIV* operator is defined only in the weak sens. We expect in the future a general strong form for the *basicPDE* class.


Additional features
^^^^^^^^^^^^^^^^^^^

In the sequel, we give some additional data access::

   # ... access the Discrete functional space
   PDE.space
   # ... get the shape of the associated system
   PDE.shape
   # ... access the Mass operator (returns None if not defined)
   PDE.mass
   # ... access the Stiffness operator (returns None if not defined)
   PDE.stiffness
   # ... access the advection operator (returns None if not defined)
   PDE.advection
   # ... access the advection^T operator (returns None if not defined)
   PDE.tadvection
   # ... access the second-deriv operator (returns None if not defined)
   PDE.D2
   # ... apply the PDE operator to a given field F
   Y = PDE.dot(F)
   # ... returns the unknown including the Dirichlet part
   PDE.unknown_dirichlet
   # ... assembles the function func as a rhs and includes the boundary contributions 
   PDE.update(func)
   # ... computes the norm of the field U (or the error if u is not None)
   PDE.norm(exact=u) 
   # ... interpolates the function u and sets the field U with its Splines coefficients
   PDE.interpolate(u, field=U)


.. rubric:: Footnotes

.. [#f1] Only scalar rhs are treated for the moment 

.. Local Variables:
.. mode: rst
.. End:
