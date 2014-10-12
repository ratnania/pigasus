.. role:: envvar(literal)
.. role:: command(literal)
.. role:: file(literal)
.. _gallery.nonlinpde:


Non-linear Partial Differential Equations
*****************************************

In order to keep things general, *Pigasus* offers a specific class for Non-linear PDEs. Using again the *basicPDE* class, the functions :math:`A, b, v, w, D2` can depend on any unknown variable of your problem, in addition to the spatial coordinates. Let us, for instance, suppose that the differential operators and right hand terms :math:`\{ \left(\mathcal{L}_i,f_i \right), ~ i=1,\cdots,r \}` defines your problem, modulo a specific time scheme. The operators :math:`\left(\mathcal{L}_i \right)_{i=1,\cdots,r}` can be described by the *basicPDE* class even in a non-linear case. In fact, the functions :math:`\mathcal{L}_i.A, \mathcal{L}_i.b, \mathcal{L}_i.v, \mathcal{L}_i.w, \mathcal{L}_i.D2` can use any global variable, which makes thing easier.

A basic example
^^^^^^^^^^^^^^^

In the sequel, we shall consider the following problem:

*For given functions F,g and k , find* :math:`u` *such that:*

.. math::
  :label: nonlin-eq

  \mathcal{L} u &= F ( [u], \mathbf{x} ),~~~~ \Omega 
  \\
  u &= g, ~~~~ \Gamma_D 
  \\
  \left( \mathcal{L}.A \right) \nabla u \cdot \mathbf{n} &= k, ~~~~ \Gamma_N  

After relaxing the equation, we get a classical Matrix system to solve:

.. math::
  :label: nonlin-linsys

  \mathcal{L} [u] = \mathcal{F}([u])

Now it is time to use your favorite Non-linear solver.


.. todo::

   Show examples using scipy and Petsc4py

As a first example, let us consider the following non-linear Poisson equation

.. math::

  - \Delta u = - a \mathrm{e}^{\beta u}

which occurs in combustion theory. It also models the electrostatic potential in a charged body.

The general form of solutions is :

.. math::

  u(x,y) = \frac{1 }{ \beta } \ln { \frac{8 C} { a \beta } } - \frac{2 }{ \beta } \ln { | (x+A)^2 + (y+B)^2 - C | }

for more solutions, we refer to :cite:`polyanin_book`.

In order to have the function :math:`u`, vanishing at the boundary, we shall take the following values of parameters:
:math:`C = - \frac{1}{2},~~~~A = B = 0,~~~~a \beta = - 4`
which gives,

.. math::

  u(x,y) =  - \frac{2 }{ \beta } \ln { | x^2 + y^2 + \frac{1}{2} | }

One can easily check that :math:`u` verifies:

.. math::
  :label: nonlin-ex1

  - \Delta u = \frac{4}{\beta} \mathrm{e}^{\beta u}

In the following test, we took :math:`\beta = -1`.

To have homogeneous Dirichlet boundary condition, the domain will be a circle of radius :math:`\frac{\sqrt{2}}{2}`, centered at :math:`0`.


Picard algorithm
^^^^^^^^^^^^^^^^


To solve iteratively :eq:`nonlin-linsys`, let us start with the Picard algorithm, which is the simplest one but also the less accurate. The Picard algorithm can be written as

* :math:`X^0` is given,

* knowing :math:`X^n`, we solve : 

.. math::

  \mathcal{S} X^{n+1} = \mathcal{F}(X^n)

where :math:`\mathcal{S}` is the *Stiffness* matrix.

You can use the *poisson_picard* class from the *poisson_nonlin* module.

.. literalinclude:: include/demo/test_nonlin_ex1_picard.py
    :linenos:
    :language: python
    :lines: 42-88

See :download:`the complete script <include/demo/test_nonlin_ex1_picard.py>`.

Which gives the numerical solution

.. image:: include/demo/test_nonlin_ex1_picard.png
   :width: 15cm
   :height: 9cm

Newton algorithm
^^^^^^^^^^^^^^^^

The Residual function is defined as:

.. math::

  \mathcal{G}(X) = \mathcal{S} X - \mathcal{F} (X)

thus :math:`[u]` is a zero of the function :math:`\mathcal{G}`. To solve :eq:`nonlin-linsys`, we use Newton's method. As :math:`J_{\mathcal{G}(X)} = \mathcal{S} - J_{\mathcal{F}(X)}`, the Newton's algorithm writes:


* :math:`X^0` is given,

* knowing :math:`X^n`, we solve : 

.. math::

  J_{\mathcal{G}(X^n)}(X^{n+1} - X^n) = - \mathcal{G}(X^n)

The algorithm is the following:

* we compute the mass matrix associated to the function : :math:`\partial_{u} \mathcal{F}`, \textit{i.e} :

.. math::  

  M_{b,b^{\prime}}^n = \int_{\Omega} \partial_{u} F (\mathbf{x},\sum_{b \in \{1, \cdots, n \}} X_b^n \varphi_b ) \varphi_b \varphi_{b^{\prime}}

* compute the term :math:`\mathcal{F}(X^n)`:

.. math::  
  [\mathcal{F}(X^n)]_{b^{\prime}} = \int_{\Omega} F (\mathbf{x},\sum_{b \in \{1, \cdots, n \}} X_b^n \varphi_b ) \varphi_{b^{\prime}}

* compute :math:`\mathcal{G}(X^n)`:

.. math::
  :label: nonlin-g-fct

  \mathcal{G}(X^n) = \mathcal{S} X^n - \mathcal{F}(X^n)

* compute :math:`J_{\mathcal{G}(X^n)}`:

.. math::  

  J_{\mathcal{G}(X^n)} = \mathcal{S} - J_{\mathcal{F}(X^n)} = \mathcal{S} - M^n

* solve :math:`J_{\mathcal{G}(X^n)}(X^{n+1} - X^n) = - \mathcal{G}(X^n)`, and then find :math:`X^{n+1}`


.. literalinclude:: include/demo/test_nonlin_ex1_newton.py
    :linenos:
    :language: python
    :lines: 42-93

The complete script can be found here

See :download:`the complete script <include/demo/test_nonlin_ex1_newton.py>`.

Which gives the numerical solution

.. image:: include/demo/test_nonlin_ex1_newton.png
   :width: 15cm
   :height: 9cm


.. Local Variables:
.. mode: rst
.. End:
