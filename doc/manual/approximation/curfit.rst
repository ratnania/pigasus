.. role:: envvar(literal)
.. role:: command(literal)
.. role:: file(literal)
.. _approximation.curfit:

Curve fitting
*************

In the sequel, we show how to construct a smooth B-spline curve given a set of input data. The first thins is to import the corresponding module::

   from pigasus.fit.curfit import curfit

*curfit* is the module you need for this purpose. The main function that we will use is the *construct* function. It can be called as the following::

   geo_f = fit.construct([xk, yk], uk=lists_uk)

As you can see, it takes as inputs the cloud points coordinates :math:`(x_k,y_k)` and the corresponding parametrization :math:`u_k`. These parameters can be computed using::

   from pigasus.fit.curfit import compute_uk
   #method = "uniform"
   #method = "centripetal"
   method = "chord"
   list_Q = zip(list_x, list_y)
   uk = compute_uk(list_Q, method=method)

As a first example, we consider the approximation of an arc of circle. We can enforce an interpolation constraint the extremeties of the curve (in the points **A** and **B**). This is done by calling the function *MakeConstraint* as the following::

   # ... Define the first point A on the face = 0
   A = [1.,0.] ; face = 0
   constraint = MakeConstraint("C0", face, A)
   constraints.append(constraint)

   # ... Define the last point B on the face = 1
   B = [0.,1.] ; face = 1
   constraint = MakeConstraint("C0", face, B)
   constraints.append(constraint)

Where the function *MakeConstraint* is

.. code-block:: python

   def MakeConstraint(cond, face=None, value=None):
       if cond.lower() == "closed":
           constraint = {}
           constraint['patch_id_m'] = 0
           constraint['face_m']     = 0
           constraint['patch_id_s'] = 0
           constraint['face_s']     = 1
       if cond.lower() == "c0":
           constraint = {}
           constraint['patch_id_m'] = 0
           constraint['face_m']     = face
           constraint['type']       = "C0"
           constraint['values']     = [value]
       if cond.lower() == "c1":
           constraint = {}
           constraint['patch_id_m'] = 0
           constraint['face_m']     = face
           constraint['type']       = "C1"
        constraint['values']     = [value]
        return constraint   

These constraints must be given when constructing the *fit* object. The user must also provide the logical domain (*patch*) which is defined as::

   from igakit.cad_geometry import line as patch
   geo = patch(n=[nx], p=[px])

The *fit* object is then constructed as the following::

   from pigasus.fit.curfit import curfit
   fit = curfit(geometry=geo, constraints=constraints, alpha=alpha)

The parameter **alpha** measures the smoothness of the curve.   

.. image:: include/approximation/ex1.png
   :width: 14cm
   :height: 14cm

The complete script can be found here :download:`script <include/approximation/ex1.py>`

In the following example we show how to approximate a regular *closed* curve.

.. image:: include/approximation/ex2.png
   :width: 14cm
   :height: 14cm

The complete script can be found here :download:`script <include/approximation/ex2.py>`

.. todo:: add C1 conditions.

In the following example we show how to approximate a singular *closed* curve (like ones with an X-point).

.. image:: include/approximation/ex3.png
   :width: 14cm
   :height: 14cm

The complete script can be found here :download:`script <include/approximation/ex3.py>`

We Have enforced an interpolation constraint the extremeties of the curve. This is done by calling the function *MakeConstraint* as the following::

   # ... Define the first point A on the face = 0
   A = [0.,0.] ; face = 0
   constraint = MakeConstraint("C0", face, A)
   constraints.append(constraint)

.. note:: As we can expect, the approximation is not good where the number of input data is not suffisant.

Using a closed curve, we can generate the corresponding *2D* description, which leads to the following mesh.

.. image:: include/approximation/ex4.png
   :width: 16cm
   :height: 14cm

The *cad_geometry* module contains a function for this purpose::

  xc = 0.65 ; yc = 0.
  geo = geo_f.polarExtrude(xyzc=[xc,yc])

The complete script can be found here :download:`script <include/approximation/ex4.py>`  

We notice that the mesh is not good where we have less control points for the boundary. This can be corrected by inserting new knots to have a better resolution.

.. Local Variables:
.. mode: rst
.. End:
