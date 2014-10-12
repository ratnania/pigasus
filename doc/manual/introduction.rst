.. role:: envvar(literal)
.. role:: command(literal)
.. role:: file(literal)
.. _introduction:


Getting Started
***************

In this section we give some implemented PDE models and classes that can be used as a starting point to write new scripts. All these models are available via the *pigasus.gallery* module. For the moment, only the **1D** and **2D** cases have been validated. The **3D** case is planned for the pigasus 0.9.9 release.

===========================    =====   =====   ====  
Models                          1D      2D      3D
===========================    =====   =====   ====
basicPDE                        Yes     Yes     No 
parabolicPDE                    Yes     Yes     No 
systemPDE                       No      No      No 
poisson                         Yes     Yes     No 
bilaplacian [#f1]_              Yes     Yes     No 
nonlinear-poisson               Yes     Yes     No 
MongeAmpere                     -       Yes     No 
GradShafranov                   -       Yes     -
AnisotropicDiffusion            -       Yes     No 
AnisotropicDiffusion-MMPDE      -       No      No 
Alignment                       -       Yes     -
ReducedMHD                      -       No      No 
Reconnection                    -       No      -
Geometric-MultiGrid [#f2]_      Yes     Yes     No
===========================    =====   =====   ====

All these models are based on the *basicPDE* class which implements a Partial Differential Operator. Understanding this class may
help you to easily write new scripts.

.. rubric:: Footnotes

.. [#f1] Only Homogeneous Dirichlet Boundary Condition
.. [#f2] Only single patch 

.. Local Variables:
.. mode: rst
.. End:
