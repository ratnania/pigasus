.. role:: envvar(literal)
.. role:: command(literal)
.. role:: file(literal)
.. _gallery.mhdreconnection:


Non-linear Magnetic Reconnection
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We are interested in the following non-linear 2D reconnection model:

.. math::

  [\phi, \psi] &= \hat{\eta} \nabla^2 \psi - \hat{E}
  \\
  [\phi, \nabla^2 \phi] &= \frac{1}{\beta} [\psi, \nabla^2 \psi] + \hat{\nu} \nabla^4 \phi

Let us introduce the variables

.. math::

  j &= \nabla^{2} \psi
  \\
  w &= \nabla^{2} \phi

The problem now writes

.. math::

  [\phi, \psi] &= \hat{\eta} \nabla^2 \psi - \hat{E}
  \\
  [\phi, w] &= \frac{1}{\beta} [\psi, j] + \hat{\nu} \nabla^2 w  
  \\
  j &= \nabla^2 \psi
  \\
  w &= \nabla^{2} \phi

Or in an equivalent form

.. math::

  \hat{\eta} \nabla^2 \psi + [\psi, \phi] &= \hat{E}
  \\
  \hat{\nu} \nabla^2 w + [w, \phi] &= \frac{1}{\beta} [j, \psi] 
  \\
  j &= \nabla^2 \psi
  \\
  w &= \nabla^{2} \phi


Picard Algorithm
________________

A simple Picard algorithm can be written as follows

.. math::

  \hat{\eta} \nabla^2 \psi^{n+1} &= [\phi^n, \psi^n] + \hat{E}
  \\
  j^{n+1} &= \nabla^2 \psi^{n+1}
  \\
  \hat{\nu} \nabla^2 w^{n+1} &= [\phi^n, w^n] + \frac{1}{\beta} [j^{n+1}, \psi^{n+1}] 
  \\
  w^{n+1} &= \nabla^{2} \phi^{n+1}


.. Local Variables:
.. mode: rst
.. End:
