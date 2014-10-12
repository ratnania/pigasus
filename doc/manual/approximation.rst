.. role:: envvar(literal)
.. role:: command(literal)
.. role:: file(literal)
.. _approximation:

Approximation
=============

**Pigasus** can be used for least square approximation. In many applications, one needs to construct smooth surfaces (or curves) from a given set of input data. The standard tools like *fitpack* only allow single patch reconstruction. Using **Pigasus**, we will show how to reconstruct a smooth surface using a domain decomposition, depending on a provided topology. 


.. include:: approximation/curfit.rst

.. include:: approximation/surfit.rst

.. Local Variables:
.. mode: rst
.. End:
