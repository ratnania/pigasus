.. role:: envvar(literal)
.. role:: command(literal)
.. role:: file(literal)
.. _figa:


Fast Solver for Tokamak geometries
==================================

A fast solver has been implemented in **pigasus**. It is based on the Kronecker product and *FFT*. Our solver treats matrices of the form

.. math::

   \Sigma = \sum_{k=1}^r \mathcal{A}_k \otimes B_k

where the matrices :math:`\left( \mathcal{A}_k \right)_{1 \leq k \leq r}` are circulants.


The latter two matrices :math:`K_\theta` and :math:`M_\theta` are circulant, which means that they 
can be both diagonalized in the same orthonormal basis corresponding to the Fourier modes.
This can be expressed by
:math:`M_\theta= P\Lambda_M P^*, ~~~~~ K_\theta= P\Lambda_K P^*,`
where :math:`\Lambda_M` and :math:`\Lambda_K` are the diagonal matrices of the eigenvalues and a multiplication by :math:`P` corresponds to the normalized Fast Fourier Transform and a multiplication by :math:`P^*` to its inverse.

This can be exploited for the fast solution of :ref:`kronsyst` bringing the solution of a linear system arising from a 2D problem to sets of smaller systems corresponding to 1D problems.
The procedure can be performed with the following algorithm:

*  Multiply system :ref:`kronsyst` on the right by :math:`P` (amounts to a 1D FFT on each line of :math:`F`). Then we get

.. math::
   :label: firststep

   K_{ar}UP\Lambda_M + M_{ar}UP\Lambda_K + M_{cr}UP\Lambda_M = FP.

* Note that a multiplication on the right by the diagonal matrix of eigenvalues corresponds to multiplying each column of the matrix by the corresponding eigenvalue, which implies that :ref:`firststep` corresponds to uncoupled problems on each of the columns of :math:`\hat{U}=UP`.

So denoting  by :math:`\hat{U}_1,\dots, \hat{U}_n` the columns of the matrix :math:`\hat{U}` and by 
:math:`\hat{F}_1,\dots, \hat{F}_n` the columns of the matrix :math:`\hat{F}=FP`, :ref:`firststep`
becomes for each column :math:`j`,

.. math::

   (\lambda_{M_j}(K_{ar}+ M_{cr}) + \lambda_{K_j}M_{ar})\hat{U}_j=\hat{F}_j 

This is a set of banded systems of size the number of points in the :math:`r` direction, that can be solved very efficiently using the LAPACK routines DPBTRF for the Cholesky factorization that is only called once at the beginning and then DPBTRS for the solution at each time step.

*  Compute :math:`U=\hat{U} P^*` by inverse FFT of the lines of :math:`\hat{U}`.

Let us now compute the cost at each time step for a :math:`N_r\times N_\theta` mesh, disregarding the cost of the Cholesky factorization for the systems in :math:`r` which needs only to be performed once for a many time steps computation.
The algorithm consists of three steps: 

1. :math:`N_r` FFTs which need :math:`O(N_\theta\log_2N_\theta)` operations,  

2. :math:`N_\theta` up and down sweeps of a Cholesky decomposed banded system which cost :math:`O(N_r)` each, 
   
3. :math:`N_r` inverse FFTs which cost :math:`O(N_\theta\log_2N_\theta)` operations. 
   
So all together the cost is :math:`O(N_r N_\theta\log_2N_\theta)` operations, which is almost optimal. This algorithm uses the structure of the system in an optimal manner and only works on dense matrices. A generic sparse systems solvers could not do this.


.. Local Variables:
.. mode: rst
.. End:
