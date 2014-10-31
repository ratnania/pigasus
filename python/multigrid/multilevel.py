# -*- coding: UTF-8 -*-
"""Generic GMG solver"""

__docformat__ = "restructuredtext en"

from warnings import warn

import scipy
import numpy


__all__ = ['multilevel_solver']

class level:
    """Stores one level of the multigrid hierarchy

    All level objects will have an 'A' attribute referencing the matrix
    of that level.  All levels, except for the coarsest level, will
    also have 'P' and 'R' attributes referencing the prolongation and
    restriction operators that act between each level and the next
    coarser level.

    Attributes
    ----------
    A :
        Problem matrix for Ax=b
    R : reduction
        Restriction matrix between levels (often R = P.T)
    P : interpolator
        Prolongation or Interpolation matrix.

    Notes
    -----
    The functionality of this class is a struct
    """
    def __init__(self, withPETSc=False):
        from pigasus.fem.matrix import matrix

        self.withPETSc = withPETSc

        self.R   = None
        self.P   = None
        self.A   = matrix()
        if self.withPETSc: # will be treated after
            self.slv = None
            self.smt = None
        else:
            from pigasus.solver.solver import solver
            from pigasus.fem.constants import SPM_SOLVER_BASIC_CG, SPM_SOLVER_BASIC_GS
            self.slv = solver(matrix=self.A, solver=SPM_SOLVER_BASIC_CG)
            self.smt = solver(matrix=self.A, solver=SPM_SOLVER_BASIC_GS)

    def set_P(self, P):
        """
        """
        self.P = P

    def set_R(self, R):
        """
        """
        self.R = R

    def set_A(self, A):
        """
        """
        self.A.set(A)

    def construct(self):
        """
        construct the current level, operators and the coarse matrix
        """
        self.P.construct()
        self.R.construct()

    def solve(self, b, maxiter=6000, tol=1.e-10):
        if self.withPETSc:
            _b = PETSc.Vec().createWithArray(b, comm=PETSc.COMM_SELF)
            _x = PETSc.Vec().createWithArray(np.zeros_like(b), comm=PETSc.COMM_SELF)

            self.ksp_slv.rtol = tol
#            self.ksp_slv.setConvergenceHistory()
            self.ksp_slv.solve(_b, _x)
            return _x.getArray()
        else:
            return self.slv.solve(b, guess=np.zeros_like(b) \
                        , maxiter=maxiter, eps=tol)

    def smoother(self, x, b, iterations=100, tol=1.e-10):
        if self.withPETSc:
            _b = PETSc.Vec().createWithArray(b, comm=PETSc.COMM_SELF)
            _x = PETSc.Vec().createWithArray(np.zeros_like(b), comm=PETSc.COMM_SELF)

            self.ksp_smt.rtol = tol
            self.ksp_smt.max_it = nu
#            self.ksp_smt.setConvergenceHistory()
            self.ksp_smt.solve(_b, _x)
            return _x.getArray()
        else:
            return self.smt.solve(b, guess=np.zeros_like(b) \
                        , maxiter=iterations, eps=tol)

class multilevel_solver:
    """Stores multigrid hierarchy and implements the multigrid cycle

    The class constructs the cycling process and points to the methods for
    coarse grid solves. A call to multilevel_solver.solve() is a typical access point.
    The class also defines methods for constructing operator, cycle, and grid complexities.

    Attributes
    ----------
    levels : level array
        Array of level objects that contain A, R, and P.
    coarse_solver : string
        String passed to coarse_grid_solver indicating the solve type

    Methods
    -------
    cycle_complexity()
        A measure of the cost of a single multigrid cycle.
    grid_complexity()
        A measure of the rate of coarsening.
    operator_complexity()
        A measure of the size of the multigrid hierarchy.
    solve()
        Iteratively solves a linear system for the right hand side.
    """

    def __init__(self, list_geometry, gamma, nu1, nu2, withPETSc=False):
        """Creates a geometric multigrid for the matrix list_A[-1]

        Parameters
        ----------
        list_A : is a list of csr_matrix or pigasus-matrix. list_A[0] is on the finest grid
        list_geometry : list of geometries [-1] -> the finest geometry and [0] -> the coarse
        nlevels : the number of subdomains levels

        Returns
        -------
        mg : the geometric multigrid

        Examples
        --------
        >>> from scipy.sparse import csr_matrix
        >>> from pigasus.gallery import poisson

        See Also
        --------
        TODO
        """
        # TODO : for the moment we only treate 1-patch geometries

        self.withPETSc      = withPETSc
        self.geometries     = list_geometry
        self.dim            = self.geometries[0].dim
        self.nlevels        = len(self.geometries)-1
        self.gamma          = gamma
        self.nu1            = nu1
        self.nu2            = nu2
        self.nloop          = 1
        self.levels         = []
        for i in range(0, len(self.geometries)):
            self.levels.append(level())

        self.list_allresiduals  = [] # list of residuals for each step
        self.list_coarseresiduals = []
    #-----------------------------------

    #-----------------------------------
    def initialize(self, list_A):
        from scipy.sparse import identity as Id
        from pigasus.multigrid.operators import interpolation, reduction

        self.A              = list_A[-1]
        self.list_A         = list_A

        n,m = self.A.shape

        ilvl = 0
        lvl = self.levels[ilvl]
        lvl.set(self.A, Id(n), Id(n))

        geometries  = self.geometries[::-1]
        list_A      = self.list_A[::-1]
        for (geo_h, geo_H) in zip(geometries[:-1], geometries[1:]):
            ilvl += 1
            lvl = self.levels[ilvl]

            # ... interpolator
            P = interpolation(geo_H, geo_h)
            # ... restriction
            R = reduction(geo_H, geo_h)
            # ... the coarse system
            try:
                A_H = list_A[i].get()
            except:
                print("Galerkin coarse grid operator has been initialized")
                A_H = coarse_matrix(geo_H, geo_h, DirFaces=self.DirFaces)
                A_h = self.levels[ilvl-1].A.get()
                A_H.construct(A_h)
#                print A_h.shape, A_H.shape
#                A_H = A_H.tocsr()

            lvl.set_P(P)
            lvl.set_R(R)
            lvl.set_A(A)

        self.levels = self.levels[::-1]
    #-----------------------------------

    #-----------------------------------
    def interpolation(self, level, vH):
        P = self.levels[level].P
        vh = P.apply(vH)
        return vh
    #-----------------------------------

    #-----------------------------------
    def restriction(self, level, vh):
        R = self.levels[level].R
        vH = R.apply(vh)
        return vH
    #-----------------------------------

    #-----------------------------------
    def mgcyc(self, k, gamma, ukm, fk, nu1, nu2 \
              , smoother=None, coarse_solver=None):
        """
        this routine will retrurn uk_{m+1} using ukm
        """
        if smoother is None:
            smoother = self.smoother
        if coarse_solver is None:
            coarse_solver = self.coarse_solver

        nlevels = self.nlevels + 1

        lvl = self.levels[::-1][nlevels-k]
        lvl1 = self.levels[::-1][nlevels-k-1]

        Rk  = lvl.R
        Pk  = lvl.P
        Lk  = lvl1.A
        Lk1 = lvl.A

        # ----------------------------
        # presmoothing
        # ----------------------------
        ukm_s = lvl1.smoother(ukm, fk, nu1)
#        ukm_s = smoother(nu1, ukm, Lk, fk)
        # ----------------------------

        # ----------------------------
        # coarse grid correction
        # ----------------------------
        # Compute the defect
        dkm = fk - Lk.dot(ukm_s)
        # Restrict the defect
        dk1m = Rk.dot(dkm)
        # Compute an approximate solution vk1m of the defect equation on Omega_{k-1}
        # if k = 1, use a direct or fast iterative solver, by calling
        if k == 1:
            # TODO : enlever le guess
            guess = np.zeros_like(dk1m)
            vk1m = lvl.solve(dk1m)
#            vk1m = coarse_solver(Lk1, guess, dk1m)
        if k > 1:
            a = np.zeros_like(dk1m)
            vk1m_ = dk1m
            for i in range(0, gamma):
                dk1m_ = vk1m_
                vk1m_, err_ = self.mgcyc(k-1, gamma, a, dk1m_, nu1, nu2 \
                                  , smoother=smoother \
                                   , coarse_solver=coarse_solver)
            vk1m = vk1m_

        # Interpolate the correction
#        print "vk1m : ", vk1m.__class__.__name__, vk1m.shape
#        print "Pk : ", Pk.__class__.__name__, Pk.shape
        vkm = Pk.dot(vk1m)
        # Compute the corrected approximation
        ukm += vkm
        # ----------------------------

        # ----------------------------
        # postsmoothing
        # ----------------------------
        ukp1m = lvl1.smoother(ukm, fk, nu2)
#        ukp1m = smoother(nu2, ukm, Lk, fk)
        # ----------------------------

        err = residual_norm(Lk, ukp1m, fk)

        return ukp1m, err

    #-----------------------------------

    def solve(self, b, x0=None, tol=1e-5, maxiter=100, cycle='V', residuals=None):
        """Main solution call to execute multigrid cycling.

        Parameters
        ----------
        b : array
            Right hand side.
        x0 : array
            Initial guess.
        tol : float
            Stopping criteria: relative residual r[k]/r[0] tolerance.
        maxiter : int
            Stopping criteria: maximum number of allowable iterations.
        cycle : {'V','W','F'}
            Type of multigrid cycle to perform in each iteration.
        residuals : list
            List to contain residual norms at each iteration.

        Returns
        -------
        x : array
            Approximate solution to Ax=b

        See Also
        --------
        aspreconditioner

        Examples
        --------
        """

        if x0 is None:
            x = np.zeros_like(b)
        else:
            x = np.array(x0)  # copy

        # Scale tol by normb
#        normb = norm(b)
#        if normb != 0:
#            tol = tol * normb

        residuals.append(residual_norm(self.A, x, b))

        self.first_pass = True

        self.nloop = 0
        while len(residuals) <= maxiter and residuals[-1] > tol:
            x, err = self.mgcyc(self.nlevels, self.gamma, x, b, self.nu1, self.nu2)

            residuals.append(err)

            self.first_pass = False

            self.nloop += 1

        return x

    #-----------------------------------


    def __repr__(self):
        """Prints basic statistics about the multigrid hierarchy.
        """
        from pyamg.util.linalg import condest

        levels = self.levels[::-1]

        output = 'multilevel_solver\n'
        output += 'Conditioning Number of the matrix:     %d\n' % condest(self.A)
        output += 'Number of Levels:     %d\n' % len(levels)
        output += 'Operator Complexity: %6.3f\n' % self.operator_complexity()
        output += 'Grid Complexity:     %6.3f\n' % self.grid_complexity()
#        output += 'Coarse Solver:        %s\n' % self.coarse_solver.name()

        total_nnz = sum([level.A.nnz for level in levels])

        output += '  level   unknowns     nonzeros\n'
        for n, level in enumerate(levels):
            A = level.A
            output += '   %2d   %10d   %10d [%5.2f%%]\n' %\
                    (n, A.shape[1], A.nnz,\
                    (100 * float(A.nnz) / float(total_nnz)))

        return output

    def operator_complexity(self):
        """Operator complexity of this multigrid hierarchy

        Defined as:
            Number of nonzeros in the matrix on all levels /
            Number of nonzeros in the matrix on the finest level

        """
        levels = self.levels[::-1]

        return sum([level.A.nnz for level in levels]) /\
                float(levels[0].A.nnz)

    def grid_complexity(self):
        """Grid complexity of this multigrid hierarchy

        Defined as:
            Number of unknowns on all levels /
            Number of unknowns on the finest level

        """
        levels = self.levels[::-1]

        return sum([level.A.shape[0] for level in levels]) /\
                float(levels[0].A.shape[0])

