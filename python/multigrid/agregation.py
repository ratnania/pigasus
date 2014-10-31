# -*- coding: UTF-8 -*-
__all__ = ['agregation']

import numpy as np
import scipy
from scipy.sparse import csr_matrix, coo_matrix, isspmatrix_csr, isspmatrix_csc
from pyamg.relaxation import gauss_seidel
#from pyamg.util.linalg import residual_norm
# ...
try:
    from petsc4py import PETSc
    importPETSc=True
except ImportError:
    importPETSc=False
# ...
from scipy.sparse.linalg import LinearOperator

norm = np.linalg.norm

# ...
def residual_norm(A, x, b):
    """Compute ||b - A*x||"""

    return norm(b - A.dot(x))
# ...


class level:
    """Stores one level of the multigrid hierarchy

    All level objects will have an 'A' attribute referencing the matrix
    of that level.  All levels, except for the coarsest level, will
    also have 'P' and 'R' attributes referencing the prolongation and
    restriction operators that act between each level and the next
    coarser level.

    Attributes
    ----------
    A : csr_matrix
        Problem matrix for Ax=b
    R : csr_matrix
        Restriction matrix between levels (often R = P.T)
    P : csr_matrix
        Prolongation or Interpolation matrix.

    Notes
    -----
    The functionality of this class is a struct
    """
    def __init__(self, withscipy=False, withPETSc=False):
        from pigasus.fem.solver import solver
#        from pigasus.solver.solver import solver
        if withscipy:
            from pigasus.fem.matrix import matrix
            self.R = None
            self.A = None
            self.P = None
            self._A = matrix()
            self.slv = solver(matrix=self._A)
        else:
            from pigasus.fem.matrix import matrix
            from pigasus.fem.constants import SPM_SOLVER_BASIC_CG, SPM_SOLVER_BASIC_GS
            self.R   = matrix()
            self.A   = matrix()
            self.P   = matrix()
            #
            if withPETSc:
                self.slv = None
                self.smt = None
#                self.smt = solver(matrix=self.A, solver=SPM_SOLVER_BASIC_GS)
            else:
                slvInfo = {}
                slvInfo['solver'] = SPM_SOLVER_BASIC_CG
                self.slv = solver(matrix=self.A, solverInfo=slvInfo)

                smtInfo = {}
                smtInfo['solver'] = SPM_SOLVER_BASIC_GS
                self.smt = solver(matrix=self.A, solverInfo=smtInfo)

        self.withscipy = withscipy
        self.withPETSc = withPETSc

    def set(self, A, R, P):
        """
        sets the scipy matrices A, R, P into pigasus
        """
        if self.withscipy:
            self.R =R
            self.A =A
            self.P =P
            self._A.set(A)
        else:
            self.R.set(R)
            self.A.set(A)
            self.P.set(P)
            self._A = A # TODO: to remove. needed for the moment for the smoother
            if self.withPETSc and importPETSc:
                self.A_petsc = PETSc.Mat().createAIJ(size=A.shape,csr=(A.indptr,A.indices,A.data))
                # ...

                # Initialize ksp solver.
                self.ksp_slv = PETSc.KSP().create()
                self.ksp_slv.setOperators(self.A_petsc)
                self.ksp_slv.setFromOptions()
                self.ksp_slv.setType(PETSc.KSP.Type.CG)

#                # Initialize ksp smoother.
#                self.ksp_smt = PETSc.KSP().create()
#                self.ksp_smt.setOperators(self.A_petsc)
#                self.ksp_smt.setFromOptions()
#                self.ksp_smt.setType(PETSc.KSP.Type.RICHARDSON)

    def solve(self, b, maxiter=6000, tol=1.e-14):
        if self.withPETSc and importPETSc:
            _b = PETSc.Vec().createWithArray(b, comm=PETSc.COMM_SELF)
            _x = PETSc.Vec().createWithArray(np.zeros_like(b), comm=PETSc.COMM_SELF)

            self.ksp_slv.rtol = tol
#            self.ksp_slv.setConvergenceHistory()
            self.ksp_slv.solve(_b, _x)
            return _x.getArray()
        else:
            return self.slv.solve(b, guess=np.zeros_like(b) \
                        , maxiter=maxiter, eps=tol)

    def smoother(self, x, b, nu, tol=1.e-10):
        if self.withscipy:
            gauss_seidel(self.A, x, b, iterations=nu)
            return x
        else:
#            if self.withPETSc and importPETSc:
            if False:
                _b = PETSc.Vec().createWithArray(b, comm=PETSc.COMM_SELF)
                _x = PETSc.Vec().createWithArray(np.zeros_like(b), comm=PETSc.COMM_SELF)

                self.ksp_smt.rtol = tol
                self.ksp_smt.max_it = nu
    #            self.ksp_smt.setConvergenceHistory()
                self.ksp_smt.solve(_b, _x)
                return _x.getArray()
            else:
                gauss_seidel(self._A, x, b, iterations=nu) # TODO: to remove
                return x
#            return self.smt.solve(b, guess=np.zeros_like(b) \
#                        , maxiter=nu, eps=tol)

class agregation(object):

    def __init__(self, list_geometry, gamma, nu1, nu2 \
                 , smoother=None, coarse_solver=None \
                 , withscipy=False, withPETSc=False):
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
        >>> from pigasus.gallery import EllipticPDE
            nrb_H = self.geometry_H[0]
            nrb_h = self.geometry_h[0]


            knots_H = nrb_H.knots[0]
            knots_h = nrb_h.knots[0]

            n = nrb_H.shape[0]
            p = nrb_H.degree[0]

            list_r = [r for r in knots_h if r not in knots_H]

        See Also
        --------
        TODO
        """
        # TODO : for the moment we only treate 1-patch geometries
        self.withscipy      = withscipy
        self.withPETSc      = withPETSc

        self.geometries     = list_geometry
        self.dim            = self.geometries[0].dim
        self.nlevels        = len(self.geometries)-1
        self.coarse_solver  = coarse_solver
        self.smoother       = smoother
        self.gamma          = gamma
        self.nu1            = nu1
        self.nu2            = nu2
        self.nloop          = 1
        self.levels         = []
        for i in range(0, len(self.geometries)):
            self.levels.append(level(withscipy=withscipy))

        self.list_allresiduals  = [] # list of residuals for each step
        self.list_coarseresiduals = []
    #-----------------------------------

    #-----------------------------------
    def initialize(self, A, list_A=None, list_DirFaces=None, GalerkinCoarseGrid=True):
        from scipy.sparse import identity as Id
        from pigasus.fem.common_obj import isMatrix

        if isMatrix(A):
            self.A              = A.get()
        else:
            self.A              = A

        if list_A is not None:
            self.list_A     = list_A
            list_A          = self.list_A[::-1]

        n,m = self.A.shape

        ilvl = 0
        lvl = self.levels[ilvl]
        lvl.set(self.A, Id(n), Id(n))

        geometries  = self.geometries[::-1]
        for (geo_h, geo_H) in zip(geometries[:-1], geometries[1:]):
            ilvl += 1
            # ... interpolator
            P = self.constructInterpolationMatrix(geo_H, geo_h, list_DirFaces=list_DirFaces)
            # ... restriction
            list_H, list_h = self.compute_H(geo_H, geo_h)
            r = list_h[0]/list_H[0]
            R = P.transpose().tocsr()
            R *= r
#            print "Interpolator : ", P.shape
#            print "Restrictor : ", R.shape
            # ... the coarse system
            if GalerkinCoarseGrid:
#                print "Galerkin coarse grid operator has been initialized"
                if self.withscipy:
                    A_h = self.levels[ilvl-1].A
                else:
                    A_h = self.levels[ilvl-1].A.get()
#                print R.shape, A_h.shape, P.shape
                A_H = R * A_h * P
#                print A_h.shape, A_H.shape
#                A_H = A_H.tocsr()
            else:
                if self.withscipy:
                    A_H = list_A[i]
                else:
                    A_H = list_A[i].get()

            lvl = self.levels[ilvl]
            lvl.set(A_H, R, P)

        self.levels = self.levels[::-1]
        A = self.levels[-1].A.get()
        self.dtype = A.dtype
        self.shape = A.shape
    #-----------------------------------

    #-----------------------------------
    def interpolation(self, level, vH):
        P = self.levels[level].P
        vh = P.dot(vH)
        return vh
    #-----------------------------------

    #-----------------------------------
    def restriction(self, level, vh):
        R = self.levels[level].R
        vH = R.dot(vh)
        return vH
    #-----------------------------------

    #-----------------------------------
    def mgcyc(self, k, gamma, ukm, fk, nu1, nu2):
        """
        this routine will retrurn uk_{m+1} using ukm
        """
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
            if self.coarse_solver is None:
                vk1m = lvl.solve(dk1m)
            else:
                vk1m = self.coarse_solver(dk1m)

        if k > 1:
            a = np.zeros_like(dk1m)
            vk1m_ = dk1m
            for i in range(0, gamma):
                dk1m_ = vk1m_
                vk1m_, err_ = self.mgcyc(k-1, gamma, a, dk1m_, nu1, nu2)
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
        # ----------------------------

        err = residual_norm(Lk, ukp1m, fk)

        return ukp1m, err
    #-----------------------------------

    #-----------------------------------
    def solve(self, b, x0=None, tol=1e-7, maxiter=20, cycle='V' \
              , residuals=None, callback=None \
              , accel=None, maxiter_prec=2, tol_prec=1e-7 \
              ):
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

        cycle = str(cycle).upper()

        if accel is not None:

#            # Check for AMLI compatability
#            if (accel != 'fgmres') and (cycle == 'AMLI'):
#                raise ValueError('AMLI cycles require acceleration (accel) to be fgmres, or no acceleration')

            # Acceleration is being used
            if isinstance(accel, str):
                from pyamg import krylov
                from scipy.sparse.linalg import isolve

                if hasattr(krylov, accel):
                    accel = getattr(krylov, accel)
                else:
                    accel = getattr(isolve, accel)

            A = self.levels[-1].A.get()

            M = self.aspreconditioner(cycle=cycle, tol=tol_prec, maxiter=maxiter_prec)

            try:  # try PyAMG style interface which has a residuals parameter
                return accel(A, b, x0=x0, tol=tol, maxiter=maxiter, M=M,
                             callback=callback, residuals=residuals)[0]
            except:  # try the scipy.sparse.linalg.isolve style interface,
                     # which requires a call back function if a residual
                     # history is desired

                cb = callback
                if residuals is not None:
                    residuals[:] = [residual_norm(A, x, b)]
                    def callback(x):
                        if scipy.isscalar(x):
                            residuals.append(x)
                        else:
                            residuals.append(residual_norm(A, x, b))
                        if cb is not None:
                            cb(x)

                return accel(A, b, x0=x0, tol=tol, maxiter=maxiter, M=M,
                             callback=callback)[0]

        else:

            # Scale tol by normb
            normb = norm(b)
            if normb != 0:
                tol = tol * normb
#            print ">>>>>>>>>>>>>> tol ", tol, normb

        if residuals is None:
            residuals = []
        else:
            residuals[:] = []

        residuals.append(residual_norm(self.A, x, b))

        self.first_pass = True

        self.nloop = 0
        while len(residuals) <= maxiter and residuals[-1] > tol:
            x, err = self.mgcyc(self.nlevels, self.gamma, x, b, self.nu1, self.nu2)

            residuals.append(err)

            self.first_pass = False

            if callback is not None:
                callback(x)

            self.nloop += 1

        return x

    def aspreconditioner(self, cycle='V', maxiter=1, tol=1e-12):
        """Create a preconditioner using this multigrid cycle

        Parameters
        ----------
        cycle : {'V','W','F','AMLI'}
            Type of multigrid cycle to perform in each iteration.

        Returns
        -------
        precond : LinearOperator
            Preconditioner suitable for the iterative solvers in defined in
            the scipy.sparse.linalg module (e.g. cg, gmres) and any other
            solver that uses the LinearOperator interface.  Refer to the
            LinearOperator documentation in scipy.sparse.linalg

        See Also
        --------
        multilevel_solver.solve, scipy.sparse.linalg.LinearOperator

        Examples
        --------
        >>> from pyamg.aggregation import smoothed_aggregation_solver
        >>> from pyamg.gallery import poisson
        >>> from scipy.sparse.linalg import cg
        >>> from scipy import rand
        >>> A = poisson((100, 100), format='csr')          # matrix
        >>> b = rand(A.shape[0])                           # random RHS
        >>> ml = smoothed_aggregation_solver(A)            # AMG solver
        >>> M = ml.aspreconditioner(cycle='V')             # preconditioner
        >>> x, info = cg(A, b, tol=1e-8, maxiter=30, M=M)  # solve with CG

        """
        shape = self.shape
        dtype = self.dtype

        def matvec(b):
            return self.solve(b, maxiter=maxiter, cycle=cycle, tol=tol)

        return LinearOperator(shape, matvec, dtype=dtype)

    #-----------------------------------
    def constructInterpolationMatrix(self, geo_H, geo_h, list_DirFaces=None):
        if self.dim ==1:
            from .splineRefMat import constructCurveMatrix as constructMatrix
        if self.dim ==2:
            from .splineRefMat import constructSurfaceMatrix as constructMatrix
        if self.dim ==3:
            print("initInterpolation: Not yet implemented for 3D")

        patch_id = 0
        nrb_H = geo_H[patch_id]
        nrb_h = geo_h[patch_id]

        if list_DirFaces is None:
            DirFaces = []
        else:
            DirFaces = list_DirFaces[patch_id]

        if self.dim ==1:
            knots_H = nrb_H.knots[0]
            knots_h = nrb_h.knots[0]

            n = nrb_H.shape[0]
            p = nrb_H.degree[0]

            list_r = [r for r in knots_h if r not in knots_H]

            M = constructMatrix(list_r, p, n, knots_H \
                                , DirFaces=DirFaces)

            return M

        if self.dim ==2:
            u_H1,u_H2   = nrb_H.knots
            n_H1,n_H2   = nrb_H.shape
            p_H1,p_H2   = nrb_H.degree

            u_h1,u_h2   = nrb_h.knots

            list_r1 = [r for r in u_h1 if r not in u_H1]
            list_r2 = [r for r in u_h2 if r not in u_H2]

            M, [n,m] = constructMatrix( list_r1, list_r2 \
                                , p_H1, p_H2 \
                                , n_H1, n_H2 \
                                , u_H1, u_H2 \
                                , DirFaces=DirFaces)
            return M

    def compute_H(self, geo_H, geo_h):
        dim = geo_H.dim

        list_H = []
        list_h = []

#        for (nrb_H, nrb_h) in zip(geo_H, geo_h):
        nrb_H = geo_H[0]
        nrb_h = geo_h[0]
        H = 1.
        h = 1.
        for d in range(0,dim):
            p_H = nrb_H.degree[d]
            p_h = nrb_h.degree[d]
            H *= nrb_H.knots[d][p_H+1]-nrb_H.knots[d][0]
            h *= nrb_h.knots[d][p_h+1]-nrb_h.knots[d][0]

        list_H.append(H)
        list_h.append(h)

        return list_H, list_h

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
