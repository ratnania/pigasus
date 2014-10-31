# -*- coding: UTF-8 -*-
import numpy as np
from .pigasusObject import *

__author__ = ['ARA']
__all__ = ['solver']


class solver(pigasusObject):
    def __init__(self, matrix=None, solverInfo=None):
        pigasusObject.__init__(self)
        self.__notInfo__ = ["guess","rhs","residual","preconditioner","matrix"]
        self._initialize(matrix=matrix, solverInfo=solverInfo)

        # this must be the last thing to do
        self.id = self.com.nsolvers
        self.com.nsolvers += 1
        self.com.solvers.append(self)

    def _initialize(self, matrix=None, solverInfo=None):
        self.size = 0
        # ...
        self.matrix = matrix
        # ...

        # ...
        try:
            self.solver = solverInfo["solver"]
        except:
            from .constants import SPM_SOLVER_BASIC_CG
            self.solver = SPM_SOLVER_BASIC_CG
        # ...
        # ...
        try:
            self.rhs = solverInfo["rhs"]
        except:
            self.rhs = None
        # ...
        # ...
        try:
            self.residual = solverInfo["residual"]
        except:
            self.residual = None
        # ...
        # ...
        try:
            self.guess = solverInfo["guess"]
        except:
            self.guess = None
        # ...
        # ...
        try:
            self.precondtioner = solverInfo["preconditioner"]
        except:
            self.preconditioner = None
        # ...
        # ...
        try:
            self.maxiter = solverInfo["maxiter"]
        except:
            self.maxiter = 5000
        # ...
        # ...
        try:
            self.eps = solverInfo["eps"]
        except:
            self.eps = 1.e-7
        # ...
        # ...
        try:
            self.constraint = solverInfo["constraint"]
        except:
            self.constraint = None
        # ...

        # ...
        self.residualErr = []
        # ...

        # ...
        self.niter = 0
        self.resnorm = 0.
        self.Initialized = False
        self.Finalized = False
        # ...

    def help(self):
        print("atol : absolute tolerance")
        print("rtol : relative tolerance")
        print("maxiter : maximum iterations for iterative solver")
        print("solver : cg, gmres, pastix, ... ")
        print("guess : the initial guess for iterative sovler")
        print("preconditioner : must be callable and returns a vector for a given vector")
        print("residual : must be callable and returns a vector for a given vector. This can be the product Mv, if we want to solve MX = Y")
        print("rhs : the right hand side")

    def initialize(self):
        if not self.Initialized:
            self.com.pyfem.solver_initialize(self.id)
            self.Initialized = True
#            print 'solver ', self.id, ' has been initialized'

    def finalize(self):
        """
        must be called at the end
        """
        if not self.Finalized:
            self.com.pyfem.solver_finalize(self.id)
            self.Finalized = True
#            print 'solver ', self.id, ' has been finalized'

#    def __del__(self):
#        self.finalize()

    def solve( self, rhs, field=None \
              , guess=None, maxiter=None, eps=None \
             , c_constraint=None):
        if maxiter is not None:
            self.maxiter = maxiter
        if eps is not None:
            self.eps = eps
        # TODO set guess

        self.com.pyfem.set_solver_maxiter(self.id, self.maxiter)
        self.com.pyfem.set_solver_eps(self.id, self.eps)

        if c_constraint is None:
            try:
                c_constraint = self.solverInfo['constraint']
                withConstraint = True
            except:
                withConstraint = False
        else:
            withConstraint = True

        if withConstraint :
            return self.solve_constraint( rhs, c_constraint, field)
        else:
            return self.solve_basic( rhs, field)

    def solve_basic( self, rhs, field):
        self.set_rhs(rhs)

        # ... initialization
        self.initialize()

        # ... Solving
        self.com.pyfem.solver_run(self.id)

        self.niter = self.com.pyfem.get_solver_niter(self.id)
        self.residualErr = self.com.pyfem.get_solver_err(self.id, self.maxiter)
        self.resnorm = self.com.pyfem.get_solver_resnorm(self.id)

        return self.get_solution(field)

    def solve_constraint( self, rhs, c, field):
#        print ">> solve_constraint"
        from scipy.sparse.linalg import cg, LinearOperator
        A = self.matrix.get()
        b = rhs
        n = A.shape[0]
        def mv(v):
            r = A*v + np.dot(c, v) / len(v)
            return r

        lin = LinearOperator((n,n), matvec=mv, dtype=np.float)
        x = cg(lin, b, tol=self.eps)[0]
#        print 'done'
        return x

    def set_rhs(self, rhs):
        from .common_obj import isField, isNumpyArray
        if isField(rhs):
            self.com.pyfem.solver_setfieldrhs(self.id, rhs.id )
        if isNumpyArray(rhs):
            self.com.pyfem.solver_setvaluesrhs(self.id, rhs )

    def get_solution(self, dest):
        from .common_obj import isField, isNumpyArray
        if isField(dest):
            self.com.pyfem.solver_getfieldsolution ( self.id, dest.id )
            return None
        else:
            if self.matrix is not None:
                self.size = self.matrix.shape[0]
            sol = self.com.pyfem.solver_getsolution ( self.id , self.size )
            return sol


if __name__ == '__main__':
    x = 0.

    sl = solver(guess=0.1, maxiter=100)
#    print(sl)
