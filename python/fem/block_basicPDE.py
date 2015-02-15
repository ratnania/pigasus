# -*- coding: UTF-8 -*-
#! /usr/bin/python

# ...
import numpy                as np
from pigasus.fem.basicPDE import *
from pigasus.utils.blockdata import *
#from pigasusClass import pigasus
from scipy.sparse.linalg import cg

# ...
class block_basicPDE():

    def __init__(self, size, dict_testcases=None, dict_PDE=None, geometry=None, same_space=True):
        """
        a block of basicPDEs. size must a list of two integers, describing the
        number of rows and columns of the block PDEs
        """
#        pigasus.__init__(self)

        # TODO PDEs with different spaces
        if not same_space:
            print("block_basicPDE: NOT YET IMPLEMENTED when same_space=False")
        self._size          = size
        self._dict_PDE      = {}
        self._list_PDE      = []
        self._system        = None
        self._dict_system   = {}
        self._list_rhs      = []
        self._list_unknown  = []
        self._list_norm     = []
        self._dict_testcases= dict_testcases
        self._geometry      = geometry

        if dict_PDE is not None:
            self._dict_PDE  = dict_PDE

        if (dict_testcases is None) and (dict_PDE is None):
            raise ("You should specify either dict_testcases or dict_PDE")

        if (dict_testcases is not None) and (dict_PDE is not None):
            raise ("You should specify either dict_testcases or dict_PDE")

        if dict_testcases is not None:
            # ... first, we create PDEs correponsing to the testcases dictionary
            #     the result is also a dictionary of PDEs
            iteration = 0
            PDE = None
            V = None
            for keys, tc in dict_testcases.items():
                i = keys[0] ; j = keys[1]
                if (iteration == 0) or not same_space:
                    PDE = basicPDE(geometry=geometry, testcase=tc)
                    V   = PDE.V
                else:
                    PDE = basicPDE(geometry=geometry, testcase=tc, V=V)
                self._dict_PDE[i,j] = PDE
                iteration += 1
            # ...

        # ... second, we create a list of PDEs, and initialize all PDE to None
        for i in range(0,self.size[0]):
            line = []
            for j in range(0,self.size[1]):
                line.append(None)
            self._list_PDE.append(line)
        # ...

        # ... then, initialize the PDEs (double) list with the PDEs dictionary
        for keys, PDE in self._dict_PDE.items():
            i = keys[0] ; j = keys[1]
            self._list_PDE[i][j] = PDE
        # ...

        # ... create the fem reference
        for j in range(0,self.size[1]):
            self._fem = None
            # ... find a non-None PDE in the same block column as M
            for i in range(0,self.size[0]):
                PDE = self._list_PDE[i][j]
                if PDE is not None:
                    self._fem = PDE.fem
                    break
            if self._fem is not None:
                break
        # ...

        # ... create the list of rhs and unknowns
        for j in range(0,self.size[1]):
            rhs = None
            # ... find a non-None PDE in the same block column as M
            for i in range(0,self.size[0]):
                PDE = self._list_PDE[i][j]
                if PDE is not None:
                    U   = PDE.unknown
                    rhs = PDE.rhs
                    break

            self._list_rhs.append(rhs)
        # ...

        # ... create the list of unknowns and their norms
        for i in range(0,self.size[0]):
            U   = None
            N_U = None
            # ... find a non-None PDE in the same block line as M
            for j in range(0,self.size[j]):
                PDE = self._list_PDE[i][j]
                if PDE is not None:
                    U   = PDE.unknown
                    N_U = PDE.N_U
                    break

            self._list_unknown.append(U)
            self._list_norm.append(N_U)
        # ...


    @property
    def fem(self):
        return self._fem

    @property
    def geometry(self):
        return self._geometry

    @property
    def size(self):
        return self._size

    @property
    def system(self):
        return self._system

    @property
    def rhs(self):
        return self._list_rhs

    @property
    def unknowns(self):
        return self._list_unknown

    @property
    def PDE(self, i, j):
        return self._list_PDE[i][j]

    def norms(self, list_exact=None):
        if list_exact is not None:
            for (N_U, exact) in zip(self._list_norm, list_exact):
                N_U.set_func(exact)

        self.fem.assembly(norms=self._list_norm)
        return [N_U.get() for N_U in self._list_norm]

    def assembly(self):
        # TODO take into account the case where a PDE is not specified and we
        # have to put a zero matrix
        self._dict_system = {}
        for keys, PDE in self._dict_PDE.items():
#            print ">>> Assembly PDE ", keys
            PDE.assembly()
            i = keys[0] ; j = keys[1]
            self._dict_system[i,j] = PDE.system

        matrices = []
        for i in range(0,self.size[0]):
            line = []
            for j in range(0,self.size[1]):
                line.append(None)
            matrices.append(line)

        for keys, system in self._dict_system.items():
            i = keys[0] ; j = keys[1]
            matrices[i][j] = system.get()

        list_rhs = [rhs.get() for rhs in self._list_rhs]

        # ... create block matrix
#        print ">>> create block matrix"
        self._system = BlockMatrix(matrices)
#        print "<<< done."
#        print ">>> assembly block matrix"
        self._system.assembly()
#        print "<<< done."

        # ... create block vector
#        print ">>> create block vector"
        self._rhs= BlockVector(list_rhs)
#        print "<<< done."
#        print ">>> assembly block vector"
        self._rhs.assembly()
#        print "<<< done."

    def free(self):
        for keys, PDE in self._dict_PDE.items():
            PDE.free()



    #-----------------------------------
    def norm(self, exact=None):
        """
        Computing the L2 norm
        """
        if exact is not None:
            self.N_U.set_func(exact)

        self.fem.assembly(norms=[self.N_U])

        return self.N_U.get()
    #-----------------------------------

    #-----------------------------------
    def solve(self, rhs=None):
        """
        solves the system assemblied after calling the assembly function.
        The solution is therefor stored in self.U_V (Homogeneous Dirichlet bc)
        or self.U_W (Dirichlet bc). The user can get it using the function get

        rhs:
           Can be either a list of field objects, a list of numpy arrays or a
           numpy array.
           If not given, self.rhs
           will be used as rhs.
        """
#        from pigasus.utils.utils import process_diagnostics
#        print ">> solve"
#        process_diagnostics(30)
        # ...
        lpr_rhs = None
        if rhs is None:
            lpr_rhs = self.rhs.get()
        else:
            # ... TODO concatenation using is not optimal => use numpy
            from common_obj import isNumpyArray, isField
            if isNumpyArray(rhs[0]):
                L = []
                for F in rhs:
                    L += list(F)
                lpr_rhs = np.asarray(L)
            if isField(rhs[0]):
                L = []
                for F in rhs:
                    L += list(F.get())
                lpr_rhs = np.asarray(L)
            if isNumpyArray(rhs):
                lpr_rhs = rhs
        # ...

        tol = 1.e-7
        maxiter = 1000
        A = self.system.get()
#        print A.shape, lpr_rhs.shape
        X = cg(A, lpr_rhs, tol=tol, maxiter=maxiter)[0]

#        print "###"
#
#        print "*****"
#        print _b
#        print "*****"
#        print X
#        print "*****"

        ind_b = 0
        for U in self.unknowns:
            U.set(X[ind_b:ind_b+U.size])
            ind_b += U.size
        # ...


    #-----------------------------------


if __name__ == "__main__":
    from caid.cad_geometry import square as domain
    from numpy import pi, sin
    nx = 15 ; ny = 15
    px = 2 ; py = 2
    geo = domain(n=[nx,ny],p=[px,px])


    #-----------------------------------
    def testcase():
        # implicit part
        tc = {}

        kx = 2. * pi
        ky = 2. * pi

        # ... exact solution
        u = lambda x,y : [sin ( kx * x ) * sin ( ky * y )]
        # ... rhs
        f = lambda x,y : [( kx**2 + ky**2 ) * sin ( kx * x ) * sin ( ky * y )]

        A = lambda x,y : [  1., 0. \
                          , 0., 1. ]

        tc['u']  = u
        tc['f']  = f
        tc['A']  = A

        tc['AllDirichlet'] = True

        return tc
    #-----------------------------------

    size = [2,2]

    dict_testcases = {}
    dict_testcases[0,0] = testcase()
#    dict_testcases[0,1] = testcase()
#    dict_testcases[1,0] = testcase()
    dict_testcases[1,1] = testcase()

    PDEs = block_basicPDE(size, dict_testcases=dict_testcases, geometry=geo)
    PDEs.assembly()

#    rhs = PDEs.rhs.get()
#    system = PDEs.system
#    matrix = system.get()
#
#    Y = matrix.dot(rhs)
#    print Y.shape

    rhs = PDEs.rhs
    PDEs.solve(rhs)

    # ... example of multiplication with a list of fields
    unknowns = PDEs.unknowns
    system = PDEs.system
    X = system.dot(unknowns)
    print([np.linalg.norm(x.get()-r.get()) for (x,r) in zip(X,rhs)])
    # ...

    # ... example of multiplication with a list of numpy arrays
    unknowns = [U.get() for U in PDEs.unknowns]
    system = PDEs.system
    X = system.dot(unknowns)
    print([np.linalg.norm(x-r.get()) for (x,r) in zip(X,rhs)])
    # ...


#    for U in PDEs.unknowns:
#        print U.get()

    print(PDEs.norms())


    PDEs.free()
