# -*- coding: UTF-8 -*-
import numpy as np
import pigasus.fem.common_obj as _com
from .solver import *
from scipy.sparse.linalg import splu as slv

__author__ = ['ARA']
__all__ = ['superLU']

class superLU(solver):
    def __init__(self, *args, **kwargs):
        solver.__init__(self, *args, **kwargs)

        # get the matrix in csr format
        self.A = self.matrix.get().tocsc()
        self.A_op = slv ( self.A )

    def solve(self, rhs, **kwargs):
        if _com.isField(rhs):
            lpr_rhs = rhs.get()
        if _com.isNumpyArray(rhs):
            lpr_rhs = rhs

        try:
            F = kwargs["field"]
        except:
            import pigasus.fem.field as fi
            F = fi.field.__new__(fi.field)
            F.id = None

        F.set(self.A_op.solve(lpr_rhs))

        return F

if __name__ == '__main__':
    import caid.cad_geometry  as cg
    from caid.cad_geometry import bilinear
    from pigasus.gallery.EllipticPDE import *
    import pylab                as pl
    import numpy                as np

    # ...
    # SINGLE PRECISION
    atol    = 1e-8
    rtol    = 1e-7
    # DOUBLE PRECISION
    #atol    = 1e-15
    #rtol    = 1e-14

    niter 	= 5000
    # ...

    # ...
    sin = np.sin ; cos = np.cos ; exp = np.exp ; sqrt = np.sqrt ; pi = np.pi
    # ...

    # ...
    testcase = {}
    testcase['list_DirFaces'] = [[1,2,3,4]]

    testcase['stiffness']  = lambda x,y : [1., 0., 0., 1.]

    kx = 2. * pi
    ky = 2. * pi

    # ...
    # exact solution
    # ...
    u = lambda x,y : sin ( kx * x ) * sin ( ky * y )
    testcase['u'] = u
    # ...

    # ...
    # rhs
    # ...
    f = lambda x,y : ( kx**2 + ky**2 ) * sin ( kx * x ) * sin ( ky * y )
    testcase['f'] = f
    # ...

    # ...
    points = np.asarray([[[0.,0.],[0.,1.]],[[1.,0.],[1.,1.]]])
    nrb = bilinear(points)
    geo = cg.cad_geometry(geo=nrb)
    geo.refine(id=0,list_p=[2, 2])
    tx = np.linspace(0.,1.,7)[1:-1]
    ty = np.linspace(0.,1.,7)[1:-1]
    geo.refine(id=0, list_t=[tx, ty])
    # ...

    PDE = EllipticPDE(geometry=geo, testcase=testcase)

    PDE.assembly()

    # ...
    sl = superLU(matrix=PDE.S_V)
    print(sl)

    size = PDE.V.size
    rhs = np.random.random(size)
    F = sl.solve(rhs)
    print(F.get())
    # ...
