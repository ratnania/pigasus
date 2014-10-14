# -*- coding: UTF-8 -*-
import numpy as np
from caid.cad_geometry import square as domain
from caid.core.bspline import bsp

# ------------------------------
def findSpan(list_x, list_t, list_p, list_n):
    x = list_x[0]
    y = list_x[1]

    tx = list_t[0]
    ty = list_t[1]

    px = list_p[0]
    py = list_p[1]

    nx = list_n[0]
    ny = list_n[1]

    leftmkx = bsp.FindSpan(px,tx,x) - px
    leftmky = bsp.FindSpan(py,ty,y) - py

    return leftmkx,leftmky
# ------------------------------

# ------------------------------
def evalSplines(list_x, list_t, list_p, list_n):
    x = list_x[0]
    y = list_x[1]

    tx = list_t[0]
    ty = list_t[1]

    px = list_p[0]
    py = list_p[1]

    nx = list_n[0]
    ny = list_n[1]

    Nx = np.zeros(px+1)
    Ny = np.zeros(py+1)

    Nx[:] = bsp.EvalBasisFuns(px,tx,x)
    Ny[:] = bsp.EvalBasisFuns(py,ty,y)

    return Nx,Ny
# ------------------------------

# ------------------------------
def getSites(Nx,Ny):
    # TODO : a changer pour le SL
    """
    generates the current position of our sites
    """
    x = np.linspace(0., 1., Nx)
    y = np.linspace(0., 1., Ny)
    X,Y = np.meshgrid(x,y)
    return X,Y
# ------------------------------

# ------------------------------
def assembleM(X, Y, list_t, list_p, list_n):
    """
    """
    n,m = X.shape
    nx = list_n[0]
    ny = list_n[1]

    px = list_p[0]
    py = list_p[1]

    # initialize the matrix M of size (nx ny, n m)
    M = np.zeros((n*m, nx*ny))

    for i in range(0,n):
        for j in range(0,m):
            I = j * n + i
            x = X[i,j] ; y = Y[i,j]
            Nx, Ny              = evalSplines([x,y], list_t, list_p, list_n)
            leftmkx, leftmky    = findSpan   ([x,y], list_t, list_p, list_n)
            for ip in range(0, px+1):
                for jp in range(0, py+1):
#                    J = (jp + leftmky) * nx + ip + leftmkx
                    J = (ip + leftmkx) * ny + jp + leftmky
                    M[I,J] = Nx[ip] * Ny[jp]

    return M
# ------------------------------

# ------------------------------
def getCoefC1Constraints(face, list_t, list_p, list_n):
    # ...
    # p : degree of the slave domain current curve
    # d : denominator for the slave domain current curve
    if face in [1,3]:
        axis = 1
    if face in [2,4]:
        axis = 0

    if face in [1,2]:
        p = list_p[axis]
        d = list_t[axis][p+1]
    if face in [3,4]:
        p = list_p[axis]
        d = list_t[axis][-1] - list_t[axis][-(p+1+1)]

    c = p / d
    return c
    # ...
# ------------------------------

# ------------------------------
def getC1Constraints(face, list_t, list_p, list_n):
    c = getCoefC1Constraints(face, list_t, list_p, list_n)
    nx = list_n[0]
    ny = list_n[1]
    list_v = []
    for i in range(0, nx):
        V = np.zeros((nx,ny))
        if face == 1:
            V[i,0] = -c
            V[i,1] =  c
            list_v.append(V.reshape(nx*ny))
        if face == 3:
            V[i,-1] = -c
            V[i,-2] =  c
            list_v.append(V.reshape(nx*ny))

    for j in range(0, ny):
        V = np.zeros((nx,ny))
        if face == 2:
            V[0,j] = -c
            V[1,j] =  c
            list_v.append(V.reshape(nx*ny))
        if face == 4:
            V[-1,j] = -c
            V[-2,j] =  c
            list_v.append(V.reshape(nx*ny))
    return list_v
# ------------------------------

#-----------------------------------
class surfint(object):
    def __init__(self, geometry, space=None, constraints=[]):
        """
        initialize the surfit object
        PDE is the Differential operator to use for smoothing (usually a 2nd
        order)
        constraints is a list of dictionaries that must be of the following form
        constraints[i] is  {'patch_id_m', 'face_m', 'patch_id_s', 'face_s',
        'type'}
        patch_id_m is the master patch id
        face_m     is the face id in the master patch
        patch_id_s is the slave  patch id
        face_s     is the face id in the slave  patch
        type       is the constraint's type: C1, C2, ... (default: C1)
        ib         is the starting index in the face element (default:0 )
        ie         is the ending   index in the face element (default:-1 )
        """
        self.geometry = geometry
        self.postAssembly = False
        self.nConstraints = 0
        self.ConstIndices = []
        self.ConstVals = []
        self.constraints = constraints

        # space may be None
        self.V = space

#    @property
#    def system(self):
#        return self.PDE.system.get()
#
#    @property
#    def space(self):
#        return self.PDE.space

    def AssembleLocalMatrix(self, nrb):
        list_t = nrb.knots
        list_p = nrb.degree
        list_n = nrb.shape

        nx = list_n[0]
        ny = list_n[1]

        n,m = list_n
        X,Y = getSites(n,m)

        from scipy.sparse import csr_matrix
        A = csr_matrix(assembleM(X, Y, list_t, list_p, list_n))

        return A

    def interpolate(self, f):
        nrb = self.geometry[0]
        list_t = nrb.knots
        list_p = nrb.degree
        list_n = nrb.shape

        nx = list_n[0]
        ny = list_n[1]

        n,m = list_n
        X,Y = getSites(n,m)

        M = self.AssembleLocalMatrix(nrb)

        list_faces = []
        list_v = []
        for face in list_faces:
            list_v += getC1Constraints(face, list_t, list_p, list_n)

        F = f(X,Y).reshape(n*m)

        from scipy.sparse import csr_matrix
        from scipy.sparse.linalg import splu

        # ...
        A = np.zeros((nx*ny,nx*ny))
        A[:n*m, :nx*ny] = M.todense()
        nConstraints = len(list_v)

        assert(nx*ny-n*m==nConstraints)
        for i,v in enumerate(list_v):
            A[n*m+i, :] = v
        # ...

        # ...
        B = np.zeros(nx*ny)
        B[:n*m] = F
        # ...

        A_ = csr_matrix(A)
        A_op = splu(A_.tocsc())
        y = A_op.solve(B)

        return y.reshape((nx,ny))
# ----------------------------------------------------------


if __name__ == '__main__':
    sin = np.sin; cos = np.cos ; pi = np.pi
#    f   = lambda x,y : sin ( 2*pi*x ) * sin ( 4*pi*y )
#    dxf = lambda x,y : 2*pi * cos ( 2*pi*x ) * sin ( 4*pi*y )
#    dyf = lambda x,y : 4*pi * sin ( 2*pi*x ) * cos ( 4*pi*y )

    f   = lambda x,y : 0.5 * (x**2 + y**2)
    dxf = lambda x,y : x
    dyf = lambda x,y : y

    from caid.cad_geometry import square as domain
    geo = domain(n=[31,31], p=[2,2])

    interpolator = surfint(geo)
    y   = interpolator.interpolate(f)


#    list_t = nrb.knots
#    list_p = nrb.degree
#    list_n = nrb.shape
#
#    nx = list_n[0]
#    ny = list_n[1]
#
#    n,m = list_n
##    n = nx-1 ; m = ny
#    X,Y = getSites(n,m)
#
#    M       = assembleM(X, Y, list_t, list_p, list_n)
#
##    list_faces = [1]
#    list_faces = []
#    list_v = []
#    for face in list_faces:
#        list_v += getC1Constraints(face, list_t, list_p, list_n)
#
#    F = f(X,Y).reshape(n*m)
#
#    from scipy.sparse import *
#    from scipy.sparse.linalg import spsolve
#
#    # ...
#    A = np.zeros((nx*ny,nx*ny))
#    A[:n*m, :nx*ny] = M
#    nConstraints = len(list_v)
#
#    print "shape M = ", M.shape
#    print "shape A = ", A.shape
#    print "nConstraints = ", nConstraints
#
#    print nx*ny-n*m
#    print nConstraints
#    assert(nx*ny-n*m==nConstraints)
#    for i,v in enumerate(list_v):
#        A[n*m+i, :] = v
#    # ...
#
#    # ...
#    b = dyf(X[:,0], Y[:,0])
#    print "b.shape = ", b.shape
#    print "expected = ", nConstraints
#    # ...
#
#    # ...
#    B = np.zeros(nx*ny)
#    B[:n*m] = F
#    try:
#        B[n*m:] = b # TODO a mettre les valeurs des derivees
#    except:
#        pass
#    # ...
#
#    A_ = csr_matrix(A)
#    y = spsolve(A_, B)

    import pylab as pl
    pl.contourf(y) ; pl.colorbar() ; pl.show()
