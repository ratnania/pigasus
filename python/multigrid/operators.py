# -*- coding: UTF-8 -*-
import numpy as np

__all__ = ['constructCurveMatrix, constructSurfaceMatrix']


# ...
def alphaFct(i, k, t, n, p, knots):
    if   (i <= k-p):
        return 1.
    elif (k-p < i) and ( i <= k):
        return (t - knots[i]) / (knots[i+p] - knots[i])
    else:
        return 0.
# ...

# ...
def interv(knots, t):
    k = -1
    if t < knots[0]:
        return k

    if t > knots[-1]:
        return k

    # now knots[0] < t < knots[-1]
    for u in knots:
        if u>t:
            return k
        k += 1
# ...

# ...
def compute_H(geo_H, geo_h):
    dim = geo_H.dim
    npatchs = geo_H.npatchs

    list_H = []
    list_h = []

    for i in range(0, npatchs):
        nrb_H = geo_H[i]
        nrb_h = geo_h[i]
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
# ...

from pigasus.fem.pigasusObject import *
def operator(pigasusObject):
    """
    the generic implementation of an operator. used for both interpolation and
    reduction
    """
    def __init__(self, geometry_H, geometry_h, useDecouple=False):
        self.spl = splineRefMat(self.dim, useDecouple=useDecouple)

        self.useDecouple    = useDecouple
        self.dim            = geometry_H[0].dim
        self.npatchs        = geometry_H.nptachs
        self.geometry_H     = geometry_H
        self.geometry_h     = geometry_h

        # in list_csr_M we will store either the 1D-matrices (useDecouple)
        # or the fulle matrix
        self.list_csr_M = []

        from pigasus.fem.matrix import matrix
        if self.useDecouple:
            # ... create pigasus matrices
            self.list_M = []
            for i in range(0,self.npatchs):
                list_M = []
                for d in range(0, self.dim):
                    M = matrix()
                    list_M.append(M)
                self.list_M.append(list_M)
        else:
            # ... there is one matrix
            self.M = matrix()

    def construct(self):
        if self.useDecouple:
            i = 0
            for list_M in self.list_M:
                list_csr_M = self.list_csr_M[i]
                j = 0
                for M in list_M:
                    csr_M = list_csr_M[j]
                    M.set(csr_M) ; j += 1
                i += 1
        else:
            self.M.set(list_csr_M[0])

    def apply(self, list_F):
        """
        apply the operator on list_F a list of numpy d-arrays if useDecouple
        or a 1d-numpy array
        """
        operator.apply(self, list_F)

        if self.useDecouple:
            list_U = []

            if self.dim ==1:
                for i, F in enumerate(list_F):
                    list_M = self.list_M[i]
                    M = list_M[0]
                    U = M.dot(F)
                    list_U.append(U)

            if self.dim ==2:
                for i, F in enumerate(list_F):
                    list_M = self.list_M[i]
                    M1 = list_M[0] ; M2 = list_M[1]
                    a = F.data ; ja = F.indices ; ia = F.indptr
                    # M1 * X * M2
                    U = self.com.pyfem.compute_mxmxm(M1.id, a, ja, ia, M2.id)
                    list_U.append(U)

            return list_U
        else:
            return self.M.dot(list_F)

def interpolation(operator):
    def __init__(self, geometry_H, geometry_h):
        operator.__init__(self, geometry_H, geometry_h)

        for i in range(0, self.npatchs):
            nrb_H = self.geometry_H[i]
            nrb_h = self.geometry_h[i]

            if self.dim ==1:
                knots_H = nrb_H.knots[0]
                knots_h = nrb_h.knots[0]

                n = nrb_H.shape[0]
                p = nrb_H.degree[0]

                list_r = [r for r in knots_h if r not in knots_H]

                M = self.spl.construct(list_r, p, n, knots_H)

                list_csr_M = [M]

            if self.dim ==2:
                u_H1,u_H2   = nrb_H.knots
                n_H1,n_H2   = nrb_H.shape
                p_H1,p_H2   = nrb_H.degree

                u_h1,u_h2   = nrb_h.knots

                list_r1 = [r for r in u_h1 if r not in u_H1]
                list_r2 = [r for r in u_h2 if r not in u_H2]

                M1, M2 = self.spl.construct( list_r1, list_r2 \
                                    , p_H1, p_H2 \
                                    , n_H1, n_H2 \
                                    , u_H1, u_H2 )
                # transpose M2
                M2 = M2.transpose()
                list_csr_M = [M1,M2]

            self.list_csr_M.append(list_csr_M)

def reduction(operator):
    def __init__(self, geometry_H, geometry_h):
        """
        P : interpolation
        """
        operator.__init__(self, geometry_H, geometry_h)

        list_H, list_h = compute_H(geo_H, geo_h)
        for i in range(0, self.npatchs):
            nrb_H = self.geometry_H[i]
            nrb_h = self.geometry_h[i]

            if self.dim ==1:
                knots_H = nrb_H.knots[0]
                knots_h = nrb_h.knots[0]

                n = nrb_H.shape[0]
                p = nrb_H.degree[0]

                list_r = [r for r in knots_h if r not in knots_H]

                M = self.spl.construct(list_r, p, n, knots_H)

                M = M.transpose().tocsr()
                r = list_h[i]/list_H[i]
                M *= r

                list_csr_M = [M]

            if self.dim ==2:
                u_H1,u_H2   = nrb_H.knots
                n_H1,n_H2   = nrb_H.shape
                p_H1,p_H2   = nrb_H.degree

                u_h1,u_h2   = nrb_h.knots

                list_r1 = [r for r in u_h1 if r not in u_H1]
                list_r2 = [r for r in u_h2 if r not in u_H2]

                M1, M2 = self.spl.construct( list_r1, list_r2 \
                                    , p_H1, p_H2 \
                                    , n_H1, n_H2 \
                                    , u_H1, u_H2 )

                M1 = M1.transpose().tocsr()
                r = list_h[i]/list_H[i]
                M1 *= r

                M2 = M2.transpose().tocsr()
                r = list_h[i]/list_H[i]
                M2 *= r

                # transpose M2
                M2 = M2.transpose()
                list_csr_M = [M1,M2]

            self.list_csr_M.append(list_csr_M)

def coarse_matrix(pigasusObject):
    """

    """
    def __init__(self, geometry_H, geometry_h, DirFaces=[]):
        self.spl = splineRefMat(self.dim, useDecouple=False, DirFaces=DirFaces)

        self.dim        = geometry_H[0].dim
        self.npatchs    = geometry_H.nptachs
        self.geometry_H = geometry_H
        self.geometry_h = geometry_h
        self.DirFaces   = DirFaces

        if self.npatchs > 1:
            print("Multipatchs not yet implemented. STOP!")
            import sys; sys.exit(0)

#        for i in range(0, self.npatchs):
        for i in range(0, 1):
            nrb_H = self.geometry_H[i]
            nrb_h = self.geometry_h[i]

            if self.dim ==1:
                knots_H = nrb_H.knots[0]
                knots_h = nrb_h.knots[0]

                n = nrb_H.shape[0]
                p = nrb_H.degree[0]

                list_r = [r for r in knots_h if r not in knots_H]

                P = self.spl.construct(list_r, p, n, knots_H)

            if self.dim ==2:
                u_H1,u_H2   = nrb_H.knots
                n_H1,n_H2   = nrb_H.shape
                p_H1,p_H2   = nrb_H.degree

                u_h1,u_h2   = nrb_h.knots

                list_r1 = [r for r in u_h1 if r not in u_H1]
                list_r2 = [r for r in u_h2 if r not in u_H2]

                P = self.spl.construct( list_r1, list_r2 \
                                    , p_H1, p_H2 \
                                    , n_H1, n_H2 \
                                    , u_H1, u_H2 )

            # ... restriction
            list_H, list_h = compute_H(geo_H, geo_h)
            r = list_h[i]/list_H[i]
            R = P.transpose().tocsr()
            R *= r

            self.P = P
            self.R = R

    def construct(self, A_h):
        self.A = self.R * A_h * self.P

#-----------------------------------

class splineRefMat(object):
    def __init__(self, dim, useDecouple=False, DirFaces=[]):
        self.dim = dim
        self.useDecouple = useDecouple
        self.DirFaces = DirFaces
        self._DirFaces = []

    def set_DirFaces(self, DirFaces):
        self._DirFaces = DirFaces

    def construct(self, *args, **kwargs):
        if self.dim ==1:
            self.set_DirFaces(self.DirFaces)
            M = self.constructCurveMatrix(*args, **kwargs)
            return M

        if self.dim ==2:
            M = self.constructSurfaceMatrix(*args, **kwargs)
            return M

        if self.dim ==3:
            print("Not done yet")
            import sys; sys.exit(1)

    # ...
    def refineSpline(self, dim, t, p, n, knots, P):
        N       = n+1
        Q       = np.zeros((N,dim))

        # WE MUST LOCATE x WITH RESPECT TO THE KNOT VECTOR
        k = interv(knots, t)

        # UPDATE Q
        # WE USE THE FACT THAT [B]_j = w_j*alpha_j/w B_j + w_[j-1]*[1-alpha_j]/w B_[j-1]
        j = 0
        alpha = alphaFct(j, k, t, n, p, knots)
        Q[j,:] = alpha * P[j,:]

        for j in range(1,n):
            alpha = alphaFct(j, k, t, n, p, knots)
            Q[j,:] = alpha * P[j,:] + (1. - alpha) * P[j-1,:]

        j = n
        alpha = alphaFct(j, k, t, n, p, knots)
        Q[j,:] = (1. - alpha) * P[j-1,:]

        return Q
    # ...

    # ...
    def constructCurveMatrix_single(self, t, p, n, knots):
        N       = n+1
        M       = np.zeros((N,n))

        # WE MUST LOCATE x WITH RESPECT TO THE KNOT VECTOR
        k = interv(knots, t)

        knots_new = np.asarray(list(knots)[:k+1] + [t] + list(knots)[k+1:])

        # UPDATE Q
        # WE USE THE FACT THAT [B]_j = w_j*alpha_j/w B_j + w_[j-1]*[1-alpha_j]/w B_[j-1]
        j = 0
        alpha = alphaFct(j, k, t, n, p, knots)
        M[0,0] = alpha

        for j in range(1,n):
            alpha = alphaFct(j, k, t, n, p, knots)
            M[j,j]      = alpha
            M[j,j-1]    = 1. - alpha

        j = n
        alpha = alphaFct(j, k, t, n, p, knots)
        M[n,n-1] = 1. - alpha

        return M, knots_new
    # ...

    # ...
    def constructCurveMatrix(self, list_t, p, n, knots):
        DirFaces = self._DirFaces

        from scipy.sparse import csr_matrix
        if len(list_t) == 0:
            return csr_matrix(np.identity(n))

        r = list_t[0]
        M, knots_new= self.constructCurveMatrix_single(r, p, n, knots)
        knots = knots_new
        n += 1
        for r in list_t[1:]:
            M_tmp, knots_new = self.constructCurveMatrix_single(r, p, n, knots)
            knots = knots_new
            n = len(knots)-p-1
            M = np.dot(M_tmp, M)

        M = csr_matrix(M)

        # ... Dirichlet boundary condition
        n,m = M.shape
        ib = 0 ; ie = n
        jb = 0 ; je = m
        if 0 in DirFaces:
            ib += 1 ; jb += 1
        if 1 in DirFaces:
            ie -= 1 ; je -= 1

        return csr_matrix(M.todense()[ib:ie,jb:je])
    # ...

    # ...
    def constructSurfaceMatrix(self, list_r1, list_r2, p1, p2, n1, n2, u1, u2):
        from scipy.sparse import csr_matrix, kron

        DirFaces1 = []
        if 0 in DirFaces:
            DirFaces1.append(0)
        if 2 in DirFaces:
            DirFaces1.append(1)

        DirFaces2 = []
        if 1 in DirFaces:
            DirFaces2.append(0)
        if 3 in DirFaces:
            DirFaces2.append(1)

        self.set_DirFaces(DirFaces1)
        M1      = self.constructCurveMatrix(list_r1, p1, n1, u1)

        self.set_DirFaces(DirFaces2)
        M2      = self.constructCurveMatrix(list_r2, p2, n2, u2)

        if self.useDecouple:
            return M1, M2
        else:
            return csr_matrix(kron(M2,M1))
    # ...

if __name__ == '__main__':

    import caid.cad_geometry  as cg
    from caid.cad_geometry import line, square

    DIM_1D = 1
    DIM_2D = 2

    # ...
    def test1D1():
        spl = splineRefMat(DIM_1D)
        list_r = list(np.random.random(20))
        for r in list_r:
            nx = 7
            px = 2
            geo = line(n=[nx], p=[px])

            nrb     = geo[0]
            knots   = nrb.knots[0]
            n       = nrb.shape[0]
            p       = nrb.degree[0]
            P       = nrb.points
            dim     = P.shape[1]

            Q = spl.refineSpline(dim, r, p, n, knots, P)
            M = spl.construct([r], p, n, knots)
            R = M.dot(nrb.points[:,0])

            geo.refine(id=0, list_t=[r])
            nrb     = geo[0]
            #print nrb.knots[0]
            Q = np.asarray(Q[:,0])
            P = np.asarray(nrb.points[:,0])
            assert(np.allclose(P,Q))
            assert(np.allclose(P,R))

        print("test1D1: OK")
    # ...

    # ...
    def test1D2():
        spl = splineRefMat(DIM_1D)
    #    list_r = list(np.random.random(20))
        list_r = [0.1,0.2,0.3]

        nx = 3
        px = 2
        geo = line(n=[nx], p=[px])

        nrb     = geo[0]
        knots   = nrb.knots[0]
        n       = nrb.shape[0]
        p       = nrb.degree[0]
        P       = nrb.points

        M = spl.construct(list_r, p, n, knots)
        from scipy.io import mmwrite
        mmwrite('M.mtx', M)
        R = M.dot(nrb.points[:,0])

        geo = line(n=[nx], p=[px])
        geo.refine(id=0, list_t=[list_r])
        nrb     = geo[0]
        P = np.asarray(nrb.points[:,0])

        assert(np.allclose(P,R))
        print("test1D2: OK")
    # ...

    # ...
    def test2D1():
        spl = splineRefMat(DIM_1D)
        list_r1 = list(np.random.random(20))
        list_r2 = list(np.random.random(20))

        nx = 10 ; ny = 15
        px = 3 ; py = 2
        geo = square(n=[nx, ny], p=[px, py])

        dim     = geo.dim
        nrb     = geo[0]

        u1,u2   = nrb.knots
        n1,n2   = nrb.shape
        p1,p2   = nrb.degree

        M1      = spl.construct(list_r1, p1, n1, u1)
        M2      = spl.construct(list_r2, p2, n2, u2)
        tM2     = M2.transpose().tocsr()

        Px      = nrb.points[:,:,0].copy()
        Py      = nrb.points[:,:,1].copy()

        geo.refine(id=0, list_t=[list_r1, list_r2])
        nrb     = geo[0]
        Qx      = np.asarray(nrb.points[:,:,0])
        Qy      = np.asarray(nrb.points[:,:,1])

        from scipy.sparse import csr_matrix

        list_P = [Px, Py]
        list_Q = [Qx, Qy]
    #    list_P = [Px]
    #    list_Q = [Qx]
        for (U,Q) in zip(list_P, list_Q):
            Us  = csr_matrix(U).dot(tM2)
            tV  = M1.dot(Us).todense()

            assert(np.allclose(tV, Q))

        print("test2D1: OK")
    # ...

    # ...
    def test2D2():
        spl = splineRefMat(DIM_2D)
        list_r1 = list(np.random.random(20))
        list_r2 = list(np.random.random(20))
    #    list_r1 = [0.1, 0.2]
    #    list_r2 = [0.9]

        nx = 20 ; ny = 31
        px = 3 ; py = 2
        geo = square(n=[nx, ny], p=[px, py])

        n = nx + px + 1 + len(list_r1)
        m = ny + py + 1 + len(list_r2)

        dim     = geo.dim
        nrb     = geo[0]

        u1,u2   = nrb.knots
        n1,n2   = nrb.shape
        p1,p2   = nrb.degree

        H = spl.construct(list_r1, list_r2, p1, p2, n1, n2, u1, u2)

        Px      = nrb.points[:,:,0].copy()
        Py      = nrb.points[:,:,1].copy()

        geo.refine(id=0, list_t=[list_r1, list_r2])
        nrb     = geo[0]
        Qx      = np.asarray(nrb.points[:,:,0])
        Qy      = np.asarray(nrb.points[:,:,1])

        list_P = [Px, Py]
        list_Q = [Qx, Qy]
    #    list_P = [Px]
    #    list_Q = [Qx]
        for (U,Q) in zip(list_P, list_Q):
            nU,mU = U.shape
            vecU    = U.transpose().reshape(nU*mU)
            vecP    = H.dot(vecU)
            P       = vecP.reshape((m,n)).transpose()

            assert(np.allclose(P, Q))

        print("test2D2: OK")
    # ...

    # ...
    def test2D3():
        spl = splineRefMat(DIM_2D, useDecouple=True)
        list_r1 = list(np.random.random(20))
        list_r2 = list(np.random.random(20))
    #    list_r1 = [0.1, 0.2]
    #    list_r2 = [0.9]

        nx = 20 ; ny = 31
        px = 3 ; py = 2
        geo = square(n=[nx, ny], p=[px, py])

        n = nx + px + 1 + len(list_r1)
        m = ny + py + 1 + len(list_r2)

        dim     = geo.dim
        nrb     = geo[0]

        u1,u2   = nrb.knots
        n1,n2   = nrb.shape
        p1,p2   = nrb.degree

        H1, H2 = spl.construct(list_r1, list_r2, p1, p2, n1, n2, u1, u2)

        assert(np.allclose(np.array(H1.shape), np.array((44,24))))
        assert(np.allclose(np.array(H2.shape), np.array((54,34))))

        print("test2D3: OK")
    # ...

    test1D1()
    test1D2()
    test2D1()
    test2D2()
    test2D3()
