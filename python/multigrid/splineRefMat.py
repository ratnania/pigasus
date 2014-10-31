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
def refineSpline(dim, t, p, n, knots, P):
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
def constructCurveMatrix_single(t, p, n, knots):
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
def constructCurveMatrix(list_t, p, n, knots, DirFaces=[]):
    from scipy.sparse import csr_matrix
    if len(list_t) == 0:
        return csr_matrix(np.identity(n))

    r = list_t[0]
    M, knots_new= constructCurveMatrix_single(r, p, n, knots)
    knots = knots_new
    n += 1
    for r in list_t[1:]:
        M_tmp, knots_new = constructCurveMatrix_single(r, p, n, knots)
        knots = knots_new
        n = len(knots)-p-1
        M = np.dot(M_tmp, M)

    M = csr_matrix(M)

    # ... Dirichlet boundary condition
    ni, nj = M.shape
    ib = 0 ; ie = ni
    jb = 0 ; je = nj
    if 0 in DirFaces:
        ib += 1 ; jb += 1
    if 1 in DirFaces:
        ie -= 1 ; je -= 1

#    print "bc-i ", ib,ie, "     bc-j ",jb,je
    return csr_matrix(M.todense()[ib:ie,jb:je])
# ...

# ...
def constructSurfaceMatrix(list_r1, list_r2, p1, p2, n1, n2, u1, u2, DirFaces=[]):
    from scipy.sparse import csr_matrix, kron

    _DirFaces = []
    if 0 in DirFaces:
        _DirFaces.append(0)
    if 2 in DirFaces:
        _DirFaces.append(1)
    M1      = constructCurveMatrix(list_r1, p1, n1, u1, DirFaces=_DirFaces)

    _DirFaces = []
    if 1 in DirFaces:
        _DirFaces.append(0)
    if 3 in DirFaces:
        _DirFaces.append(1)
    M2      = constructCurveMatrix(list_r2, p2, n2, u2, DirFaces=_DirFaces)

    H       = csr_matrix(kron(M2,M1))

    return H, [n1+len(list_r1),n2+len(list_r2)]
# ...

if __name__ == '__main__':

    import caid.cad_geometry  as cg
    from caid.cad_geometry import line, bilinear

    # ...
    def geometry_1D(nx,px):
        # ...
        nrb = line(p0=(0,0), p1=(1,0))
        geo = cg.cad_geometry(geo=nrb)
        geo.refine(id=0,list_p=[px-1])
        tx = np.linspace(0.,1.,nx+2)[1:-1]
        geo.refine(id=0, list_t=[tx])
        # ...
        return geo
    # ...

    # ...
    def geometry_2D(nx,ny,px,py):
        # ...
    #    points = np.asarray([[[0.,0.],[0.,1.]],[[1.,0.],[1.,1.]]])
        points = np.asarray([[[0.,0.],[0.,1.]],[[2.,0.],[1.,1.]]])
        nrb = bilinear(points)
        geo = cg.cad_geometry(geo=nrb)
        geo.refine(id=0,list_p=[px-1, py-1])
        tx = np.linspace(0.,1.,nx+2)[1:-1]
        ty = np.linspace(0.,1.,ny+2)[1:-1]
        geo.refine(id=0, list_t=[tx, ty])
        # ...

        return geo
    # ...

    # ...
    def test1D1():
        list_r = list(np.random.random(20))
        for r in list_r:
            nx = 7
            px = 2
            geo = geometry_1D(nx, px)

            nrb     = geo[0]
            knots   = nrb.knots[0]
            n       = nrb.shape[0]
            p       = nrb.degree[0]
            P       = nrb.points
            dim     = P.shape[1]

            Q = refineSpline(dim, r, p, n, knots, P)
            M = constructCurveMatrix([r], p, n, knots)
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
    #    list_r = list(np.random.random(20))
        list_r = [0.1,0.2,0.3]

        nx = 3
        px = 2
        geo = geometry_1D(nx, px)

        nrb     = geo[0]
        knots   = nrb.knots[0]
        n       = nrb.shape[0]
        p       = nrb.degree[0]
        P       = nrb.points

        M = constructCurveMatrix(list_r, p, n, knots)
        from scipy.io import mmwrite
        mmwrite('M.mtx', M)
        R = M.dot(nrb.points[:,0])

        geo = geometry_1D(nx, px)
        geo.refine(id=0, list_t=[list_r])
        nrb     = geo[0]
        P = np.asarray(nrb.points[:,0])

        assert(np.allclose(P,R))
        print("test1D2: OK")
    # ...

    # ...
    def test2D1():
        list_r1 = list(np.random.random(20))
        list_r2 = list(np.random.random(20))

        nx = 10 ; ny = 15
        px = 3 ; py = 2
        geo = geometry_2D(nx, ny, px, py)

        dim     = geo.dim
        nrb     = geo[0]

        u1,u2   = nrb.knots
        n1,n2   = nrb.shape
        p1,p2   = nrb.degree

        M1      = constructCurveMatrix(list_r1, p1, n1, u1)
        M2      = constructCurveMatrix(list_r2, p2, n2, u2)
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
        list_r1 = list(np.random.random(20))
        list_r2 = list(np.random.random(20))
    #    list_r1 = [0.1, 0.2]
    #    list_r2 = [0.9]

        nx = 20 ; ny = 31
        px = 3 ; py = 2
        geo = geometry_2D(nx, ny, px, py)

        dim     = geo.dim
        nrb     = geo[0]

        u1,u2   = nrb.knots
        n1,n2   = nrb.shape
        p1,p2   = nrb.degree

        H,[n,m]  = constructSurfaceMatrix(list_r1, list_r2, p1, p2, n1, n2, u1, u2)

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
    def test2D3(DirFaces):
        list_r1 = list(np.random.random(20))
        list_r2 = list(np.random.random(20))
    #    list_r1 = [0.1, 0.2]
    #    list_r2 = [0.9]

        nx = 20 ; ny = 31
        px = 3 ; py = 2
        geo = geometry_2D(nx, ny, px, py)

        dim     = geo.dim
        nrb     = geo[0]

        u1,u2   = nrb.knots
        n1,n2   = nrb.shape
        p1,p2   = nrb.degree

        print("DirFaces ", DirFaces)
        H,[n,m]  = constructSurfaceMatrix(list_r1, list_r2, p1, p2, n1, n2, u1, u2, DirFaces=DirFaces)
        print("shape ", H.shape)

        Px      = nrb.points[:,:,0].copy()
        Py      = nrb.points[:,:,1].copy()

        print("test2D3: OK")
    # ...

    test1D1()
    test1D2()
    test2D1()
    test2D2()
    test2D3([0,1,2,3])
