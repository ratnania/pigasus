# coding: utf-8
from caid.cad_geometry import linear
from caid.cad_geometry import cad_geometry, cad_nurbs
import numpy as np
from numpy import array, asarray
from igakit.cad import join
from caid.utils.intersection import intersect_crv
import sys
import os
from pigasus.utils.utils import progress
from caid.graphics.color import Color
from glob import glob
from scipy.interpolate import splprep
# ...
try:
    import Pycluster
    CLUSTER = "Pycluster"
except ImportError:
    print("Pycluster could not be found. scipy.cluster.vq.kmeans2 will be used, even if it is deprecated. ")
    from scipy.cluster.vq import *
    CLUSTER = "scipy"
# ...

# ------------------------------------------------------
def points_to_geo_interp(x, y, z=None, p=3):
    if z is None:
        z = np.zeros_like(x)

	tck, u = splprep([x,y], s=0, k=p)
    #-----------------------------------
    geo = cad_geometry()
    knots = [tck[0]]
    Px = tck[1][0]
    Py = tck[1][1]
    P = np.zeros((len(Px),3))
    P[:,0] = Px
    P[:,1] = Py

    nrb = cad_nurbs(knots, P)
    geo.append(nrb)
    #-----------------------------------

    return geo
# ------------------------------------------------------

# ------------------------------------------------------
def points_to_geo_approx(x, y, z=None, n=15, p=3, alpha=1.):
    if z is None:
        z = np.zeros_like(x)


    #-----------------------------------
    def MakeConstraint(cond, face=None, value=None):
        if cond.lower() == "closed":
            constraint = {}
            constraint['patch_id_m'] = 0
            constraint['face_m']     = 0
            constraint['patch_id_s'] = 0
            constraint['face_s']     = 1
        if cond.lower() == "c0":
            constraint = {}
            constraint['patch_id_m'] = 0
            constraint['face_m']     = face
            constraint['type']       = "C0"
            constraint['values']     = [value]
        if cond.lower() == "c1":
            constraint = {}
            constraint['patch_id_m'] = 0
            constraint['face_m']     = face
            constraint['type']       = "C1"
            constraint['values']     = [value]
        return constraint
    #-----------------------------------

    from pigasus.fit.curfit import curfit
    from pigasus.fit.curfit import compute_uk

    #-----------------------------------
    #method = "uniform"
    method = "chord"
    #method = "centripetal"
    #-----------------------------------

    #-----------------------------------
    # ...
    list_x = list(x) ; list_y = list(y)
    # ...

    # ...
    list_Q = list(zip(list_x, list_y))
    uk = compute_uk(list_Q, method=method)
    U1       = []
    U1      += list(uk)
    list_xk  = []     ; list_yk  = []
    list_xk += list_x ; list_yk += list_y

    lists_uk = [U1]
    lists_xk = [list_xk]
    lists_yk = [list_yk]
    # ...
    #-----------------------------------

    #-----------------------------------
    constraints = []

    # ... C0_COND
    constraint = MakeConstraint("C0", 0, [x[0], y[0]])
    constraints.append(constraint)

    constraint = MakeConstraint("C0", 1, [x[-1], y[-1]])
    constraints.append(constraint)
    #-----------------------------------

    #-----------------------------------
    from caid.cad_geometry import line
    geo = line(n=[n], p=[p])
    #-----------------------------------

    #-----------------------------------
    fit = curfit(geometry=geo, constraints=constraints, alpha=alpha)
    #-----------------------------------

    #-----------------------------------
    patch_id = 0
    xk = lists_xk[patch_id]
    yk = lists_yk[patch_id]

    geo = fit.construct([xk, yk], uk=lists_uk)
    #-----------------------------------

    return geo
# ------------------------------------------------------

# ------------------------------------------------------
def create_wall(R_wall, Z_wall):
    zz = np.zeros_like(R_wall)
    list_P = list(zip(R_wall, Z_wall, zz))
    list_crv = []
    for P,Q in zip(list_P[:-1], list_P[1:]):
        points = np.zeros((2,2))
        points[0,0] = P[0] ; points[0,1] = P[1]
        points[1,0] = Q[0] ; points[1,1] = Q[1]
        crv = linear(points = points)[0]
        list_crv.append(crv)

    nrb = list_crv[0]
    axis = 0
    for crv in list_crv[1:]:
        nrb = join(nrb, crv, axis)
        nrb = cad_nurbs(nrb.knots, nrb.points, weights=nrb.weights)

    geo = cad_geometry()
    geo.append(nrb)
    return geo
# ------------------------------------------------------

# ------------------------------------------------------
def multilevel_intersection(  c0, c1 \
                            , nlevel=4, rtol=1e-03, atol=1e-04 \
                            , n_notfound=2, verbose=False):
    npts = 4
    ll_condition = True
    level = 0
    i_notfound = 0

    list_P, list_t, list_s, ierr = intersect_crv(c0, c1, npts=npts)
    if len(list_t) > 0:
        t_old = np.asarray(list_t)
        s_old = np.asarray(list_s)
        level += 1
    else:
        i_notfound += 1
    npts *= npts

    while ll_condition and (level < nlevel) and (i_notfound < n_notfound):
        list_P, list_t, list_s, ierr = intersect_crv(c0, c1, npts=npts)
        if len(list_t) > 0:
            if verbose:
                print(">>> level : ", level)
            t_new = np.asarray(list_t)
            s_new = np.asarray(list_s)
            if len(t_new) == len(t_old):
                ll_condition = not np.allclose(t_new, t_old, rtol=rtol, atol=atol)
                ll_condition = ll_condition and \
                        not np.allclose(s_new, s_old, rtol=rtol, atol=atol)
            else:
                ll_condition = True

            t_old = t_new
            s_old = s_new

            level += 1
        else:
            if verbose:
                print("--- notfound : ", i_notfound)
            i_notfound += 1

        npts *= npts

    return list_P, list_t, list_s, ierr
# ------------------------------------------------------

# ------------------------------------------------------
class fluxAtWall():
    def __init__(  self, crv_wall, crv_flux\
                 , nlevel=4 \
                 , rtol=1e-03 \
                 , atol=1e-04 \
                 , npts=None \
                 , geo_flux=None \
                 , find_intersection=True):

        self._crv_wall = crv_wall
        self._crv_flux = crv_flux
        self._geo_flux = geo_flux

        if find_intersection:
            if npts is None:
                list_P, list_t, list_s, err = multilevel_intersection( crv_wall, crv_flux \
                                                                     , nlevel=nlevel \
                                                                     , rtol=rtol \
                                                                     , atol=atol)
            else:
                list_P, list_t, list_s, ierr = intersect_crv(crv_wall, crv_flux \
                                                             , npts=npts)

            self._nt     = len(list_t)
            self._list_P = list_P
            self._list_t = list_t
            self._list_s = list_s
        else:
            self._nt     = None
            self._list_P = []
            self._list_t = []
            self._list_s = []

        self._distances = {}
        self._means = {}
        self._geo_int = None

    @property
    def wall(self):
        return self._crv_wall

    @property
    def flux(self):
        return self._crv_flux

    @property
    def geo_intersection(self):
        return self._geo_int

    @property
    def patch_id(self):
        try:
            return self._geo_flux._list.index(self._crv_flux)
        except:
            return None

    @property
    def n_intersect(self):
        return self._nt

    @property
    def found_intersection(self):
        return (self._list_t is not None) \
                and (self._nt > 0)

    @property
    def list_t(self):
        return self._list_t

    @property
    def list_s(self):
        return self._list_s

    @property
    def points(self):
        return self._list_P

    @property
    def distances(self):
        return self._distances

    @property
    def means(self):
        return self._means

    @property
    def cost(self):
        if self._nt == 2:
            return self.distances[0,1]
        if self._nt == 4:
            M_opt = None
            i_opt = None
            j_opt = None

            for i in range(0, self._nt):
                for j in range(0, i):
                    M = self.means[j,i]
                    if M_opt is None:
                        M_opt = M
                        i_opt = i
                        j_opt = j
                    elif M[1] > M_opt[1]:
#                    elif M[1] < M_opt[1]:
                        M_opt = M
                        i_opt = i
                        j_opt = j
            if M_opt is None:
                print("No optimal point has been found.")
                print(M_opt, j_opt, i_opt)
                raise
            return self.distances[j_opt,i_opt]
#        if self._nt == 4:
#            M_min = None
#            i_min = None
#            j_min = None
#
#            M_max = None
#            i_max = None
#            j_max = None
#
#            for i in range(0, self._nt):
#                for j in range(0, i):
#                    M = self.means[j,i]
#                    if M_min is None:
#                        M_min = M
#                        i_min = i
#                        j_min = j
#                    elif M[1] < M_min[1]:
#                        M_min = M
#                        i_min = i
#                        j_min = j
#
#                    if M_max is None:
#                        M_max = M
#                        i_max = i
#                        j_max = j
#                    elif M[1] > M_max[1]:
#                        M_max = M
#                        i_max = i
#                        j_max = j
#
#            d_min = self.distances[j_min,i_min]
#            d_max = self.distances[j_max,i_max]
##            return d_min + 1./d_max
#            return np.log(d_min) - np.log(d_max)

    @property
    def center(self):
        if self._nt == 2:
            return self.means[0,1]
        if self._nt == 4:
            M_opt = None
            i_opt = None
            j_opt = None

            for i in range(0, self._nt):
                for j in range(0, i):
                    M = self.means[j,i]
                    if M_opt is None:
                        M_opt = M
                        i_opt = i
                        j_opt = j
                    elif M[1] < M_opt[1]:
                        M_opt = M
                        i_opt = i
                        j_opt = j

            return M_opt

    def compute_distance(self, eps=None):
        """
        computes the distance between all points (that are intersection of the
        wall and the flux surfaces).
        if eps is set, we only consider distances greater than eps
        """
        if not self.found_intersection:
            return None

        if self._nt < 2:
            return None

        if eps is None:
            eps = 0.0

        distances = {}
        for i in range(0, self._nt):
            A = self._list_P[i]
            for j in range(0, i):
                B = self._list_P[j]

                my_distance = (A[0]-B[0])**2 + (A[1]-B[1])**2 + (A[2]-B[2])**2
                if my_distance < eps:
                    my_distance = None
                distances[j,i] = my_distance

        self._distances = distances
        return distances

    def compute_means(self):
        """
        computes the mean all points (that are intersection of the
        wall and the flux surfaces).
        """
        if self._nt < 2:
            return None

        means = {}
        for i in range(0, self._nt):
            A = self._list_P[i]
            for j in range(0, i):
                B = self._list_P[j]

                mean = [.5*(a+b) for (a,b) in zip(A,B)]
                means[j,i] = asarray(mean)

        self._means = means
        return means

    def intersect(self):
        if not self.found_intersection:
            self._geo_int = None
            return None

        geo_out = cad_geometry()

        list_t = asarray(self.list_t); list_t.sort()
        list_s = asarray(self.list_s); list_s.sort()

        axis = 0

        geo = cad_geometry()
        geo.append(self.wall)
        for i,t in enumerate(list_t):
            geo.split(i,t,axis)
        for i_nrb in range(0, geo.npatchs):
            geo_out.append(geo[i_nrb])

        geo = cad_geometry()
        geo.append(self.flux)
        for i,s in enumerate(list_s):
            geo.split(i,s,axis)

        for i_nrb in range(0, geo.npatchs):
            geo_out.append(geo[i_nrb])

        self._geo_int = geo_out
        return geo_out

    def save(self, name, filename, color=Color("red").rgb):
        if self.geo_intersection is not None:
            geo = self.geo_intersection
            geo.set_attribut("name", name)
            geo.set_attribut('color', str(color))
            # write intersection points
            geo.set_attribut('n_intersection', self.n_intersect)
            for j in range(0, self.n_intersect):
                P_str = str(self.points[j][0])
                for x in self.points[j][1:]:
                    P_str += " , " + str(x)
                geo.set_attribut('intersection_P'+str(j), P_str)
                geo.set_attribut('intersection_t'+str(j), self.list_t[j])
                geo.set_attribut('intersection_s'+str(j), self.list_s[j])

            geo.save(filename)

    def load(self, filename):
        geo = cad_geometry(filename)
        self._geo_int = geo
        self._nt      = int(geo.get_attribut("n_intersection"))
        self._list_P = []
        for j in range(0, self.n_intersect):
            P_str = geo.get_attribut('intersection_P'+str(j))
            P = [float(x) for x in P_str.split(',')]

            t = geo.get_attribut('intersection_t'+str(j))
            s = geo.get_attribut('intersection_s'+str(j))

            P = asarray(P)
            t = float(t)
            s = float(s)

            self._list_P.append(P)
            self._list_t.append(t)
            self._list_s.append(s)

# ------------------------------------------------------

# ------------------------------------------------------
def clustering(x,y,cost,ngroup=2):
    if CLUSTER == "scipy":
        z = whiten(cost)

        # let scipy do its magic (k==3 groups)
        res, labels = kmeans2(array(list(zip(x,y,z))),ngroup)

    if CLUSTER == "Pycluster":
        points = np.zeros((x.shape[0], 2))
        points[:,0] = x
        points[:,1] = y

#        labels, error, nfound = Pycluster.kcluster(points, ngroup, weights=cost)
        labels, error, nfound = Pycluster.kcluster(points, ngroup)

    return labels
# ------------------------------------------------------

# ------------------------------------------------------
class container():
    def __init__(self, geo_flux, geo_wall):
        self._geo_flux      = geo_flux
        self._geo_wall      = geo_wall
        self._list_fw       = None
        self._n_xpoint      = None
        self._ngroup        = None
        self._dmin          = None
        self._dmax          = None
        self._list_geo_opt  = None

        axis = 0
        list_geo = []
        # including the wall
        crv_wall = geo_wall[0]

        list_id_flux = list(range(0, geo_flux.npatchs))

        # ...
        prog = progress(len(list_id_flux), "Creating Fluxes objects")
        list_fw = []
        for rank, id_patch in enumerate(list_id_flux):
            prog.log(rank)

            crv_flux = geo_flux[id_patch]

            fw = fluxAtWall(crv_wall, crv_flux \
                            , geo_flux=geo_flux \
                            , find_intersection=False)

            list_fw.append(fw)

        prog.finalize()
        # ...

        self.set_list_fw(list_fw)

    @property
    def geo_flux(self):
        return self._geo_flux

    @property
    def geo_wall(self):
        return self._geo_wall

    @property
    def list_fw(self):
        return self._list_fw

    @property
    def list_geo_opt(self):
        return self._list_geo_opt

    @property
    def n_xpoint(self):
        if self._n_xpoint is None:
            print("Warning n_xpoint is None")

        return self._n_xpoint

    @property
    def ngroup(self):
        if self._ngroup is None:
            print("Warning ngroup is None")
        return self._ngroup

    @property
    def dmin(self):
        if self._dmin is None:
            print("Warning dmin is None")
        return self._dmin

    @property
    def dmax(self):
        if self._dmax is None:
            print("Warning dmax is None")
        return self._dmax

    @property
    def centers(self):
        for fw in self.list_fw:
            if len(fw.means) == 0:
                fw.compute_means()
        return [fw.center for fw in self.list_fw]

    def set_list_fw(self, list_fw):
        self._list_fw = list_fw

    def gather_intersections(self, n_intersection=2 \
                    , ngroup=2 \
                    , dmin=False \
                    , dmax=True \
                    , verbose=False):

        if verbose:
            print(" ")
            print("gen_centers called with ngroup=" \
                    + str(ngroup) + " and n_intersection="+ str(n_intersection))

        if dmax and dmin:
            raise("dmax and dmin can not be True.")

        list_fw = self.list_fw

        _list_fw = [fw for fw in list_fw if fw.n_intersect == n_intersection]
        for fw in _list_fw:
            fw.compute_means()

        list_means = [(fw.center,fw.patch_id) \
                     for fw in _list_fw]

        list_M     = [fw.center for fw in _list_fw]
        list_cost  = [fw.cost for fw in _list_fw]
        x,y,z = list(zip(*list_M))
        x = asarray(x)
        y = asarray(y)
        cost = np.asarray(list_cost)

        list_group = clustering(x,y,cost,ngroup=ngroup)
        if verbose:
            print(" ")
            print("---------------")
            print("centers ", list_M)
            print("cost    ", cost)
            print("groups  ", list_group)

        # for each group we put the flux surface that minimize the intersection
        # distance
        list_fw_opt = []

        for i_group in range(0, ngroup):
            list_fw_group = [fw \
                             for (id_g, fw) in zip(list_group, _list_fw) \
                             if (id_g == i_group) and (fw.cost is not None)]

            list_dist = [fw.cost for fw in list_fw_group]
            if verbose:
                print("\n>>> i_group ", i_group)
                print("list_dist     ", list_dist)

            if dmin:
                d_opt = np.min(np.asarray(list_dist))
                i_opt = np.argmin(np.asarray(list_dist))
            if dmax:
                d_opt = np.max(np.asarray(list_dist))
                i_opt = np.argmax(np.asarray(list_dist))

            if verbose:
                print("d_opt ", d_opt)
                print("i_opt " , i_opt)

            fw_opt = list_fw_group[i_opt]
            list_fw_opt.append(fw_opt)

        return list_fw_opt

    def wall_intersection(self \
                          , npts=10 \
                          , nlevel=4 \
                          , rtol=1e-03 \
                          , atol=1e-04 \
                          , eps= 0.0 \
                          , verbose=False):
        geo_flux = self.geo_flux
        geo_wall = self.geo_wall

        axis = 0
        list_geo = []
        # including the wall
        crv_wall = geo_wall[0]

        list_id_flux = list(range(0, geo_flux.npatchs))

        # ...
        prog = progress(len(list_id_flux), "Computing intersections with the Wall")
        list_fw = []
        for rank, id_patch in enumerate(list_id_flux):
            prog.log(rank)

            crv_flux = geo_flux[id_patch]
            if verbose:
                print("\n<<< patch " + str(id_patch) + " >>>")

            fw = fluxAtWall(crv_wall, crv_flux \
                            , nlevel=nlevel \
                            , rtol=rtol \
                            , atol=atol \
                            , npts=npts \
                            , geo_flux=geo_flux)

            if verbose:
                print("n_intersection ", fw.n_intersect)
            fw.intersect()
            fw.compute_distance(eps=eps)

            list_fw.append(fw)

        prog.finalize()
        # ...

        self.set_list_fw(list_fw)

    def gather_fluxes(self, n_intersection \
                     , n_xpoint=None \
                     , list_ngroup=None \
                     , ngroup=None \
                     , dmin=None \
                     , dmax=None \
                     , reset=False):
        if reset:
            self._list_geo_opt = []

        if self._list_geo_opt is None:
            self._list_geo_opt = []

        if n_xpoint is None:
            n_xpoint = self.n_xpoint

        if (list_ngroup is None) and (ngroup is None):
            list_ngroup = [self.ngroup]
        elif (list_ngroup is None):
            list_ngroup = [ngroup]

        list_fw = self.list_fw

        # ...
        if list_ngroup is None:
            ngroup, _dmin, _dmax = self.info(n_intersection)
            list_ngroup = [ngroup]
        else:
            # ntmp will not be used as ngroup has been specified
            ntmp, _dmin, _dmax = self.info(n_intersection)

        if dmin is None:
            dmin = _dmin

        if dmax is None:
            dmax = _dmax

        txt  =  "Gathering fluxes with " \
                        +str(n_intersection)+ \
                        " intersections"
        txt  +=  " into "+str(list_ngroup)+" groups"
        prog_g = progress(len(list_ngroup), txt)
        for rank_g,ngroup in enumerate(list_ngroup):
            prog_g.log(rank_g)

            list_fw_min = self.gather_intersections(n_intersection=n_intersection \
                                          , ngroup=ngroup \
                                          , dmin=dmin \
                                          , dmax=dmax \
                                          , verbose=False)
            if list_fw_min is not None:
                for fw in list_fw_min:
                    self._list_geo_opt.append(fw.geo_intersection)

##            if list_fw_min is None:
#        try:
#            list_fw_min
#        except None as e:
#            print "Re-run with less groups"
#            print sys.exc_info()[0]
#            raise

        prog_g.finalize()
        # ...

        return self._list_geo_opt

    def save(self, dirname, color=Color("red").rgb):
        _dirname = dirname + "/intersection"
        os.system("mkdir -p " + _dirname)
        basename = _dirname + "/geo_flux"
        if self.list_fw is not None:
            for i,fw in enumerate(self.list_fw):
                name = "Intersection " + str(i)
                filename = basename + str(i) + ".xml"
                fw.save(name, filename, color=Color("red").rgb)

        _dirname = dirname + "/optimal"
        os.system("mkdir -p " + _dirname)
        basename = _dirname + "/geo_flux"
        if self.list_geo_opt is not None:
            for i in range(0, len(self.list_geo_opt)):
                geo = self.list_geo_opt[i]
                geo.set_attribut("name", "Optimal Fluxes "+str(i))
                geo.set_attribut('color', str(color))
                geo.save(basename+str(i)+".xml")

    def load(self, dirname, eps=1.e-2):
        _dirname = dirname + "/intersection"
        pattern = 'geo_flux*.xml'
        wildcard = os.path.join(_dirname, pattern)
        filenames = glob(wildcard)
        filenames.sort()
        pattern = dirname + "/intersection/geo_flux"
        for filename in filenames:
            i = int(filename.split(pattern)[-1].split(".")[0])

            fw = self.list_fw[i]
            fw.load(filename)
            fw.compute_distance(eps=eps)
            fw.compute_means()

# ------------------------------------------------------

# ------------------------------------------------------
class iter(container):
    def __init__(self, *args, **kwargs):
        container.__init__(self, *args, **kwargs)

        self._n_xpoint  = 1

    def info(self, n_intersection):
        if n_intersection == 2:
            self._dmin=False
            self._dmax=True
            self._ngroup = 2
        elif n_intersection == 4:
            self._dmin=True
            self._dmax=False
            self._ngroup = 1

        return self.ngroup, self.dmin, self.dmax
# ------------------------------------------------------

# ------------------------------------------------------
class west(container):
    def __init__(self, *args, **kwargs):
        container.__init__(self, *args, **kwargs)

        self._n_xpoint  = 2

    def info(self, n_intersection):
        if n_intersection == 2:
            self._dmin=False
            self._dmax=True
            self._ngroup = 2
        elif n_intersection == 4:
            self._dmin=True
            self._dmax=False
            self._ngroup = 1

        return self.ngroup, self.dmin, self.dmax
# ------------------------------------------------------


