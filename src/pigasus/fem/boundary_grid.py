# -*- coding: UTF-8 -*-
#! /usr/bin/python

# To change this template, choose Tools | Templates
# and open the template in the editor.


__author__="ARA"
__all__ = ['boundary_grid']
__date__ ="$Jan 13, 2012 11:48:35 AM$"

# -------------------------------------------
def get_nel(nrb):
    list_nel = []
    nel = 1
    li_dim = nrb.dim
    import numpy as np
    for li_d in range(0, li_dim):
        # we test if the knot vector is periodic
        li_n = nrb.shape[li_d]
        li_p = nrb.degree[li_d]
        li_un = len(np.unique(nrb.knots[li_d]))
        if li_un == li_n + li_p + 1:
            li_nel = li_n - li_p
        else:
            li_nel = li_un - 1
        list_nel.append(li_nel)
        nel *= li_nel
    return nel, list_nel
# -------------------------------------------

# -------------------------------------------
def get_bounds(nrb,ai_d):
    import numpy as np
    li_n = nrb.shape[ai_d]
    li_p = nrb.degree[ai_d]
    li_un = len(np.unique(nrb.knots[ai_d]))
    ll_isperiodic = (li_un == li_n + li_p + 1)
    if ll_isperiodic :
        list_u = np.array([nrb.knots[ai_d][i] for i in range(li_p,li_n+1) ])
    else:
        list_u = np.unique(nrb.knots[ai_d])
    tmin = min(list_u)
    tmax = max(list_u)
    return [tmin,tmax]
# -------------------------------------------

from quadratures import *
from pigasusObject import *
class boundary_grid(pigasusObject):
    def __init__(self,ai_current_grids, ai_id, nrb \
    , api_k=None, tensorlevel=1, as_type="lobatto", faces=None):
        pigasusObject.__init__(self)

        if api_k is None:
            print "boundary_grid : api_k must be given"
            import sys
            sys.exit(2)

        self.current_grids = ai_current_grids

        self.id = ai_id
        self.patch = nrb
        self.tensorlevel = tensorlevel
        self.dim = self.patch.dim
        self.Rd = self.patch.points.shape[-1]
        self.maxnpts = 0
        self.dirmaxnpts = 0

        if faces is None:
            self.faces = range(0, 2*self.dim)
        else:
            self.faces = faces

        li_nel, lpi_nel = get_nel(nrb)

        # in the case of the boundary grid, the nbr of elements is dim * lpi_nel[d]
        if self.dim==1:
            self.nel = 2
#        if self.dim==2:
#            self.nel = 2 * sum(lpi_nel)
#        if self.dim==3:
#            self.nel = 4 * sum(lpi_nel)
        if self.dim > 1:
            self.nel = 0
            for face in self.faces:
                axis = face % self.dim
                self.nel += lpi_nel[axis]

        if self.tensorlevel in [1]:
            self.gen_grid(api_k, as_type)
        else:
            print "boundary_grid : Not yet implemented"
            import sys
            sys.exit(2)

    def gen_vector(self, ai_d, ai_k, as_type):
        import numpy as np
        qd = quadratures()
        li_n = self.patch.shape[ai_d]
        li_p = self.patch.degree[ai_d]
        li_un = len(np.unique(self.patch.knots[ai_d]))
        ll_isperiodic = (li_un == li_n + li_p + 1)
        if ll_isperiodic :
            list_u = np.array([self.patch.knots[ai_d][i] for i in range(li_p,li_n+1) ])
            li_nel = li_n - li_p
        else:
            list_u = np.unique(self.patch.knots[ai_d])
            li_nel = len(list_u)-1

        li_maxnpts = ai_k + 1

        lpr_local = np.zeros((li_nel,li_maxnpts+1), dtype=np.double)
        lpr_localw = np.ones((li_nel,li_maxnpts), dtype=np.double)

        [lpr_x,lpr_w] = qd.generate(np.asarray(list_u), ai_k, as_type)
        for li_i in range(0,li_nel):
            lpr_local [li_i,0] = ai_k + 1
            lpr_local [li_i,1:li_maxnpts+1] = lpr_x[li_i,:]
            lpr_localw[li_i,:] = lpr_w[li_i,:]

#        for li_i in range(0,li_nel):
#            lr_a = list_u[li_i]
#            lr_b = list_u[li_i+1]
#            [lpr_x,lpr_w] = qd.generate(lr_a, lr_b, 2, ai_k, as_type)
#
#            lpr_local [li_i,0] = ai_k + 1
#            lpr_local [li_i,1:li_maxnpts+1] = lpr_x[0,:]
#            lpr_localw[li_i,:] = lpr_w[0,:]

        return li_nel, lpr_local, lpr_localw

    def gen_grid(self, api_k, as_type):
        if self.dim==1:
            self.gen_grid_1D(api_k, as_type)
        if self.dim==2:
            self.gen_grid_2D(api_k, as_type)
        if self.dim==3:
            self.gen_grid_3D(api_k, as_type)
        if self.dim not in [1,2,3]:
            print "gen_grid : dimension not done yet"

    def save_grid(self):
        if self.dim==1:
            self.save_grid_1D()
        if self.dim==2:
            self.save_grid_2D()
        if self.dim==3:
            self.save_grid_3D()
        if self.dim not in [1,2,3]:
            print "save_grid : dimension not done yet"

    def gen_grid_1D(self, ai_k, as_type):
        print "boundary_grid : gen_grid_1D not done yet"

    def gen_grid_2D(self, api_k, as_type):

        # x-direction
        li_d = 0; li_k = api_k[li_d]
        li_nx, lpr_localx, lpr_localwx = self.gen_vector(li_d, li_k, as_type)
        # y-direction
        li_d = 1; li_k = api_k[li_d]
        li_ny, lpr_localy, lpr_localwy = self.gen_vector(li_d, li_k, as_type)

        li_maxnptsx = api_k[0] + 1
        li_maxnptsy = api_k[1] + 1

        self.maxnpts = max (li_maxnptsx , li_maxnptsy)
        self.dirmaxnpts = max(li_maxnptsx , li_maxnptsy)

        self.maxnptsx = [li_maxnptsx,li_maxnptsy]
        self.n  = [li_nx,li_ny]

#        print "lpr_localx=", lpr_localx
#        print "lpr_localwx=", lpr_localwx
#        print "lpr_localy=", lpr_localy
#        print "lpr_localwy=", lpr_localwy

        self.npts = 2 * ( li_nx + li_ny )

        li_d = 0
        lr_xmin, lr_xmax = get_bounds(self.patch,li_d)
        li_d = 1
        lr_ymin, lr_ymax = get_bounds(self.patch,li_d)

        list_pts = []                           ; list_w = []

        list_pts.append([lr_ymin, lpr_localx])   ; list_w.append([1.0, lpr_localwx])
        list_pts.append([lr_xmin, lpr_localy])   ; list_w.append([1.0, lpr_localwy])
        list_pts.append([lr_ymax, lpr_localx])   ; list_w.append([1.0, lpr_localwx])
        list_pts.append([lr_xmax, lpr_localy])   ; list_w.append([1.0, lpr_localwy])

        self.pts = list_pts
        self.w  = list_w

    def gen_grid_3D(self, api_k, as_type):
        print "boundary_grid : gen_grid_3D not done yet"

    def save_grid_1D(self):
        print "boundary_grid : save_grid_1D not done yet"

    def save_grid_2D(self):
        li_curelt = 0
        for face in self.faces:
            axis = face % self.dim
            n,m = self.pts[face][1].shape

            li_curelt = self.com.pyfem.pyfem_set_unidirection_tensor_2d ( self.current_grids \
            , self.id, axis, li_curelt  \
            , self.pts[face][0], self.w[face][0]  \
            , self.pts[face][1], self.w[face][1]  \
            , n, m)

#    def save_grid_2D(self):
#        li_curelt = 0
#        # ....... x-direction
#        li_direction = 0
#        # 1st bound
#        li_bound = 0
#        li_curelt = self.com.pyfem.pyfem_set_unidirection_tensor_2d ( self.current_grids \
#        , self.id, li_direction, li_curelt  \
#        , self.pts[li_bound][0], self.w[li_bound][0]  \
#        , self.pts[li_bound][1], self.w[li_bound][1]  \
#        , self.n[li_direction], self.maxnptsx[li_direction]+1 )
#        # 2nd bound
#        li_bound = 1
#        li_curelt = self.com.pyfem.pyfem_set_unidirection_tensor_2d ( self.current_grids \
#        , self.id, li_direction, li_curelt  \
#        , self.pts[li_bound][0], self.w[li_bound][0]  \
#        , self.pts[li_bound][1], self.w[li_bound][1]  \
#        , self.n[li_direction], self.maxnptsx[li_direction]+1 )
#        # ....... y-direction
#        li_direction = 1
#        # 1st bound
#        li_bound = 2
#        li_curelt = self.com.pyfem.pyfem_set_unidirection_tensor_2d ( self.current_grids \
#        , self.id, li_direction, li_curelt  \
#        , self.pts[li_bound][0], self.w[li_bound][0]  \
#        , self.pts[li_bound][1], self.w[li_bound][1]  \
#        , self.n[li_direction], self.maxnptsx[li_direction]+1 )
#        # 2nd bound
#        li_bound = 3
#        li_curelt = self.com.pyfem.pyfem_set_unidirection_tensor_2d ( self.current_grids \
#        , self.id, li_direction, li_curelt  \
#        , self.pts[li_bound][0], self.w[li_bound][0]  \
#        , self.pts[li_bound][1], self.w[li_bound][1]  \
#        , self.n[li_direction], self.maxnptsx[li_direction]+1 )

    def save_grid_3D(self):
        print "boundary_grid : save_grid_3D not done yet"

    def get_real_elts_1D(self):
        # default
        print "boundary_grid : get_real_elts_1D : not done yet"

    def get_real_elts_2D(self):
        # default
        import numpy as np

        # x-direction
        li_d = 0
        li_n = self.patch.shape[li_d]
        li_p = self.patch.degree[li_d]
        li_un = len(np.unique(self.patch.knots[li_d]))
        ll_isperiodic = (li_un == li_n + li_p + 1)
        if ll_isperiodic :
            print "get_real_elts_2D : is incompatible with periodic obundary condition"
        else:
            list_u = np.unique(self.patch.knots[li_d])
            li_nel = len(list_u)-1
        list_ux = np.array(list_u)

        # y-direction
        li_d = 1
        li_n = self.patch.shape[li_d]
        li_p = self.patch.degree[li_d]
        li_un = len(np.unique(self.patch.knots[li_d]))
        ll_isperiodic = (li_un == li_n + li_p + 1)
        if ll_isperiodic :
            print "get_real_elts_2D : is incompatible with periodic obundary condition"
        else:
            list_u = np.unique(self.patch.knots[li_d])
            li_nel = len(list_u)-1
        list_uy = np.array(list_u)

        li_nex = len(list_ux)-1
        li_ney = len(list_uy)-1

        list_real_elts = []

        # bottom boundary
        if 0 in self.faces:
            li_j = 0
            for li_i in range(0,li_nex):
                li_index = li_j * li_nex + li_i
                list_real_elts.append(li_index+1)

        # upper boundary
        if 2 in self.faces:
            li_j = li_ney - 1
            for li_i in range(0,li_nex):
                li_index = li_j * li_nex + li_i
                list_real_elts.append(li_index+1)

        # left boundary
        if 1 in self.faces:
            li_i = 0
            for li_j in range(0,li_ney):
                li_index = li_j * li_nex + li_i
                list_real_elts.append(li_index+1)

        # right boundary
        if 3 in self.faces:
            li_i = li_nex - 1
            for li_j in range(0,li_ney):
                li_index = li_j * li_nex + li_i
                list_real_elts.append(li_index+1)

        return np.array(list_real_elts)


    def get_real_elts_3D(self):
        # default
        print "boundary_grid : get_real_elts_3D : not done yet"

    def get_real_elts(self):
        if self.dim==1:
            return self.get_real_elts_1D()
        if self.dim==2:
            return self.get_real_elts_2D()
        if self.dim==3:
            return self.get_real_elts_3D()
        if self.dim not in [1,2,3]:
            print "boundary_grid : get_real_elts : dimension not done yet"


if __name__ == '__main__':
    from igakit.cad_geometry import square as domain
    geo = domain(n=[3,2], p=[1,1])
    k = [2,2]
    nrb = geo[0]
    grids_id = 0
    id = 0
    grid = boundary_grid(grids_id, id, nrb, k, faces=[0,1])
