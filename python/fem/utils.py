# -*- coding: UTF-8 -*-
import numpy as np

asarray = np.asarray
zeros_like = np.zeros_like

__all__ = ['assembler']

def evaluator(dim, sites, EVALUATE, APPEND, nx=None, npts=None, elts=None, fmt=True):

    if elts is None:
        # ... evaluate on a single patch
        list_values = []
        list_x = []

        for x in sites:
            array_x = x.reshape(x.size)
            list_x.append(array_x)
        list_S = EVALUATE(*list_x)
        np.set_printoptions(precision=3)
#        print ">>> list_S"
#        print "shape ", list_S.shape
#        print list_S[:,:,0].transpose()[::-1,:]
        if APPEND:
            list_S = list_S[np.newaxis, ...]
        shp = list_S.shape

        n = 1
        for (_n,_npts) in zip(nx, npts):
            n *= _n * _npts
        Nv = shp[0]
        Nx = shp[-1]
        list_values = np.zeros((Nx,Nv,n))

        for iv in range(0, Nv):
            for ix in range(0,Nx):
                values = list_S[iv,...,ix]
                Z = []
                # ...
                if fmt:
                    if dim == 1:
                        for i in range(0,nx[0]):
                                ib = i * npts[0] ; ie = (i+1) * npts[0]
                                d = values[ib:ie]
                                Z += list(d.reshape(d.size))
                    if dim == 2:
                        for j in range(0,nx[1]):
                            for i in range(0,nx[0]):
                                ib = i * npts[0] ; ie = (i+1) * npts[0]
                                jb = j * npts[1] ; je = (j+1) * npts[1]
                                d = values[ib:ie,jb:je]
                                Z += list(d.transpose().reshape(d.size))
                    if dim == 3:
                        # TODO  to validate
                        for k in range(0,nx[2]):
                            for j in range(0,nx[1]):
                                for i in range(0,nx[0]):
                                    ib = i * npts[0] ; ie = (i+1) * npts[0]
                                    jb = j * npts[1] ; je = (j+1) * npts[1]
                                    kb = k * npts[2] ; ke = (k+1) * npts[2]
                                    d = values[ib:ie,jb:je,kb:ke]
                                    Z += list(d.transpose().reshape(d.size))

                    Z = asarray(Z)
                    list_values[ix, iv, :] = Z
        return list_values

    else:
        values = []
        # ... evaluate on a single element
        if dim == 1:
            x = sites
            for elt in elts:
                i = elt
                S = EVALUATE(x[i])
                if fmt:
                    c = S[...,0]
                    values += list(c.reshape(c.size))
                else:
                    values.append(S[...,0])
        if dim == 2:
            x,y = sites
            n = len(x)
            m = len(y)
            for elt in elts:
                j = elt / n
                i = elt - n * j
                S = EVALUATE(x[i], y[j])
                if fmt:
                    c = S[...,0]
                    values += list(c.transpose().reshape(c.size))
                else:
                    values.append(S[...,0])
        if dim == 3:
            # TODO  to validate
            x,y,z = sites
            n = len(x)
            m = len(y)
            p = len(z)
            for elt in elts:
                k = elt / m * n
                j = ( elt - m * n * k ) / n
                i = elt - ( m * n * k + n * j )
                S = EVALUATE(x[i], y[j], z[k])
                if fmt:
                    c = S[...,0]
                    values += list(c.transpose().reshape(c.size))
                else:
                    values.append(S[...,0])
        # ...

    return values


class assembler(object):
    def __init__(self, pdes=None, attributs={}):
        from . import fem      as fem
        self.fem = fem.fem()
        self.pdes = pdes

    def assembly(self, pdes=None, matrices=True, fields=True, norms=True):
        """
        this function assemblies a list of pdes
        """
        if pdes is None:
            pdes = self.pdes

        # prepare the assembly  process (rise boundary fcts)
        for pde in pdes:
            pde.prepare_assembly()

        list_matrices = []
        if matrices :
            for pde in pdes:
                list_matrices += pde.matrices

        list_fields = []
        if fields :
            for pde in pdes:
                list_fields += pde.fields

        list_norms = []
        if norms :
            for pde in pdes:
                list_norms += pde.norms

#        print "****************"
#        print "matrices : ", [M.id for M in list_matrices]
#        print "fields : ", [F.id for F in list_fields]
#        print "norms : ", [N.id for N in list_norms]

        self.fem.assembly(  matrices    = list_matrices \
                          , fields      = list_fields   \
                          , norms       = list_norms    )

#        print "****************"

class function(object):
    def __init__(self, func, fields=None, space=None):
        if space is not None:
            self.space = space
        else:
            if fields is not None:
                self.space = fields[0].space
            else:
                raise("Can not access to the discret vectorial space")
        self.fields = fields
        self.func = func
        self.dim = self.space.dim
#        print "========"
#        print "dim    ", self.dim
#        print "fields ", self.fields

    def __call__(self, xyz):
#        print ">>>"
#        print "__CALL__ with ", len(xyz), " arguments"
        arglist = []
        if self.fields is not None:
            for F in self.fields:
#                print "=== currentPatchID : ", F.space.currentPatchID
#                print "=== id-space ", id(F.space)
#                print "=== id-field ", id(F)

                arglist.append(F)

        for x in xyz:
            arglist.append(x)
        f = asarray(self.func(*arglist))
        _f = asarray([x[:] - x[:] + u for u in f])
#        print "<<<"
        return _f
