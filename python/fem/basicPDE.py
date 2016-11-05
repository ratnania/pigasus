# -*- coding: UTF-8 -*-
#! /usr/bin/python

import sys
import numpy as np
from scipy.sparse.linalg import cg as ConjugateGradient
from caid.cad_geometry import cad_geometry
from caid.cad_geometry import entier_vers_couple as face_to_bc

from .constants import *
from .field import *
from .norm import *
from .grids import *
from .matrix import *
from .space import *
from .graph import *
from .oper import *
from .pigasusClass import pigasus
from .solver import *

from numpy import zeros

__all__ = ['basicPDE', 'multi_basicPDE']

class boundary_function():
    def __init__(self, geo, patch_id, face_id, func_g):
        self.geo        = geo
        self.patch      = geo[patch_id]
        self.patch_id   = patch_id
        self.face_id    = face_id
        self.func_g     = func_g

    def __call__(self, u):
#        print "* evaluate *"
#        print "patch_id, face : ", self.patch_id, self.face_id
        nrb = self.patch
        axis, side = face_to_bc(self.face_id)
        nrb_bnd = nrb.extract_face(axis, side).clone()
        sgn = nrb_bnd.orientation[0]
        P = nrb_bnd.evaluate_deriv(u)
        x  = P[0,:,0]
        y  = P[0,:,1]
        dx = P[1,:,0]
        dy = P[1,:,1]
        list_g = self.func_g(x,y)
        r = np.sqrt(dx**2 + dy**2)
        list_n = [ sgn * dy/r, -sgn *dx/r]
        r = list_g[0]*list_n[0]
        for (g,n) in zip(list_g[1:], list_n[1:]):
            r += g * n
        return [r]

class basicPDE(pigasus):
    """
    A multidimentional basicPDE class solver.
    This class supports the following 5 differential operators: MASS, ADVECTION,
    transpose(ADVECTION), STIFFNESS and
    SECOND_DERIV. It provides the assembling function, which assemblies the
    linear system, and the solve, plot, norm and get functions, as explained
    next.
    We start by importing the following modules
        >>> import caid.cad_geometry  as cg
        >>> from caid.cad_geometry import line as domain
        >>> import pylab                as pl
        >>> import numpy                as np

    creation of the geometry
        >>> geo = domain(n=nx, p=px)

    creation of the testcase: it must be a dictionary
        >>> testcase = {}

    we need to specify the list of faces for which we apply Homogeneous Dirichlet boundary condition
        >>> testcase['list_DirFaces'] = [[1,2]]

    define of differentials operators by specifying their "functions"
        >>> # testcase['mass']       = lambda x : [0.]
        >>> # testcase['adv']        = lambda x : [0.]
        >>> testcase['stiffness']  = lambda x : [1.]

    define of the exact solution (if known) and the right hand side  (source/load term).
        >>> kx = 2. * pi
        >>> # exact solution
        >>> u = lambda x : sin ( kx * x )
        >>> testcase['u'] = u
        >>> # rhs
        >>> f = lambda x : ( kx**2) * sin ( kx * x )
        >>> testcase['f'] = f

    specify the values of the solution on the boundary.
        >>> # values of u at the boundary
        >>> g1 = lambda x: u(x)
        >>> g2 = lambda x: u(x)
        >>> testcase['list_faces'] = [[1,2]]
        >>> testcase['list_g'] = [[g1, g2]]

    Finally, the script for solving the PDE becomes:
        >>> PDE = basicPDE(geometry=geo, testcase=testcase)
        >>> PDE.assembly()
        >>> PDE.solve()
        >>> PDE.plot(); pl.show()
        >>> PDE.norm()
        >>> U = PDE.get() # get the unknown field
        >>> u = U.get()   # get the B-splines/NURBS coeffs as a vector
        >>> U.tomatrix()  # get the B-splines/NURBS coeffs as a list of matrices (for each patch)

    """

    #: Doc comment for class attribute gallery.basicPDE.
    #: It can have multiple lines.
    def __init__(self, *args, **kwargs):
        """Creates an Elliptic PDE solver. arguments are the same as pigasus.__init__

        V :
            if specified, this tells pigasus to use the space V as self.V

        U_V :
            if specified, this tells pigasus to use this given field as (the unknown).

        F_V :
            if specified, this tells pigasus to use this given field as (the rhs).

        M_V :
            if specified, this tells pigasus to use this given matrix as mass.

        A_V :
            if specified, this tells pigasus to use this given matrix as advection.

        At_V :
            if specified, this tells pigasus to use this given matrix as transpose Advection.

        S_V :
            if specified, this tells pigasus to use this given matrix as stiffness.

        B_V :
            if specified, this tells pigasus to use this given matrix as second_deriv.

        Returns:
           A PDE object.

        .. note::

            See also: :ref:`fem.pigasus`.

        """
#        verbose = True
        verbose = False

        pigasus.__init__(self, *args, **kwargs)

        # ...
        # non homogeneous Dirichlet Boundary conditions
        # ...
        try:
            self.bc_dirichlet = self.testcase['bc_dirichlet']
            self.Dirichlet = True
        except:
            self.bc_dirichlet = {}
            self.Dirichlet = False
        # ...

        # ...
        # neumann Boundary conditions
        # ...
        try:
            self.bc_neumann = self.testcase['bc_neumann']
            self.Neumann = True
        except:
            self.bc_neumann = {}
            self.Neumann = False
        # ...

        # ...
        try:
            func_mass = self.testcase['b']
            self.withMass = True
        except:
            self.withMass = False
        # ...

        # ...
        try:
            func_adv = self.testcase['v']
            self.withAdvection = True
        except:
            self.withAdvection = False
        # ...

        # ...
        try:
            func_tadv = self.testcase['w']
            try:
                func_dw   = self.testcase['dw']
            except:
                print("Warning: You should give the (div w) term. Pigasus is not yet smart enough!")
                # ...
                if self.dim == 1:
                    def func_dw(x):
                        eps = 1.e-3
                        return [(w(x+eps)[0]-w(x)[0])/eps]
                if self.dim == 2:
                    def func_dw(x,y):
                        eps = 1.e-3
                        return [  (w(x+eps,y)[0]-w(x,y)[0])/eps \
                                + (w(x,y+eps)[1]-w(x,y)[1])/eps]
                if self.dim == 3:
                    def func_dw(x,y,z):
                        eps = 1.e-3
                        return [  (w(x+eps,y,z)[0]-w(x,y,z)[0])/eps \
                                + (w(x,y+eps,z)[1]-w(x,y,z)[1])/eps \
                                + (w(x,y,z+eps)[2]-w(x,y,z)[2])/eps]
                # ...
            self.withTAdvection = True
            self.withMass       = True
        except:
            self.withTAdvection = False
        # ...

        # ...
        try:
            func_stiff = self.testcase['A']
            self.withStiffness = True
        except:
            self.withStiffness = False
        # ...

        # ...
        try:
            func_D2 = self.testcase['D2']
            self.withD2 = True
        except:
            self.withD2 = False
        # ...

        # ...
        self.withMetric = False
        try:
            Metric = self.testcase['metric']
            if Metric is not None:
                self.withMetric = True
        except:
            pass
        # ...

        # ...
        withAllDirichlet = False
        try:
            if self.testcase['AllDirichlet'] is not None:
                withAllDirichlet = self.testcase['AllDirichlet']
        except:
            pass
        # ...

        # ...
        try:
            solverInfo = self.testcase['solverInfo']
            if solverInfo is not None:
                self.solverInfo = {}
        except:
            self.solverInfo = None
        # ...

        # ... set geometry
        if self.geometry is None:
            try:
                V  = kwargs['V']
                self.geometry = V.geometry
            except:
                raise("Unable to find a geometry for the current basicPDE")
        # ...

        self.withAllDirichlet = withAllDirichlet

        # ...
        list_DirFaces = []
        for i in range(0, self.geometry.npatchs):
            list_DirFaces.append([])
        if withAllDirichlet:
            list_extFaces = self.geometry.external_faces
            for extFaces in list_extFaces:
                patch_id    = extFaces[0]
                face_id     = extFaces[1]
                list_DirFaces[patch_id].append(face_id)
        else:
            try:
                list_DirFaces = self.testcase['Dirichlet']
            except:
                pass
        # ...

        self.list_DirFaces = list_DirFaces

        # ...
        self.UseDuplicateFaces          = False
        list_DuplicatedFaces            = []
        list_DuplicataFaces             = []
        list_DuplicatedFacesPeriodic    = []
        try:
            list_connectivity = self.testcase['connectivity']
        except:
            list_connectivity = self.geometry.connectivity

        for dict_con in list_connectivity:
            list_DuplicatedFaces.append(dict_con['original'])
            list_DuplicataFaces.append(dict_con['clone'])
            try:
                list_DuplicatedFacesPeriodic.append(dict_con['periodic'])
            except:
                list_DuplicatedFacesPeriodic.append(False)

        if len(list_DuplicatedFaces) > 0:
            self.UseDuplicateFaces = True
        # ...

        self.list_DuplicataFaces            = list_DuplicataFaces
        self.list_DuplicatedFaces           = list_DuplicatedFaces
        self.list_DuplicatedFacesPeriodic    = list_DuplicatedFacesPeriodic

        # ...
        self.meanConstraint = False
        nExtFaces           = len(self.geometry.external_faces)
        nBCNeumann          = len(self.bc_neumann)
        nDuplicatedFaces    = len(self.list_DuplicataFaces)
        if (nExtFaces == nBCNeumann) \
           or (nExtFaces == nDuplicatedFaces) \
           or (nExtFaces == nBCNeumann+nDuplicatedFaces) \
           or (nExtFaces == 0):
            self.meanConstraint = True
#            print (" self.meanConstraint ", self.meanConstraint)
        # ...

        #-----------------------------------
        nrb = self.geometry[0]
        list_n = nrb.shape
        list_p = nrb.degree

        # ...
        try:
            n_gauss  = kwargs['n_gauss']
            lpi_ordergl = n_gauss
        except:
            lpi_ordergl = list_p

        _system = matrix()
        _slv    = solver(matrix=_system, solverInfo=self.solverInfo)
        #-----------------------------------

        #-----------------------------------
        if self.dim == 1:
            func_zero  = lambda x        : [ 0. ]
            func_one   = lambda x        : [ 1. ]
        if self.dim == 2:
            func_zero  = lambda x,y      : [ 0. ]
            func_one   = lambda x,y      : [ 1. ]
        if self.dim == 3:
            func_zero  = lambda x,y,z    : [ 0. ]
            func_one   = lambda x,y,z    : [ 1. ]
        #-----------------------------------

        # ...
        try:
            V  = kwargs['V']
        except:
            #-----------------------------------
            # space declaration
            #-----------------------------------
            list_DirFaces_V = list_DirFaces
            if self.Dirichlet:
                list_DirFaces_V = []
                for i in range(0,self.geometry.npatchs):
                    list_DirFaces_V.append([])
                list_extFaces = self.geometry.external_faces
                for extFaces in list_extFaces:
                    patch_id    = extFaces[0]
                    face_id     = extFaces[1]
                    list_DirFaces_V[patch_id].append(face_id)
            self.list_DirFaces_V = list_DirFaces_V

            if verbose:
                print("* DirFaces_V      ", self.list_DirFaces_V)
                print("* DuplicatedFaces ", list_DuplicatedFaces)
                print("* DuplicataFaces  ", list_DuplicataFaces)

#            print("* DuplicatedFaces ", list_DuplicatedFaces)
#            print("* DuplicataFaces  ", list_DuplicataFaces)
#            print("* DuplicatedFacesPeriodic ", list_DuplicatedFacesPeriodic)

            V = space(geometry=self.geometry)
            V.dirichlet(faces=list_DirFaces_V)
            if self.UseDuplicateFaces:
                V.duplicate(  faces_base = list_DuplicatedFaces \
                            , faces      = list_DuplicataFaces  \
                            , isPeriodic = list_DuplicatedFacesPeriodic)
            V.set_boundary_conditions()
            if self.withMetric:
                V.create_grids(type="legendre", k=lpi_ordergl, metric=Metric)
            else:
                V.create_grids(type="legendre", k=lpi_ordergl)
            #-----------------------------------
        # ...

        # ...
        try:
            GrV  = kwargs['GrV'] # [V,V] graph
        except:
            #-----------------------------------
            # graph declaration
            #-----------------------------------
            GrV = graph(spaces=[V,V])
            #-----------------------------------
        # ...

        # ...
        #-----------------------------------
        # Matrix declaration
        #-----------------------------------
        Matrix_V  = matrix(graph=GrV)
        MatProj_V = matrix(graph=GrV)
        #-----------------------------------
        # ...

        # ...
        if self.Dirichlet:
            try:
                W  = kwargs['W']
            except:
                #-----------------------------------
                # space declaration
                #-----------------------------------
                W = space(geometry=self.geometry)
                W.dirichlet(faces=[[]] * self.geometry.npatchs)
                if self.UseDuplicateFaces:
                    W.duplicate(  faces_base = list_DuplicatedFaces \
                                , faces      = list_DuplicataFaces  )
                W.set_boundary_conditions()
#                W.create_grids(type="legendre", k=lpi_ordergl)
                W.grids = V.grids
                #-----------------------------------
        # ...

        # ...
        if self.Dirichlet:
            try:
                GrVW  = kwargs['GrVW']
            except:
                #-----------------------------------
                # graph declaration
                #-----------------------------------
                GrVW = graph(spaces=[V,W])
                #-----------------------------------
        # ...

        # ...
        if self.Dirichlet:
            #-----------------------------------
            # Matrix declaration
            #-----------------------------------
            Matrix_VW = matrix(graph=GrVW)
            #-----------------------------------
        # ...

        # ...
        try:
            func_u = self.testcase['u']
        except:
            func_u = func_zero

        try:
            U_V  = kwargs['U_V']
        except:
            U_V = field(space=V, func = func_u)

        if self.Dirichlet:
            try:
                U_W  = kwargs['U_W']
            except:
                U_W = field(space=W, func = func_u)
        # ...

        # ...
        try:
            func_f = self.testcase['f']
        except:
            func_f = func_zero

        try:
            F_V  = kwargs['F_V']
        except:
            F_V = field(space=V, func = func_f)

        if self.Dirichlet:
            func_w = func_one
            try:
                G_W  = kwargs['G_W']
            except:
                G_W = field(space=W, func = func_w)
        # ...

        # ... Temp field
        T_V = field(space=V, func = func_zero)
        # ...

        # ...
        if self.Dirichlet:
            try:
                U_W  = kwargs['U_W']
            except:
                U_W = field(space=W, func = func_u)
        # ...

        # ...
        try:
            Projector_V  = kwargs['Projector_V']
        except:
            if V.dim == 1:
                func_one  = lambda x        : [ 1. ]
            if V.dim == 2:
                func_one  = lambda x,y      : [ 1. ]
            if V.dim == 3:
                func_one  = lambda x,y,z    : [ 1. ]
            Projector_V  = oper(spaces=[V, V], type=MASS        , func=func_one)
            Projector_V.addto(MatProj_V)
        # ...

        # ...
        try:
            trial = kwargs['trial']
        except:
            trial = V
        # ...

        # ...
        try:
            M_V  = kwargs['M_V']
        except:
            if self.withMass:
                M_V  = oper(spaces=[V, trial], type=MASS        , func=func_mass)
                M_V.addto(Matrix_V)
                # ...
                if self.withTAdvection:
                    Mdw_V  = oper(spaces=[V, trial], type=MASS  , func=func_dw)
                    Mdw_V.addto(Matrix_V)
                # ...

        # TODO what to do with the trial space for VW?
        try:
            M_VW  = kwargs['M_VW']
        except:
            if self.withMass and self.Dirichlet:
                M_VW  = oper(spaces=[V, W], type=MASS      , func=func_mass)
                M_VW.addto(Matrix_VW)
        # ...

        # ...
        try:
            A_V  = kwargs['A_V']
        except:
            if self.withAdvection:
                A_V  = oper(spaces=[V, trial], type=ADVECTION   , func=func_adv)
                A_V.addto(Matrix_V)

        try:
            A_VW  = kwargs['A_VW']
        except:
            if self.withAdvection and self.Dirichlet:
                A_VW  = oper(spaces=[V, W], type=ADVECTION   , func=func_adv)
                A_VW.addto(Matrix_VW)
        # ...

        # ...
        try:
            At_V  = kwargs['At_V']
        except:
            if self.withTAdvection:
                At_V  = oper(spaces=[V, trial], type=ADVECTION  , func=func_tadv, transpose=True)
                At_V.addto(Matrix_V)

        try:
            At_VW  = kwargs['At_VW']
        except:
            if self.withTAdvection and self.Dirichlet:
                At_VW  = oper(spaces=[V, W], type=ADVECTION  , func=func_tadv, transpose=True)
                At_VW.addto(Matrix_VW)
        # ...

        # ...
        try:
            S_V  = kwargs['S_V']
        except:
            if self.withStiffness:
                S_V  = oper(spaces=[V, trial], type=STIFFNESS   , func=func_stiff)
                S_V.addto(Matrix_V)

        try:
            S_VW  = kwargs['S_VW']
        except:
            if self.withStiffness and self.Dirichlet:
                S_VW  = oper(spaces=[V, W], type=STIFFNESS   , func=func_stiff)
                S_VW.addto(Matrix_VW)
        # ...

        # ...
        try:
            B_V  = kwargs['B_V']
        except:
            if self.withD2:
                B_V  = oper(spaces=[V, trial], type=SECOND_DERIV, func=func_D2)
                B_V.addto(Matrix_V)

        try:
            B_VW  = kwargs['B_VW']
        except:
            if self.withD2 and self.Dirichlet:
                B_VW  = oper(spaces=[V, W], type=SECOND_DERIV, func=func_D2)
                B_VW.addto(Matrix_VW)
        # ...

        # ...
        self.list_G_V_BC_faces = []
        try:
            V_BC  = kwargs['V_BC']
        except:
            if self.Neumann:
                list_V_BC = []
                list_G_V_BC = []
                #-----------------------------------
                # boundary space declaration
                #-----------------------------------
                for key, func_g in self.bc_neumann.items():
                    patch_id    = int(key[0])
                    face_id     = int(key[1])

                    self.list_G_V_BC_faces.append([patch_id,face_id])

                    nrb = self.geometry[patch_id]

                    axis, side = face_to_bc(face_id)

                    nrb_bnd = nrb.extract_face(axis, side)
                    geo = cad_geometry()
                    geo.append(nrb_bnd)

                    lpi_ordergl = nrb_bnd.degree

                    # ...
                    V_BC = space(geometry=geo)
                    V_BC.dirichlet(faces=[[]])
                    V_BC.set_boundary_conditions()
                    V_BC.create_grids(type="legendre", k=lpi_ordergl)
                    # ...

                    # ...
                    g = boundary_function(self.geometry\
                                          , patch_id\
                                          , face_id\
                                          , func_g)
                    G_V_BC  = field(space=V_BC, pfunc=g)
                    # ...

                    # ...
                    list_V_BC.append(V_BC)
                    list_G_V_BC.append(G_V_BC)
                    # ...
                #-----------------------------------
        # ...

        # ...
        try:
            Mean_V  = kwargs['Mean_V']
        except:
            if self.meanConstraint:
                if self.dim == 1:
                    func_one  = lambda x        : [ 1. ]
                if self.dim == 2:
                    func_one  = lambda x,y      : [ 1. ]
                if self.dim == 3:
                    func_one  = lambda x,y,z    : [ 1. ]
                Mean_V = field(space=V, func = func_one)
        # ...

        # ...
        try:
            N_U  = kwargs['N_U']
        except:
            #-----------------------------------
            # NORM
            #-----------------------------------
            if self.Dirichlet:
                N_U = norm(field=U_W, type=NORM_L2, exact=func_u)
            else:
                N_U = norm(field=U_V, type=NORM_L2, exact=func_u)
            #-----------------------------------
        # ...

        #-----------------------------------
        # Save access for data
        #-----------------------------------
        self._slv           = _slv
        self._system        = _system
        self.Matrix_V       = Matrix_V
        self.MatProj_V      = MatProj_V
        self.GrV            = GrV

        self.spaces         = []
        self.norms          = []
        self.fields         = []
        self.graphs         = [self.GrV]
        self.matrices       = [self.Matrix_V]

        self.V              = V
        self.spaces        += [self.V]

        if self.withMetric:
            self.Metric     = Metric

        if self.Dirichlet:
            self.W          = W
            self.GrVW       = GrVW
            self.spaces    += [self.W]
            self.graphs    += [self.GrVW]

        self.N_U            = N_U
        self.norms         += [self.N_U]

        # ... rhs
        self.F_V            = F_V
        self.fields        += [self.F_V]
        # ... L2 projector
        self.Projector_V    = Projector_V

        if self.Dirichlet:
            self.G_W        = G_W
            # insteed fo assembling this field, we will rise it from the
            # boundary
#            self.fields    += [self.G_W]
        # ...

        # ... unknown
        self.U_V            = U_V
        self.fields        += [self.U_V]
        if self.Dirichlet:
            self.U_W        = U_W
#            self.fields    += [self.U_W]
        # ...

        # ... Temp field
        self.T_V            = T_V
        # ...

        # ... Mass operator
        if self.withMass:
            self.M_V        = M_V
            if self.withTAdvection:
                self.Mdw_V  = Mdw_V
            if self.Dirichlet:
                self.M_VW   = M_VW
        # ...

        # ... Advection operator
        if self.withAdvection:
            self.A_V        = A_V
            if self.Dirichlet:
                self.A_VW   = A_VW
        # ...

        # ... Transpose Advection operator
        if self.withTAdvection:
            self.At_V       = At_V
            if self.Dirichlet:
                self.At_VW  = At_VW
        # ...

        # ... Stiffness operator
        if self.withStiffness:
            self.S_V        = S_V
            if self.Dirichlet:
                self.S_VW   = S_VW
        # ...

        # ... Second Derivative operator
        if self.withD2:
            self.B_V        = B_V
            if self.Dirichlet:
                self.B_VW   = B_VW
        # ...

        # ...
        if self.Dirichlet:
            self.Matrix_VW  = Matrix_VW
            self.matrices  += [self.Matrix_VW]
        # ...

        # ...
        if self.Neumann:
            self.list_V_BC  = list_V_BC
            for S in list_V_BC:
                self.spaces+= [S]

            self.list_G_V_BC= list_G_V_BC
            for G in list_G_V_BC:
                self.fields+= [G]
        # ...

        # ... Mean Constraint in the case of Pure Neumann and Periodic bc
        if self.meanConstraint:
            self.Mean_V     = Mean_V
            self.fields    += [self.Mean_V]
        # ...

        self.forceAssembly  = False
        self.Assembled      = False

        #-----------------------------------

    #-----------------------------------

    #-----------------------------------
    def __del__(self):
        self.fem.__del__()
    #-----------------------------------

    #-----------------------------------
    def free(self):
        self.fem.__del__()
    #-----------------------------------

    #-----------------------------------
    @property
    def dim(self):
        """
        returns the parametric dimension
        """
        return self.geometry.dim
    #-----------------------------------

    #-----------------------------------
    @property
    def size(self):
        """
        returns the vectorial space dimension
        """
        return self.V.size
    #-----------------------------------

    #-----------------------------------
    @property
    def shape(self):
        """
        returns the shape of the linear system
        """
        return self._system.shape
    #-----------------------------------

    #-----------------------------------
    @property
    def space(self):
        """
        returns the vectorial space
        """
        return self.V
    #-----------------------------------

    #-----------------------------------
    @property
    def mass(self):
        """
        returns the MASS operator if defined, None otherwise.
        """
        if self.withMass:
            return self.M_V
        else:
            return None
    #-----------------------------------

    #-----------------------------------
    @property
    def stiffness(self):
        """
        returns the STIFFNESS operator if defined, None otherwise.
        """
        if self.withStiffness:
            return self.S_V
        else:
            return None
    #-----------------------------------

    #-----------------------------------
    @property
    def advection(self):
        """
        returns the ADVECTION operator if defined, None otherwise.
        """
        if self.withAdvection:
            return self.A_V
        else:
            return None
    #-----------------------------------

    #-----------------------------------
    @property
    def tadvection(self):
        """
        returns the transpose of ADVECTION operator if defined, None otherwise.
        """
        if self.withTAdvection:
            return self.At_V
        else:
            return None
    #-----------------------------------

    #-----------------------------------
    @property
    def D2(self):
        """
        returns the SECOND_DERIV operator if defined, None otherwise.
        """
        if self.withD2:
            return self.B_V
        else:
            return None
    #-----------------------------------

    #-----------------------------------
    @property
    def system(self):
        """
        returns the linear system as a pigasus matrix object
        """
        return self._system
    #-----------------------------------

    #-----------------------------------
    def plot(self, vmin=None, vmax=None, fastplot=True):
        """
        Plots the solution of the PDE, using the fast_plot function
        """
        geo = self.geometry
        if self.Dirichlet:
            U = self.U_W
        else:
            U = self.U_V

        withpcolor = False
        if self.geometry.npatchs > 1:
            withpcolor = True
            if vmin is None:
                vmin = np.min(U.get())
            if vmax is None:
                vmax = np.max(U.get())
        if fastplot:
            for patch_id in range(0, geo.npatchs):
                U.fast_plot(ai_patch_id=patch_id \
                                   , withpcolor=withpcolor \
                                  , vmin=vmin \
                                  , vmax=vmax)
        else:
            U.plot(withpcolor)
    #-----------------------------------

    #-----------------------------------
    def dot(self, X):
        """
        returns the linear system as a pigasus matrix object

        X:
            may be a field or a numpy array
        """
        return self._system.dot(X)
    #-----------------------------------

    #-----------------------------------
    @property
    def unknown(self):
        """
        Returns the solution of the PDE.
        """
        return self.U_V
    #-----------------------------------

    #-----------------------------------
    @property
    def unknown_dirichlet(self):
        """
        Returns the solution of the PDE.
        """
        if self.Dirichlet:
            return self.U_W
        else:
            return None
    #-----------------------------------


    #-----------------------------------
    def get(self):
        """
        Returns the solution of the PDE.
        """
        if self.Dirichlet:
            return self.U_W
        else:
            return self.U_V
    #-----------------------------------

    #-----------------------------------
    def set(self, X):
        """
        sets the solution to the field or numpy array X.
        """
        self.U_V.set(X)
        if self.Dirichlet:
            self._fromVtoW ( self.U_V, self.U_W )
    #-----------------------------------

    #-----------------------------------
    @property
    def source(self):
        """
        returns the source term ie int f phi_i
        """
        return self.F_V
    #-----------------------------------

    #-----------------------------------
    @property
    def rhs(self):
        """
        returns the source term (rhs) including the boundary terms
        """
        return self.F_V
    #-----------------------------------

    #-----------------------------------
    @property
    def solver(self):
        """
        returns the solver instance
        """
        return self._slv
    #-----------------------------------

    #-----------------------------------
    @property
    def neumann_terms(self):
        """
        returns the Neumann integrals parts
        """
        if self.Neumann:
            return self.T_V
        else:
            return None
    #-----------------------------------

    #-----------------------------------
    def update(self, func=None):
        """
        assembly the right hand side
        """
        if func is not None:
            self.F_V.set_func(func)

        self.fem.assembly(fields=[self.F_V])
        self.post_assembly()
    #-----------------------------------

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
    def interpolate(self, u, field=None):
        if field is None:
            U = self.U_V
        else:
            U = field
        from pigasus.interpolate.interpolation import surfint
        interpolator = surfint(self.geometry)
        c_ini   = interpolator.interpolate(u)
        U.frommatrix(0, c_ini)

    #-----------------------------------
    def project(self, u, field=None):
        """
        computes the L2 projection of the function u
        """
        if field is None:
            U = self.U_V
        else:
            U = field

        u_old = self.U_V.func
        U.set_func(u)
        self.fem.initialize()
        self.prepare_assembly()
        #-----------------------------------
        # assembling everything
        #-----------------------------------
        U.reset()
        self.fem.assembly(matrices=[self.Matrix_V], fields=[U])
        #-----------------------------------
        _A = self.Matrix_V.get()
        _b = U.get()
        from scipy.sparse.linalg import cg as conj
        u_v = conj(_A,_b)[0]
        U.set(u_v)

        if self.Dirichlet:
            self._fromVtoW ( U, self.U_W )

        U.set_func(u_old)
    #-----------------------------------

    #-----------------------------------
    def assembly(self):
        """
        Assemblies the linear system to be solved.
        """
#        print "Enter poisson.assembly"
        self.fem.initialize()
#        print "%"
        self.prepare_assembly()
#        print "%%"
        #-----------------------------------
        # assembling everything
        #-----------------------------------
        self.fem.assembly(matrices=self.matrices, fields=self.fields)
        #-----------------------------------
#        print "%%%"
        M = self.Matrix_V.get()
#        print "%%%%"
#        _M = BlockMatrix([[M]])
#        _M.assembly()
#        self._system.set(_M.get())
        self._system.set(M)
#        print "%%%%%"
        # ...
        self.post_assembly()
#        print "%%%%%%"
        # ...
        self.Assembled = True
#        print "Leave poisson.assembly"
    #-----------------------------------

    #-----------------------------------
    def post_assembly(self):
        """
        """
        self.update_dirichlet_terms()
        self.update_neumann_terms()
    #-----------------------------------

    #-----------------------------------
    def prepare_assembly(self):
        """
        """
        if self.Dirichlet:
            # ...
            # Define a surface that interpolates the boundary conditions
            # ...
            G_W = self.G_W
            G_W.set_initialize()
            for i in range(0, self.geometry.npatchs):
                nrb = self.geometry[i]
                list_faces  = []
                list_g      = []
                for key, func in self.bc_dirichlet.items():
                    patch_id = key[0] ; face = key[1]
                    if i == patch_id:
                        list_faces.append(face) ; list_g.append(func)

                mat_g = G_W.riseFromBoundary(i, list_faces=list_faces, list_g=list_g)
                G_W.frommatrix(i, mat_g)
            G_W.set_finalize()
            # ...

#            print "-------"
#            print "prepare_assembly"
#            import matplotlib.pyplot as plt
#            G_W.plot(withpcolor=True) ; plt.colorbar() ; plt.show()
#            print "-------"

#        if self.withMetric:
#            V = self.V
#            lpr_parampts = V.get_sites()
#            print "============"
#            print lpr_parampts
#            self.Metric.update(apr_points=lpr_parampts)
    #-----------------------------------

    #-----------------------------------
    def load(self, filename):
        """
        load Matrix_V from file
        """
        self.Matrix_V.load(filename)
    #-----------------------------------

    #-----------------------------------
    def save(self, filename):
        """
        save Matrix_V from file
        """
        self.Matrix_V.save(filename)
    #-----------------------------------

    #-----------------------------------
    def update_dirichlet_terms(self, rhs=None):
        """
        updates the term F_V
        """
        if self.Dirichlet:
            if rhs is None:
                lpr_rhs = self.F_V.get()
            else:
                if isNumpyArray(rhs):
                    lpr_rhs = rhs
                if isField(rhs):
                    lpr_rhs = rhs.get()

            lpr_rhs -= self.Matrix_VW.dot(self.G_W.get())
            if rhs is None:
                self.F_V.set(lpr_rhs)
            else:
                if isNumpyArray(rhs):
                    rhs = lpr_rhs
                if isField(rhs):
                    rhs.set(lpr_rhs)
    #-----------------------------------

    #-----------------------------------
    def update_neumann_terms(self, rhs=None):
        """
        updates the term F_V and save the Neumann contribution into T_V
        """
        from .common_obj import isNumpyArray, isField
        if self.Neumann:
            if rhs is None:
                lpr_rhs = self.F_V.get()
            else:
                if isNumpyArray(rhs):
                    lpr_rhs = rhs
                if isField(rhs):
                    lpr_rhs = rhs.get()

            # ...
#            T_V = self.T_V
#            T_V.set_initialize()
#            for (G, pf) in zip(self.list_G_V_BC, self.list_G_V_BC_faces):
#                patch_id = pf[0]
#                face     = pf[1]
#                bnd_patch_id = 0
#
#                mat_rhs = self.F_V.tomatrix(patch_id)
#
#                vec_g = G.tomatrix(bnd_patch_id)
#                mat_g = np.zeros_like(mat_rhs)
#                if self.dim == 2:
#                    if face == 0:
#                        mat_g[:,0]  = vec_g[:]
#                    if face == 1:
#                        mat_g[0,:]  = vec_g[:]
#                    if face == 2:
#                        mat_g[:,-1] = vec_g[:]
#                    if face == 3:
#                        mat_g[-1,:] = vec_g[:]
#
#                T_V.frommatrix(patch_id, mat_g)
#            T_V.set_finalize()
#            lpr_g = T_V.get()
            # ...

            # ...
            lpr_g = np.zeros_like(lpr_rhs)
            for (G, pf) in zip(self.list_G_V_BC, self.list_G_V_BC_faces):
                patch_id = pf[0]
                face     = pf[1]
                bnd_patch_id = 0

                mat_rhs = self.F_V.tomatrix(patch_id)

                vec_g = G.tomatrix(bnd_patch_id)
                mat_g = np.zeros_like(mat_rhs)
                #Â TODO: must take the weighted sum of the corners
                # for the B-spline basis function that appears in two adjacent
                # faces
#                print "patch ", patch_id, " face ", face

                if self.dim == 2:
                    if face == 0:
                        mat_g[:,0]  = vec_g[:]
                    if face == 1:
                        mat_g[0,:]  = vec_g[:]
                    if face == 2:
                        mat_g[:,-1] = vec_g[:]
                    if face == 3:
                        mat_g[-1,:] = vec_g[:]

                self.T_V.reset()
                self.T_V.frommatrix(patch_id, mat_g)
                vec_g = self.T_V.get()
#                np.savetxt("vec_g_p"+str(patch_id)+"_f"+str(face)+".txt", vec_g)
                lpr_g += vec_g
#            np.savetxt("g.txt", lpr_g)
            self.T_V.set(lpr_g)
            # ...

            lpr_rhs += lpr_g
            if rhs is None:
                self.F_V.set(lpr_rhs)
            else:
                if isNumpyArray(rhs):
                    rhs = lpr_rhs
                if isField(rhs):
                    rhs.set(lpr_rhs)
    #-----------------------------------

    #-----------------------------------
    def solve(self, rhs=None):
        """
        solves the system assemblied after calling the assembly function.
        The solution is therefor stored in self.U_V (Homogeneous Dirichlet bc)
        or self.U_W (Dirichlet bc). The user can get it using the function get

        rhs:
           Can be either a field object or a numpy array. If not given, self.F_V
           will be used as rhs.
        """
#        from pigasus.utils.utils import process_diagnostics
#        print ">> solve"
#        process_diagnostics(30)
        # ...
        if rhs is None:
            lpr_rhs = self.F_V.get()
#            print "%"
        else:
            from .common_obj import isNumpyArray, isField
            if isNumpyArray(rhs):
                lpr_rhs = rhs
            if isField(rhs):
                lpr_rhs = rhs.get()

#        self._rhs    = BlockVector([[lpr_rhs]])
#        self._rhs.assembly()
#        _b = self._rhs.get().transpose()
        _b = lpr_rhs
#        print "%%"

        # ...
#        print "meanConstraint : ",  self.meanConstraint

        if self.meanConstraint:
            u_v = self._slv.solve(_b, guess=np.zeros_like(_b) \
                                 , c_constraint=self.Mean_V.get())
#            print "$$$"
        else:
            u_v = self._slv.solve(_b, guess=np.zeros_like(_b))
#            print "###"

#            print "*****"
#            print _b
#            print "*****"
#            print u_v
#            print "*****"

        self.U_V.set(u_v)
        # ...

        if self.Dirichlet:
            self._fromVtoW ( self.U_V, self.U_W )
    #-----------------------------------


    #-----------------------------------

    #-----------------------------------


    #-----------------------------------
    def _fromVtoW(self, X_V, X_W):
        """
        converts X_V to X_W
        """
        X_W.set_initialize()
        for i in range(0, self.geometry.npatchs):
            mat_v  = X_V.tomatrix(i)
            mat_w  = self.G_W.tomatrix(i)
            mat_w += mat_v
            X_W.frommatrix(i, mat_w)
        X_W.set_finalize()
    #-----------------------------------

    #-----------------------------------
    def _fromWtoV(self, X_W, X_V):
        """
        converts X_W to X_V
        """
        X_V.set_initialize()
        for i in range(0, self.geometry.npatchs):
            mat_w = X_W.tomatrix(i)
            mat_g = self.G_W.tomatrix(i)
            mat_v = mat_w - mat_g
            X_V.frommatrix(i, mat_v)
        X_V.set_finalize()
    #-----------------------------------


class multi_basicPDE(object):
    """
    A multidimentional basicPDE class.
        >>> import caid.cad_geometry  as cg

    """

    #: Doc comment for class attribute gallery.poisson.
    #: It can have multiple lines.
    def __init__(self, geometry, list_tc, bc_dirichlet=None, bc_neumann=None,
                 AllDirichlet=None, Dirichlet=None, metric=None, solverInfo=None):
        """Creates an poisson PDE solver. arguments are the same as pigasus.__init__

        geometry:
            The geometry must be an object cad_geometry.

        list_tc:
            contains the testcases list, needed for each implicit/explicit operator. The user must first specify the explicit terms, then the implicit ones.

        Returns:
           A PDE object.

        .. note::

            See also: :ref:`fem.basicPDE`.

        """

        # ...
        self.dim = geometry.dim
        self.nPDEs = len(list_tc)
        # ...

        # ...
        self.list_PDE = []
        for tc in list_tc:
            testcase = {}

            try:
                testcase['A'] = tc['A']
            except:
                pass

            try:
                testcase['v'] = tc['v']
            except:
                pass

            try:
                testcase['w'] = tc['w']
            except:
                pass

            try:
                testcase['b'] = tc['b']
            except:
                pass

            try:
                testcase['D2'] = tc['D2']
            except:
                pass

            try:
                testcase['u'] = tc['u']
            except:
                pass

            try:
                testcase['f'] = tc['f']
            except:
                pass

            if bc_dirichlet is not None:
                testcase['bc_dirichlet'] = bc_dirichlet

            if bc_neumann is not None:
                testcase['bc_neumann'] = bc_neumann

            if AllDirichlet is not None:
                testcase['AllDirichlet']  = AllDirichlet

            if Dirichlet is not None:
                testcase['Dirichlet']  = Dirichlet

            if metric is not None:
                testcase['metric']  = metric

            if solverInfo is not None:
                testcase['solverInfo']  = solverInfo

            # ...
            PDE = basicPDE(geometry=geometry, testcase=testcase)
            self.list_PDE.append(PDE)
            # ...
        # ...
    #-----------------------------------

    #-----------------------------------
    def __del__(self):
        for PDE in self.list_PDE:
            PDE.__del__()
    #-----------------------------------

    #-----------------------------------
    def free(self):
        for PDE in self.list_PDE:
            PDE.free()
    #-----------------------------------

    #-----------------------------------
    def initialize(self, list_u0=None):
        _list_u0        = [None]*self.nPDEs
        if list_u0 is not None:
            _list_u0    = list_u0

        for u0,PDE in zip(_list_u0,self.list_PDE):
            U = PDE.unknown
            if u0 is None:
                U.set(np.zeros(U.size))
            else:
                PDE.interpolate(u0, field=U)
    #-----------------------------------

    #-----------------------------------
    def assembly(self, list_f=None, list_update=None):
        _list_f         = [None]*self.nPDEs
        if list_f is not None:
            _list_f     = list_f

        _list_update    = [False]*self.nPDEs
        if list_update is not None:
            _list_update = list_update

        for f,update,PDE in zip(_list_f, _list_update, self.list_PDE):
            if f is not None:
                self.F_V.set_func(f)

            PDE.forceAssembly = update

            if PDE.forceAssembly or not PDE.Assembled:
                PDE.assembly()
    #-----------------------------------

    #-----------------------------------
    def plot(self):
        PDE = self.list_PDE[-1]
        PDE.plot()
    #-----------------------------------

if __name__ == '__main__':
    import caid.cad_geometry  as cg
    import pylab                as pl
    import numpy                as np
    from .utils import assembler

    # ...
    sin = np.sin ; cos = np.cos ; exp = np.exp ; sqrt = np.sqrt ; pi = np.pi
    # ...

    def testcase_line():
        # ...
        testcase = {}

    #    testcase['mass']       = lambda x : [0.]
    #    testcase['adv']        = lambda x : [0.]
        testcase['stiffness']  = lambda x : [1.]

        kx = 2. * pi

        # ...
        # exact solution
        # ...
        u = lambda x : sin ( kx * x )
        testcase['u'] = u
        # ...

        # ...
        # rhs
        # ...
        f = lambda x : ( kx**2) * sin ( kx * x )
        testcase['f'] = f
        # ...

        # ...
        # values of u at the boundary
        # ...
        g1 = lambda x: u(x)
        g2 = lambda x: u(x)

        testcase['list_faces'] = [[1,2]]
        testcase['list_g'] = [[g1, g2]]
        # ...

        return testcase

    def testcase_square_Dirichlet():
        # ...
        testcase = {}

    #    testcase['mass']       = lambda x,y : [0.]
    #    testcase['adv']        = lambda x,y : [0., 0.]
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
        # values of u at the boundary
        # ...
        g1 = lambda x,y : u(x,y)
        g2 = lambda x,y : u(x,y)
        g3 = lambda x,y : u(x,y)
        g4 = lambda x,y : u(x,y)

        testcase['list_faces'] = [[1,2,3,4]]
        testcase['list_g'] = [[g1, g2, g3, g4]]
        # ...

        return testcase

    def testcase_annulus_Dirichlet():
        # ...
        testcase = {}

    #    testcase['mass']       = lambda x,y : [0.]
    #    testcase['adv']        = lambda x,y : [0., 0.]
        testcase['stiffness']  = lambda x,y : [1., 0., 0., 1.]

        kx = 2. * pi
        ky = 2. * pi

        rmin = 0.5 ; rmax = 1.0

        # ...
        # exact solution
        # ...
        u = lambda x,y : [sin ( kx * ( x**2 + y**2 - rmin**2 ) / ( rmax**2 - rmin**2 ) )]
        testcase['u'] = u
        # ...

        # ...
        # rhs
        # ...
        def f ( x , y ) :
            r2 = x**2  + y**2
            w  = ( r2 - rmin**2 ) / ( rmax**2 - rmin**2 )
            c  = cos ( kx * w )
            s  = sin ( kx * w )
            coef = 2.0 * kx / ( rmax**2 - rmin**2 )
            return [- ( 2.0 * coef * c - coef**2 * r2 * s )]
        testcase['f'] = f
        # ...

        # ...
        # values of u at the boundary
        # ...
        g1 = lambda x,y : u(x,y)
        g2 = lambda x,y : u(x,y)
        g3 = lambda x,y : u(x,y)
        g4 = lambda x,y : u(x,y)

        testcase['list_faces'] = [[1,2,3,4]]
        testcase['list_g'] = [[g1, g2, g3, g4]]
        # ...

        return testcase

    # ...
    def func_n (x, y):
        list_nx = []
        list_ny = []
        for (u,v) in zip(x,y):
            nx = 0. ; ny = 0.
            if np.allclose(u, 0.):
                nx = -1.
                ny = 0.
            if np.allclose(u, 1.):
                nx = 1.
                ny = 0.

            if np.allclose(v, 0.):
                nx = 0.
                ny = -1.
            if np.allclose(v, 1.):
                nx = 0.
                ny = 1.

            list_nx.append(nx)
            list_ny.append(ny)

        nx = np.asarray(list_nx)
        ny = np.asarray(list_ny)

        return nx, ny
    # ...

    def testcase_square_Neumann():
        # ...
        testcase = {}
        testcase['Neumann'] = True
        testcase['list_DirFaces'] = [[]]

    #    testcase['mass']       = lambda x,y : [0.]
    #    testcase['adv']        = lambda x,y : [0., 0.]
        testcase['stiffness']  = lambda x,y : [1., 0., 0., 1.]

        xc = 0.5
        yc = 0.5

        kx = 1. * pi
        ky = 1. * pi
        # ...
        # exact solution
        # ...
        u = lambda x,y : sin(kx*(x-xc)) * sin(ky*(y-yc))
        testcase['u'] = u
        # ...

        # ...
        # rhs
        # ...
        f = lambda x,y : (kx**2 + ky**2) * sin(kx*(x-xc)) * sin(ky*(y-yc))
        testcase['f'] = f
        # ...

        # ...
        # values of gradu.n at the boundary
        # ...
        gradu   = lambda x,y : [  kx * cos(kx*(x-xc)) * sin(ky*(y-yc)) \
                                    , ky * sin(kx*(x-xc)) * cos(ky*(y-yc)) ]

        def g (x,y) :
            du  = gradu (x, y)
            nx, ny = func_n (x, y)

            return nx * du[0] + ny * du[1]

        testcase['g'] = g
        # ...

        return testcase
    #-----------------------------------

    #-----------------------------------
    # SINGLE PRECISION
    atol    = 1e-8
    rtol    = 1e-7
    # DOUBLE PRECISION
    #atol    = 1e-15
    #rtol    = 1e-14

    niter 	= 5000
    nx      = 15
    ny      = 15
    px      = 2
    py      = 2
    #-----------------------------------

    #-----------------------------------
    # ...
    from caid.cad_geometry import line
    geo1 = line(n=[nx], p=[px])

    testcase1 = testcase_line()
    PDE1 = basicPDE(geometry=geo1, testcase=testcase1)
    # ...

    # ...
#    testcase2 = testcase_square_Neumann()
#    from caid.cad_geometry import square as domain
#    geo2 = domain(n=[nx,ny], p=[px,py])

    testcase2 = testcase_square_Dirichlet()
    from caid.cad_geometry import square as domain
    geo2 = domain(n=[nx,ny], p=[px,py])
##    PDE1 = basicPDE(geometry=geo2, testcase=testcase2)

#    from caid.cad_geometry import annulus as domain
#    geo2 = domain(n=[nx+1,ny+1], p=[px,py])
#    testcase2 = testcase_annulus_Dirichlet()

    PDE2 = basicPDE(geometry=geo2, testcase=testcase2)
    # ...

    # ...
    MyAssembler = assembler()
    # ...

    # ...
    PDE1.assembly()
    PDE2.assembly()
#    MyAssembler.assembly([PDE1, PDE2])
    # ...

    # ...
    from scipy.io import mmwrite
    PDE1.solve()
    PDE2.solve()
    # ...

    # ...
    PDE1.plot()  ; pl.show()
    PDE2.plot()  ; pl.show()
    # ...

    # ...
    normU1 = PDE1.norm()
    print("norm U-1D   = ", normU1)
    # ...

    # ...
    normU2 = PDE2.norm()
    print("norm U-2D   = ", normU2)
    # ...



    # ...
#    U1 = PDE1.get()
#    u1 = U1.get()       #; print u1
#    c1 = U1.tomatrix(0) #; print c1
#
#    U2 = PDE2.get()
#    u2 = U2.get()       #; print u2
#    c2 = U2.tomatrix(0) #; print c2
    # ...
