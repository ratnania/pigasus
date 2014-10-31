# -*- coding: UTF-8 -*-
#! /usr/bin/python

# To change this template, choose Tools | Templates
# and open the template in the editor.

__author__="ARA"
__all__ = ['fem']
__date__ ="$Jan 11, 2012 6:43:48 PM$"

from . import constants as _cst
import numpy as _np
from .pigasusObject import *

from .common_obj import _singleton
@_singleton
class fem(pigasusObject):

    def __init__(self):
        pigasusObject.__init__(self)

        self.stdoutput      = True
#        self.detail         = 3
        self.detail         = 0
        self.Initialized    = False
        self.IJVInitialized = False
        self.deleted        = False

#        print "==== CREATE NEW FEM ===="

    def _set_nvalues(self):
        """
        Global initialization. Must be called after all fields, matrices, norms declarations
        """
        # ...
        self.com.pyfem.set_nspaces(self.com.nspaces)
        self.com.pyfem.set_nmappings(self.com.nmappings)
        self.com.pyfem.set_nmetrics(self.com.nmetrics)
        self.com.pyfem.set_ngraphs(self.com.ngraphs)
        self.com.pyfem.set_nmatrices(self.com.nmatrices)
        self.com.pyfem.set_noperators(self.com.noperators)
        self.com.pyfem.set_nfields(self.com.nfields)
        self.com.pyfem.set_ngrids(self.com.ngrids)
        self.com.pyfem.set_nnorms(self.com.nnorms)
        self.com.pyfem.set_nsolvers(self.com.nsolvers)
        # ...

        # ...
        try:
            li_maxnpatchs = max( [G.npatchs for G in self.com.grids])
        except:
            li_maxnpatchs = 0
        self.com.pyfem.set_maxnpatchs(li_maxnpatchs)
        # ...

        # ...
        try:
            li_maxnparams_operators = max( [O.nparam for O in self.com.operators])
        except:
            li_maxnparams_operators = 0

        self.com.pyfem.set_maxnparams_operators(li_maxnparams_operators)
        # ...

        # ...
        try:
            li_maxnparams_fields = max( [F.nparam for F in self.com.fields])
        except:
            li_maxnparams_fields = 0

        self.com.pyfem.set_maxnparams_fields(li_maxnparams_fields)
        # ...

        # ...
        try:
            li_maxnaddto = max( [len(O.list_addto) for O in self.com.operators])
        except:
            li_maxnaddto = 0

        self.com.pyfem.set_maxnaddto(li_maxnaddto)
        # ...

        # ...
        try:
            li_maxnparams_norms = max( [N.nparam for N in self.com.norms])
        except:
            li_maxnparams_norms = 0

        self.com.pyfem.set_maxnparams_norms(li_maxnparams_norms)
        # ...

        # ...
        try:
            li_maxder_field_op = max([F.nderiv for F in self.com.fields])
        except:
            li_maxder_field_op = 0
        self.com.pyfem.set_maxder_field_op(li_maxder_field_op)
        # ...

    def _set_infoFields(self):
        # ...
        # getting informations for each field
        # ...
        for F in self.com.fields:
            li_parameval = 0
            if F.parameval :
                li_parameval = 1
            self.com.pyfem.set_field(F.id, F.ndof, F.space.id, F.size, F.loc_id \
            , F.type, F.operator, F.operande_id, li_parameval, F.nparam)
        # ...

    def _set_infoNorms(self):
        # ...
        # getting informations for each norm
        # ...
        for N in self.com.norms:
            self.com.pyfem.set_norm(N.id, N.field.id, N.type, N.loc_id)
        # ...

    def _set_infoOperators(self):
        # ...
        # getting informations for each matrix
        # ...
        for O in self.com.operators:

            li_parameval = 0
            if O.parameval :
                li_parameval = 1

            li_transpose = 0
            if O.transpose :
                li_transpose = 1

            if O.spaces is None:
                space1_id = -1
                space2_id = -1
            else:
                space1_id = O.spaces[0].id
                space2_id = O.spaces[1].id

#            print "current operator id : ", O.id
#            print "params :", O.id, O.type, space1_id, space2_id, O.loc_id, O.nparam
            self.com.pyfem.set_operator(  O.id, O.type, space1_id, space2_id \
                                      , O.loc_id, O.nparam, li_parameval \
                                      , li_transpose)
        # ...

    def _set_infomatrices(self):
        # ...
        # getting informations for each matrix
        # ...
        for M in self.com.matrices:

            if M.graph is None:
                graph_id = -1
            else:
                graph_id = M.graph.id

            if M.IJVAssembly :
                IJVAssembly = 1
            else:
                IJVAssembly = 0

#            print "current matrix id : ", M.id
#            print "params :", M.id, M.graph.id
            self.com.pyfem.set_matrix(  M.id, M.type, graph_id, IJVAssembly)
        # ...

    def _set_infoSpacesPart1(self):
        # ...
        # getting informations for each space - part 1
        # ...
        for S in self.com.spaces:
            li_extmap = 0
            li_mapping_id = -1
            if S.exterior_map :
                li_extmap = 1
                li_mapping_id = S.metric.mapping.id

            li_storeddata = 0
            if S.storeddata :
                li_storeddata = 1

            if S.grids is None :
                li_grids_id = -1
            else:
                li_grids_id = S.grids.id

            li_composed_space = 0
            if S.type == "space_vect" :
                li_composed_space = 1

            li_composed_nspaces = 0
            if S.type == "space_vect" :
                li_composed_nspaces = len(S.spaces)

            self.com.pyfem.set_space(S.id, li_extmap, li_mapping_id, S.tensorlevel \
            , li_grids_id, li_storeddata, S.ndof, S.size, S.maxnen    \
            , li_composed_space, li_composed_nspaces)
        # ...

    def _set_infoGraphs(self):
        # ...
        # getting informations for each graph
        # ...
        for G in self.com.graphs:
            self.com.pyfem.set_graph(G.id, G.spaces[0].id, G.spaces[1].id)
        # ...

    def _set_infoSolvers(self):
        # ...
        # getting informations for each solver
        # ...
        for S in self.com.solvers:
            if S.matrix is None :
                li_matrix_id = -1
            else:
                li_matrix_id = S.matrix.id

            if S.residual is None :
                li_residual_id = -1
            else:
                li_residual_id = S.residual.id

            li_solver  = S.solver
            if li_solver is None:
                li_solver = -1

            self.com.pyfem.set_solver(S.id, li_matrix_id, li_residual_id, li_solver)
        # ...

    def _set_infoGridsPart1(self):
        # ...
        # getting informations for each grids
        # ...
        for G in self.com.grids:

            li_usemetric = 0
            li_metric_id = -1
            if G.usemetric:
                li_usemetric = 1
                li_metric_id = G.metric_id

            li_nfields = len(G.get_fields_id())
            li_noperators = len(G.get_operators_id())
            li_nnorms = len(G.get_norms_id())
            self.com.pyfem.set_grids(G.id, G.npatchs, G.dim, G.Rd, G.ndof \
            , li_nfields, li_noperators, li_usemetric, li_metric_id, li_nnorms)
#            print "G.id, G.npatchs, G.dim, G.Rd, G.ndof, li_nfields, li_nmatrices, li_usemetric, li_metric_id", \
#            G.id, G.npatchs, G.dim, G.Rd, G.ndof, li_nfields, li_nmatrices, li_usemetric, li_metric_id
        # ...

    def _set_infoSpacesPart2(self):
        # ...
        # getting informations for each space - part 2
        # ...
        for V in self.com.spaces:
            if V.exterior_map :
                M = self.com.mappings[V.metric.mapping.id]

                li_storeddata = 0
                if M.storeddata :
                    li_storeddata = 1
                self.com.pyfem.set_mapping(M.id, V.tensorlevel, li_storeddata, V.id)
#                print "M.id, li_tensor, V.grids.npatchs, V.grids.id, li_storeddata : " \
#                , M.id, li_tensor, V.grids.npatchs, V.grids.id, li_storeddata
        # ...

    def _set_infoSpacesPart3(self):
        # ...
        # getting informations for each space - part 3
        # ...
        for V in self.com.spaces:
            if V.type == "space_vect" :
                li_position = 0
                for loc_V in V.spaces:
                    li_position += 1
                    self.com.pyfem.pyfem_set_idspaces_composedspace(V.id,li_position,loc_V.id)

            V.add_patchs()
            V.add_connectivity()
            if V.exterior_map :
                M = self.com.mappings[V.metric.mapping.id]
                M.add_patchs()
        # ...

    def _set_infoGridsPart2(self):
        # ...
        # getting informations for each grids - part 2
        # ...
        for G in self.com.grids:
#            print "==============="
#            print " Current Grid is ", G.id
#            print "==============="
#            print "G.maxnpts="+str(G.maxnpts)
#            print "G.dirmaxnpts="+str(G.dirmaxnpts)
            for li_id in range(0,G.npatchs):
#                print "------"
#                print " Current patch is", li_id
#                print "------"
#                print "G.list_grid[li_id].maxnpts="+str(G.list_grid[li_id].maxnpts)
#                print "G.list_grid[li_id].dirmaxnpts="+str(G.list_grid[li_id].dirmaxnpts)
#                print "G.list_grid[li_id].nel="+str(G.list_grid[li_id].nel)

                self.com.pyfem.set_patch(G.id, li_id    \
                , G.list_grid[li_id].nel    \
                , G.list_grid[li_id].maxnpts \
                , G.list_grid[li_id].tensorlevel \
                , G.list_grid[li_id].dirmaxnpts )
        # ...

    def _set_infoGridsPart3(self):
        # ...
        # getting informations for each grids - part 3
        # ...
        for G in self.com.grids:
            G.save_grids()
        # ...

    def initialize(self, dbasis=None, stdoutput=None, detail=None):
        if self.Initialized:
            return
        # ...
        if detail is not None:
            self.detail = detail
        if stdoutput is not None:
            self.stdoutput = stdoutput
        li_stdoutput = 0
        if not self.stdoutput :
            li_stdoutput = 1
        self.com.pyfem.init(li_stdoutput,self.detail)
        # ...
        # ...
        self._set_nvalues()
        # ...
        # ...
#        print "part 1"
        self.com.pyfem.create_partone()
#        print 'done.'
        # ...

        # ...
        self._set_infoFields()
        self._set_infoNorms()
        self._set_infoOperators()
        self._set_infomatrices()
        self._set_infoSpacesPart1()
        self._set_infoGraphs()
        self._set_infoSolvers()
        self._set_infoGridsPart1()
        self._set_infoSpacesPart2()
        # ...

        # ...
#        print "part 2"
        self.com.pyfem.create_parttwo()
#        print 'done.'
        # ...

        # ...
        self._set_infoSpacesPart3()
        self._set_infoGridsPart2()
        # ...

        # ...
#        print "part 3"
        self.com.pyfem.create_partthree()
#        print 'done.'
        # ...

        # ...
        self._set_infoGridsPart3()
        # ...

        # ...
#        print "part 4"
        self.com.pyfem.create_partfour()
#        print 'done.'
        # ...

        # ...
        # generate dbasis for each grids
        # ...
#        maxnderiv_m = max([M.nderiv for M in self.com.matrices])
#        maxnderiv_f = max([F.nderiv for F in self.com.fields])
#        maxnderiv = max([maxnderiv_f, maxnderiv_m])
#        for G in self.com.grids:
#            G.gen_dbasis(dbasis=dbasis, nderiv=maxnderiv)

        # ...

        # ...
        self.setInfoData()
        for V in self.com.spaces:
            V.setInfoData()
        for O in self.com.operators:
            O.setInfoData()
        for M in self.com.matrices:
            M.setInfoData()
        for F in self.com.fields:
            F.setInfoData()
        for N in self.com.norms:
            N.setInfoData()
#        for G in self.com.grids:
#            G.setInfoData()
        # ...
        # ...
        self.Initialized = True
        self.deleted     = False
        # ...

    def InitializeIJVAssembly(self):
        # ...
        # do this only if we use ijv assembly on this space
        for V in self.com.spaces:
            IJVAssembly = False
            for M in self.com.matrices:
                if  M.IJVAssembly:
                    if (M.graph.spaces[0].id == V.id) or (M.graph.spaces[1].id == V.id):
                        IJVAssembly = True
            if IJVAssembly:
                V.EM.InitializeIDsupports(V.connectivity.ID_loc, V.connectivity.ID)
        # ...
        self.IJVInitialized = True

    def __del__(self):
        if not self.deleted :
            if self.com.usefiga:
                self.com.pyfem.free_figa()
            self.com.pyfem.pyfem_free()
            self.Initialized = False
            self.com.reset()

            self.deleted = True

    def setInfoData(self):
        """
        prints informations about the current matrix
        """
        self.infoData['nfields'] = str(self.com.nfields)
        self.infoData['nmatrices'] = str(self.com.nmatrices)
        self.infoData['noperators'] = str(self.com.noperators)
        self.infoData['nspaces'] = str(self.com.nspaces)
        self.infoData['nmappings'] = str(self.com.nmappings)
        self.infoData['ngrids'] = str(self.com.ngrids)

    def reset_toassembly(self):
        self.com.pyfem.reset_operators_toassembly()
        self.com.pyfem.reset_matrices_toassembly()
        self.com.pyfem.reset_fields_toassembly()
        self.com.pyfem.reset_norms_toassembly()
        self.com.pyfem.reset_spaces_toassembly()

    def _set_patch_toassembly(self, grids_id, list_patchs):
        for P in list_patchs:
            self.com.pyfem.set_patch_toassembly(grids_id, P)
        self.com.pyfem.pyfem_set_patchs_toassembly(grids_id)

    def _set_elts_toassembly(self, grids_id, list_elts):
        self.com.pyfem.pyfem_set_elts_toassembly(grids_id, list_elts)

    def _set_matrix_toassembly(self, list_matrices):
        for M in list_matrices:
            self.com.pyfem.pyfem_reset_matrix(M.id)
            self.com.pyfem.set_matrix_toassembly(M.id)

    def _set_operator_matrices_toassembly(self, list_operators, list_matrices):
        for O in list_operators:
            list_OM = [data[0] for data in O.list_addto]
            list_M = [M.id for M in list_matrices if M in list_OM]
            self.com.pyfem.pyfem_set_operator_matrices_toassembly( O.id, list_M, len(list_M))

    def _set_field_toassembly(self, list_fields):
        for F in list_fields:
            self.com.pyfem.set_field_toassembly(F.id)

    def _set_norm_toassembly(self, list_norms):
        for N in list_norms:
            self.com.pyfem.set_norm_toassembly(N.id)

    def _set_space_toassembly(self,list_spaces):
        for V in list_spaces:
            self.com.pyfem.set_space_toassembly(V.id)

    def assembly(self, matrices=[], fields=[], norms=[], patchs=None, elts=None):
        """
        This assembles the desired fields, matrices and norms over the computational grid
        """
        # ...
        if not self.Initialized:
            self.initialize()
        # ...

#        from time import clock
#        start = clock()
#        self.evalfunc(matrices, fields, norms, patchs, elts)
#        elapsed = (clock() - start)
#        print ("CPU time for evaluating Fields and Matrices is :" + str (elapsed))
#        start = clock()

#        self.reset_toassembly()

        operators_id = []
        for M in matrices:
            operators_id += [O.id for O in M.operators]
        operators_id = list(_np.unique(_np.asarray(operators_id)))
        operators = [self.com.operators[id] for id in operators_id]

        list_grids_id_o = [self.com.spaces[O.space.id].grids.id for O in operators]
        list_grids_id_f = [self.com.spaces[F.space.id].grids.id for F in fields]
        list_grids_id_n = [self.com.spaces[N.field.space.id].grids.id for N in norms]
        list_grids_id   = list(set(list_grids_id_o) | set(list_grids_id_f) | set(list_grids_id_n))
        list_grids_id   = list(_np.unique(_np.asarray(list_grids_id)))

        for grids_id in list_grids_id:
            # ...
            # ...

            Gr = self.com.grids[grids_id]
            self.com.pyfem.reset_patchs_toassembly(grids_id)

            Gr.set_patchs_toassembly(patchs)
            Gr.set_elts_toassembly(elts)

            list_patchs     = Gr.get_patchs_toassembly()
            list_elts       = Gr.get_elts_toassembly()
            list_fields     = [F for F in fields if self.com.spaces[F.space.id].grids.id==grids_id]
            list_norms      = [N for N in norms if self.com.spaces[N.field.space.id].grids.id==grids_id]
            list_operators  = [O for O in operators if self.com.spaces[O.space.id].grids.id==grids_id]
            list_matrices = []
            for M in matrices:
                addit = False
                # TODO a optimiser
                for O in M.operators:
                    if O in list_operators:
                        addit = True
                if addit :
                    list_matrices.append(M)

            list_spaces_o   = [O.spaces[0].id for O in list_operators] + [O.spaces[1].id for O in list_operators]
            list_spaces_f   = [F.space.id for F in list_fields]
            list_spaces_n   = [N.field.space.id for N in list_norms]
            list_all_spaces = list(set(list_spaces_o) | set(list_spaces_f) | set(list_spaces_n))
            list_all_spaces = list(_np.unique(_np.asarray(list_all_spaces)))
            list_spaces     = [self.com.spaces[V_id] for V_id in list_all_spaces]

#            self.evalfunc(list_operators, list_fields, list_norms, list_patchs, list_elts)

#            print "---------------------------------"
#            print "Current Grid : ", grids_id
#            print "list_patchs : ", list_patchs
#            print "list_elts : ", list_elts
#            print "list_matrices : ", [M.id for M in list_matrices]
#            print "list_operators : ", [O.id for O in list_operators]
#            print "list_fields : ", [M.id for M in list_fields]
#            print "list_norms : ", [M.id for M in list_norms]
#            print "list_spaces : ", [M.id for M in list_spaces]
#            print "---------------------------------"

            self._set_patch_toassembly(grids_id, list_patchs)
            self._set_elts_toassembly(Gr.id, list_elts)
            self._set_matrix_toassembly(list_matrices)
            self._set_operator_matrices_toassembly(list_operators, list_matrices)
            self._set_field_toassembly(list_fields)
            self._set_norm_toassembly(list_norms)
            self._set_space_toassembly(list_spaces)

#            self.com.pyfem.pyfem_save_terms_toassembly()
            list_spaces_id = [S.id for S in list_spaces]
            self.evalfunc(list_operators, list_fields, list_norms \
                          , list_spaces_id, list_patchs, list_elts)
#            self.com.pyfem.pyfem_load_terms_toassembly()

            self.com.pyfem.pyfem_assembly(grids_id)
            self.reset_toassembly()



#        elapsed = (clock() - start)
#        print ("CPU time for assembling Fields and Matrices is :" + str (elapsed))

# ----------------------------------------------------------------------------------
    def assembly_nodelist(self, matrices=[], fields=[], norms=[], nodelist=None):
        """
        This assembles the desired fields, matrices and norms over the computational grid
        we suppose here that we use one space
        """
        # TODO : upgrade to nbr of spaces > 1

        # ...
        if not self.Initialized:
            self.initialize()
        # ...
        # ...
        if not self.IJVInitialized:
            self.InitializeIJVAssembly()
        # ...

        V = self.com.spaces[0]
        V.EM.setNodeList(nodelist)
        list_patchs = V.EM.GetLocalPatchs()
        for ipatch in list_patchs:
            list_elts = V.EM.GetLocalElements(ipatch)
#            print "-----"
#            print 'patch : ', ipatch
#            print 'elts   : ', list_elts
            self.assembly(  matrices=matrices, fields=fields, norms=norms \
                          , patchs=[ipatch], elts=list_elts)
# ----------------------------------------------------------------------------------


    def _find_fields_sp(self, fields=[]):
        if len(fields) == 0:
            return [], []
        # we store the spaces concerned by the evaluation
        list_spaces_id_f = []
        # we store for each space, the list of related fields
        list_fields_sp = []

        for S in self.com.spaces:
            list_fields_sp.append([])

        # correspondance intialization for fields
        for F in fields:
            if F.space.id not in list_spaces_id_f:
                list_spaces_id_f.append(F.space.id)
            list_fields_sp[F.space.id].append(F.id)

        return list_fields_sp, list_spaces_id_f

    def _find_operators_sp(self, operators=[]):
        if len(operators) == 0:
            return [], []
        # we store the spaces concerned by the evaluation
        list_spaces_id_o = []
        # we store for each space, the list of related operators
        list_operators_sp= []

        for S in self.com.spaces:
            list_operators_sp.append([])

        # correspondance intialization for operators
        for O in operators:
            if O.space.id not in list_spaces_id_o:
                list_spaces_id_o.append(O.space.id)
            list_operators_sp[O.space.id].append(O.id)

        return list_operators_sp, list_spaces_id_o

    def _find_norms_sp(self, norms=[]):
        if len(norms) == 0:
            return [], []
        # we store the spaces concerned by the evaluation
        list_spaces_id_n = []
        # we store for each space, the list of related norms
        list_norms_sp = []

        for S in self.com.spaces:
            list_norms_sp.append([])

        # correspondance intialization for norms
        for N in norms:
            if N.field.space.id not in list_spaces_id_n:
                list_spaces_id_n.append(N.field.space.id)
            list_norms_sp[N.field.space.id].append(N.id)

        return list_norms_sp, list_spaces_id_n


    def evalfunc(self, operators=[], fields=[], norms=[], list_spaces=None, patchs=None, elts=None):
        """
        Evaluation of the param-functions for the given fields, norms and operators over the computational grid
        """
        list_operators=[[]]
        list_fields=[[]]
        list_norms=[[]]
        ll_store = True

        list_fields_sp, list_spaces_id_f    = self._find_fields_sp(fields)
        list_operators_sp, list_spaces_id_o = self._find_operators_sp(operators)
        list_norms_sp, list_spaces_id_n     = self._find_norms_sp(norms)

        # find the union of the two lists list_spaces_id_o and list_spaces_id_f
        list_spaces_id = list(set(list_spaces_id_o) | set(list_spaces_id_f) | set(list_spaces_id_n))
        if list_spaces is None:
            list_spaces = list_spaces_id

#        print "lists : ", list_spaces_id_o, list_spaces_id_f, list_spaces_id_n

        for li_space in list_spaces:
#            print "current space : ", li_space
            V = self.com.spaces[li_space]
            grids_id = self.com.spaces[li_space].grids.id
            Gr = self.com.grids[grids_id]
            self.com.pyfem.reset_patchs_toassembly(grids_id)

            Gr.set_patchs_toassembly(patchs)
            Gr.set_elts_toassembly(elts)

            list_patchs = Gr.get_patchs_toassembly()
            list_elts = Gr.get_elts_toassembly()

            for P in list_patchs:
                self.com.pyfem.set_patch_toassembly(grids_id, P)

            self.com.pyfem.pyfem_set_elts_toassembly(Gr.id, list_elts)

            li_npatchs = V.grids.npatchs

            for li_patch in range(0,li_npatchs):
                # ...
                # get the grid
                # ...
                _G = V.grids.list_grid[li_patch]
                shape_pts = [_G.nel, _G.npts]
                shape_val = [0]+list(shape_pts)

                patch = V.geometry[li_patch]
                V.set_currentPatch(patch)

                V.computePoints(li_patch)
                lpr_parampts = V.get_sites()
#                print "---------"
#                print 'patch-id : ', li_patch
#                print "parampts : ", lpr_parampts
                lpr_pts      = V.get_points()
#                    print "pts : ", lpr_pts
                xyz = []
                for i in range(0, V.dim):
                    xyz.append(lpr_pts[i,0,:])
                lpr_pts = xyz
                # ...

                # ...
                # Fields treatment
                # ...
                if li_space in list_spaces_id_f:
                    for li_field in list_fields_sp[li_space]:
                        F = self.com.fields[li_field]
#                        print "current field :", li_field
#                        print "F.type =", F.type
#                        print "F.operator =", F.operator
                        list_val = F.evalfunc(li_patch, lpr_pts, elts)
#                        print list_val[0]
                        shape_val[0] = len(list_val)
                        lpr_val = _np.zeros(shape_val)
                        for i in range(0, shape_val[0]):
                            lpr_val[i,:,:] = list_val[i].reshape(shape_val[1:])
                        if ll_store :
                            self.com.pyfem.set_field_on_grids  (li_field, li_patch, lpr_val)
                # ...

                # ...
                # Operators treatment
                # ...
                if li_space in list_spaces_id_o:
                    for li_operator in list_operators_sp[li_space]:
#                        print "current opoerator :", li_operator
                        O = self.com.operators[li_operator]
                        list_val = O.evalfunc(li_patch, lpr_pts, elts)
                        shape_val[0] = len(list_val)
#                        print list_val
#                        print shape_val
                        lpr_val = _np.zeros(shape_val)
                        for i in range(0, shape_val[0]):
                            lpr_val[i,:,:] = list_val[i].reshape(shape_val[1:])
#                            print lpr_val
                        if ll_store :
                            self.com.pyfem.set_operator_on_grids  (li_operator, li_patch, lpr_val)
#                            print "set_operator_on_grids: done."
                # ...

                # ...
                # Norms treatment
                # ...
                if li_space in list_spaces_id_n:
                    for li_norm in list_norms_sp[li_space]:
#                        print "current norm :", li_norm
                        N = self.com.norms[li_norm]
                        # ...
                        # evaluation of the param function
                        # ...
                        list_val = N.evalfunc(li_patch, lpr_pts, elts)
                        shape_val[0] = len(list_val)
                        lpr_val = _np.zeros(shape_val)
                        for i in range(0, shape_val[0]):
                            lpr_val[i,:,:] = list_val[i].reshape(shape_val[1:])
                        if ll_store :
                            self.com.pyfem.set_norm_on_grids  (li_norm, li_patch, lpr_val)
                        # ...
                        # ...
                        # evaluation of the exact solution
                        # ...
                        list_val = N.evalfunc(li_patch, lpr_pts, elts, type="exact")
                        shape_val[0] = len(list_val)
                        lpr_val = _np.zeros(shape_val)
                        for i in range(0, shape_val[0]):
                            lpr_val[i,:,:] = list_val[i].reshape(shape_val[1:])
                        if ll_store :
                            self.com.pyfem.set_field_on_grids  (N.field.id, li_patch, lpr_val)
                        # ...
                # ...

            # must reinitialize the current patch
            V.set_currentPatch(V.geometry[0])

#    def evalfields(self, fields, list_operator, apr_points=None):
#        """
#        Evaluates the fields operators over the computational grid
#        """
#        list_fields_sp, list_spaces_id_f = self._find_fields_sp(fields)
#        list_spaces_id = list_spaces_id_f
#        # ...
#        self.com.pyfem.pyfem_save_terms_toassembly()
#        # ...
#
#        list_values = [[]]*len(fields)
#
#        for li_space in list_spaces_id:
##            print "current space :", li_space
#            self.com.pyfem.pyfem_reset_terms_toassembly(self.com.spaces[li_space].grids.id)
#            if self.com.spaces[li_space].type == "space":
#                li_npatchs = self.com.spaces[li_space].grids.npatchs
#            else :
#                # the composed case
#                li_npatchs = self.com.spaces[li_space].spaces[0].grids.npatchs
#
#            for li_patch in range(0,li_npatchs):
#                if apr_points is None:
#                    lpr_pts_new = self.com.spaces[li_space].get_points_advanced(ai_patch_id=li_patch)
##                    lpr_pts = _np.zeros(lpr_pts_new.shape[:-1])
##                    lpr_pts[:,:,:] = lpr_pts_new[:,:,:, 0]
#                    lpi_shape = lpr_pts_new.shape
#                    nprod = _np.asarray(lpi_shape[:-1]).prod()
#                    lpr_pts = lpr_pts_new[:,:,:,0].reshape(nprod)
#                else :
#                    lpr_pts = apr_points
#
#                # ...
#                if li_space in list_spaces_id_f:
#                    for li_field in list_fields_sp[li_space]:
##                        print "current field :", li_field
#                        F = self.com.fields[li_field]
##                        print "F.type =", F.type
##                        print "F.operator =", F.operator
#                        # ...
#                        # first we change the field's type and operator
#                        # ...
#                        if F.type == _cst.FIELD_OPERATOR:
#                            li_type = F.type
#                            li_operande_id = F.operande_id
#                            li_operator = F.operator
#                            li_nparam = F.nparam
#                        else:
#                            li_type = _cst.FIELD_OPERATOR
#                            li_operande_id = F.id
#                            li_operator = list_operator[[FF.id for FF in fields].index(F.id)]
#                            li_nparam = F._compute_nparam(li_type, li_operator)
#
#                        li_parameval = 0
#                        if F.parameval :
#                            li_parameval = 1
#                        self.com.pyfem.set_field(F.id, F.ndof, F.space.id, F.size, F.loc_id \
#                        , li_type, li_operator, li_operande_id, li_parameval, li_nparam)
#                        # ...
#
#                        # ...
#                        # evaluation of the param function on the grid and storage into the elements_grid
#                        # ...
#                        lpr_val = F.evalfunc(li_patch, lpr_pts)
#                        self.com.pyfem.set_field_on_grids  (F.id, li_patch, lpr_val)
#                        # ...
#
#                        # ...
#                        # assemblage : field's evaluation on the grid
#                        # ...
#                        self.com.pyfem.set_field_toassembly(F.id)
#                        # ...
#
#                    # ...
#
#                    # ...
#                    self.com.pyfem.pyfem_assembly(self.com.spaces[li_space].grids.id)
#                    # ...
#
#                    # ...
#                    for li_field in list_fields_sp[li_space]:
##                        print "current field :", li_field
#                        F = self.com.fields[li_field]
#
#                        V = self.com.spaces[F.space.id]
#                        li_nel = V.grids.list_grid[li_patch].nel
#                        li_maxnpts = V.grids.list_grid[li_patch].maxnpts
#                        li_maxder = (V.dim ** MAP_FIELD_DIMOUT[li_operator])
#                        li_k = li_maxder * F.ndof
#                        lpr_values =  self.com.pyfem.pyfem_get_fieldh ( F.id  \
#                        , li_patch, li_k, F.ndof, li_maxder, li_nel, li_maxnpts )
#
#                        list_values[[FF.id for FF in fields].index(F.id)] = lpr_values
#
#                        # ...
#                        # finally we change the field's type and operator
#                        # ...
#                        li_parameval = 0
#                        if F.parameval :
#                            li_parameval = 1
#                        self.com.pyfem.set_field(F.id, F.ndof, F.space.id, F.size, F.loc_id \
#                        , F.type, F.operator, F.operande_id, li_parameval, F.nparam)
#                        # ...
#                    # ...
#
#                    # ...
#                    self.com.pyfem.pyfem_load_terms_toassembly()
#                    # ...
#
#        return list_values

    # ...
    def write(self, filename):
        from xml.dom.minidom import Document
        # Create the minidom document
        doc = Document()

        # Create the <pyiga> base element
        ROOT_TAG = 'pigasus'
        rootElt = doc.createElement(ROOT_TAG)
        doc.appendChild(rootElt)

        list_obj = self.com.spaces \
                + self.com.fields  \
                + self.com.matrices  \
                + self.com.norms

        for obj in list_obj:
            obj.writeInfoDataXML(doc, rootElt)

        f = open(filename, 'wr')
        s = doc.toprettyxml()
        f.write(s)
        f.close()
    # ...
