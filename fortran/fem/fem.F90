!**************************************************
!
!                   FEM MODULE
!
!**************************************************
module fem_module
    use metric_module
    use grids_module
    use geometries_def
    use geometries_module
    use field_module
    use norm_module
    use graph_module
    use operator_module
    use matrix_module
    use spaces_module
    use connectivities_def
    use connectivity_module
    use solver_module
    use fem_def
    implicit none

    integer, parameter, private :: mi_dtllevel_base = 0

contains

    !---------------------------------------------------------------
    subroutine set_nspaces_fem(self, ai_nspaces)
        implicit none
        TYPE(FEM) :: self
        INTEGER :: ai_nspaces

        self % oi_nspaces = ai_nspaces
    end subroutine set_nspaces_fem
    !---------------------------------------------------------------
    subroutine set_nmappings_fem(self, ai_nmappings)
        implicit none
        TYPE(FEM) :: self
        INTEGER :: ai_nmappings

        self % oi_nmappings = ai_nmappings
    end subroutine set_nmappings_fem
    !---------------------------------------------------------------
    subroutine set_nmetrics_fem(self, ai_val)
        implicit none
        TYPE(FEM) :: self
        INTEGER :: ai_val

        self % oi_nmetrics = ai_val
    end subroutine set_nmetrics_fem
    !---------------------------------------------------------------
    subroutine set_nsolvers_fem(self, ai_n)
        implicit none
        TYPE(FEM) :: self
        INTEGER :: ai_n

        self % oi_nSolvers = ai_n
    end subroutine set_nsolvers_fem     
    !---------------------------------------------------------------
    subroutine set_ngraphs_fem(self, ai_n)
        implicit none
        TYPE(FEM) :: self
        INTEGER :: ai_n

        self % oi_nGraphs = ai_n
    end subroutine set_ngraphs_fem    
    !---------------------------------------------------------------
    subroutine set_noperators_fem(self, ai_n)
        implicit none
        TYPE(FEM) :: self
        INTEGER :: ai_n

        self % oi_noperators = ai_n
    end subroutine set_noperators_fem    
    !---------------------------------------------------------------
    subroutine set_nmatrices_fem(self, ai_nMatrices)
        implicit none
        TYPE(FEM) :: self
        INTEGER :: ai_nMatrices

        self % oi_nMatrices = ai_nMatrices
    end subroutine set_nmatrices_fem
    !---------------------------------------------------------------
    subroutine set_nFields_fem(self, ai_nFields)
        implicit none
        TYPE(FEM) :: self
        INTEGER :: ai_nFields

        self % oi_nFields = ai_nFields
    end subroutine set_nFields_fem
    !---------------------------------------------------------------
    subroutine set_ngrids_fem(self, ai_ngrids)
        implicit none
        type(FEM) :: self
        INTEGER, INTENT(INOUT) :: ai_ngrids
        ! LOCAL

        self % oi_nGrids = ai_ngrids
    end subroutine set_ngrids_fem
    !---------------------------------------------------------------
    subroutine set_nnorms_fem(self, ai_nnorms)
        implicit none
        TYPE(FEM) :: self
        INTEGER :: ai_nnorms

        self % oi_nNorms = ai_nnorms
    end subroutine set_nnorms_fem
    !---------------------------------------------------------------
    subroutine set_maxnpatchs_fem(self, ai_val)
        implicit none
        TYPE(FEM) :: self
        INTEGER, intent(in) :: ai_val

        self % oi_maxnpatchs = ai_val
    end subroutine set_maxnpatchs_fem
    !---------------------------------------------------------------
    subroutine set_maxnparams_operators_fem(self, ai_val)
        implicit none
        TYPE(FEM) :: self
        INTEGER, intent(in) :: ai_val

        self % oi_maxnparams_operators = ai_val
    end subroutine set_maxnparams_operators_fem
    !---------------------------------------------------------------
    subroutine set_maxnaddto_fem(self, ai_val)
        implicit none
        TYPE(FEM) :: self
        INTEGER, intent(in) :: ai_val

        self % oi_maxaddto = ai_val
    end subroutine set_maxnaddto_fem    
    !---------------------------------------------------------------
    subroutine set_maxnparams_fields_fem(self, ai_val)
        implicit none
        TYPE(FEM) :: self
        INTEGER, intent(in) :: ai_val

        self % oi_maxnparams_fields = ai_val
    end subroutine set_maxnparams_fields_fem
    !---------------------------------------------------------------
    subroutine set_maxnparams_norms_fem(self, ai_val)
        implicit none
        TYPE(FEM) :: self
        INTEGER, intent(in) :: ai_val

        self % oi_maxnparams_norms = ai_val
    end subroutine set_maxnparams_norms_fem
    !---------------------------------------------------------------
    subroutine set_maxder_field_op_fem(self, ai_val)
        implicit none
        TYPE(FEM) :: self
        INTEGER, intent(in) :: ai_val

        self % oi_maxder_field_op = ai_val
    end subroutine set_maxder_field_op_fem
    !---------------------------------------------------------------
    subroutine set_maxnen_fem(self)
        implicit none
        TYPE(FEM) :: self
        ! LOCAL
        INTEGER :: li_id

        self % oi_maxnen = self % opo_spaces (0) % oi_maxnen
        DO li_id = 1, self % oi_nspaces - 1
                IF (self % oi_maxnen < self % opo_spaces (li_id) % oi_maxnen) THEN
                        self % oi_maxnen = self % opo_spaces (li_id) % oi_maxnen
                END IF
        END DO

    end subroutine set_maxnen_fem
    !---------------------------------------------------------------
    subroutine get_maxnpatchs_fem(self, ai_val)
        implicit none
        TYPE(FEM) :: self
        INTEGER, intent(inout) :: ai_val

        ai_val = MAXVAL(self % opi_infoGrids(:, INFOGRIDS_NPATCHS))
    end subroutine get_maxnpatchs_fem
    !---------------------------------------------------------------
    subroutine get_maxnelts_fem(self, ai_val)
        implicit none
        TYPE(FEM) :: self
        INTEGER, intent(inout) :: ai_val
        ! LOCAL
        INTEGER :: li_grids
        INTEGER :: li_npatchs
        INTEGER :: li_val

        li_grids = 0
        li_npatchs  = self % opi_InfoGrids(li_grids, INFOGRIDS_NPATCHS)
        ai_val = MAXVAL(self % opi_InfoPatch(li_grids,0:li_npatchs-1, INFOPATCH_NEL))
        DO li_grids = 1, self % oi_nGrids - 1
           li_npatchs  = self % opi_InfoGrids(li_grids, INFOGRIDS_NPATCHS)
           li_val = MAXVAL(self % opi_InfoPatch(li_grids,0:li_npatchs-1, INFOPATCH_NEL))
           IF (li_val > ai_val) THEN
              ai_val = li_val
           END IF
        END DO

    end subroutine get_maxnelts_fem
    !---------------------------------------------------------------
    subroutine get_ngrids_fem(self, ai_nGrids)
        implicit none
        type(FEM) :: self
        INTEGER, INTENT(INOUT) :: ai_nGrids
        ! LOCAL

        ai_nGrids = self % oi_nGrids
    end subroutine get_ngrids_fem
    !---------------------------------------------------------------
    !> \todo to change : we must do  it for each grids
    subroutine get_maxnpoints_fem(self, ai_maxnpts)
        implicit none
        type(FEM) :: self
        INTEGER, INTENT(INOUT) :: ai_maxnpts
        ! LOCAL
        INTEGER :: li_grids
        INTEGER :: li_npatchs
        INTEGER :: li_maxnpts

        li_grids = 0
        li_npatchs  = self % opi_InfoGrids(li_grids, INFOGRIDS_NPATCHS)
        ai_maxnpts = MAXVAL(self % opi_InfoPatch(li_grids,0:li_npatchs-1, INFOPATCH_MAXNPTS))
        DO li_grids = 1, self % oi_nGrids - 1
           li_npatchs  = self % opi_InfoGrids(li_grids, INFOGRIDS_NPATCHS)
           li_maxnpts = MAXVAL(self % opi_InfoPatch(li_grids,0:li_npatchs-1, INFOPATCH_MAXNPTS))
           IF (li_maxnpts > ai_maxnpts) THEN
              ai_maxnpts = li_maxnpts
           END IF
        END DO

    end subroutine get_maxnpoints_fem
    !---------------------------------------------------------------
    subroutine get_operator_maxnderiv_fem(self, ai_nderiv)
        implicit none
        type(FEM) :: self
        INTEGER, INTENT(INOUT) :: ai_nderiv
        ! LOCAL
        INTEGER :: li_val

        if ( self % oi_nOperators == 0 ) then
            li_val = 1
        else
            li_val = MAXVAL(self % opi_InfoOperator(:, INFOOPERATOR_NDERIV))
        end if

        ai_nderiv = li_val

    end subroutine get_operator_maxnderiv_fem
    !---------------------------------------------------------------
    subroutine get_mapping_maxnderiv_fem(self, ai_nderiv)
        implicit none
        type(FEM) :: self
        INTEGER, INTENT(INOUT) :: ai_nderiv
        ! LOCAL
        INTEGER :: li_val

!        if ( self % oi_nmatrices == 0 ) then
!            li_val = 1
!        else
!            li_val = MAXVAL(self % opi_InfoMatrix(:, INFOMATRIX_NDERIV))
!        end if
!
!        ai_nderiv = li_val
        ai_nderiv = 2

    end subroutine get_mapping_maxnderiv_fem
    !---------------------------------------------------------------
    subroutine get_maxndof_fem(self, ai_maxndof)
        implicit none
        type(FEM) :: self
        INTEGER, INTENT(INOUT) :: ai_maxndof
        ! LOCAL
        INTEGER :: li_val

        if ( self % oi_nfields == 0 ) then
            li_val = 1
        else
            li_val = MAXVAL(self % opi_InfoField(:, INFOFIELD_NDOF))
        end if
        
        ai_maxndof = li_val

    end subroutine get_maxndof_fem
    !---------------------------------------------------------------
    subroutine reset_matrix (self, ai_id )
        implicit none
        type(FEM) :: self
        integer :: ai_id
        ! LOCAL
        integer :: ierr

        CALL SPM_MATRIXRESET(ai_id, ierr)

    end subroutine reset_matrix
!-------------------------------------------------------------------
    subroutine set_operator_matrices_toassembly ( self, ai_operator, api_values, ai_size )
    !> for eahc operator ,we specify the matrices to assembly
        implicit none
        TYPE(FEM) :: self        
        INTEGER, intent(in)  :: ai_operator
        INTEGER, intent(in)  :: ai_size
        INTEGER , dimension(ai_size), intent(in)  :: api_values
        ! LOCAL
        INTEGER :: li_i
        INTEGER :: li_matrix
        INTEGER :: li_n

!        print *, 'set_operator_matrices_toassembly: Begin'
!        print *, '---'
        self % opo_Op (ai_operator) % opi_matrices_toassembly (0) = ai_size
        self % opo_Op (ai_operator) % opi_matrices_toassembly (1:ai_size) = api_values(1:ai_size)
!        print *, '---'
!
!        print *, '---'
        li_n = self % opo_Op (ai_operator) % opi_matrices_toassembly (0)
        DO li_i = 1, li_n
!        print *, 'i=',li_i
        li_matrix = self % opo_Op (ai_operator) % opi_matrices_toassembly (li_i)
!        print *, ai_operator, li_matrix
!        print *, '---'
!        CALL reset_matrix (self, li_matrix )
!        print *, '---'
        END DO        
!        print *, 'set_operator_matrices_toassembly: End'
        
    end subroutine set_operator_matrices_toassembly    
    !---------------------------------------------------------------
    subroutine init_space_maxnen_fem(self)
        implicit none
        type(FEM) :: self
        ! LOCAL
        INTEGER :: li_id
        INTEGER :: li_maxnen
        INTEGER :: li_composed
#ifdef _TRACE
        CALL printlog("init_space_maxnen_fem : Start", ai_dtllevel = 1)
#endif

        DO li_id = 0, self % oi_nspaces - 1
            ! if it is a composed space, we do not need to defin the geometry
!            li_composed = self % opi_InfoSpace(li_id, INFOSPACE_COMPOSED)
!            IF (li_composed == 1) THEN
!                CYCLE
!            END IF

!            CALL print_geometry (self % opo_spaces(li_id) % oo_mapping % opo_geo(1))
            
!            CALL get_maxnen_geos(self % opo_spaces(li_id) % oo_mapping, li_maxnen)
!            self % opo_spaces(li_id) % oi_maxnen = li_maxnen

            self % opo_spaces(li_id) % oi_maxnen = self % opi_InfoSpace ( li_id, INFOSPACE_MAXNEN )
!            print *, 'space ', li_id, 'with maxnen = ', self % opo_spaces(li_id) % oi_maxnen

        END DO

        CALL set_maxnen_fem(self)

#ifdef _TRACE
        CALL printlog("init_space_maxnen_fem : End", ai_dtllevel = 1)
#endif

    end subroutine init_space_maxnen_fem
    !---------------------------------------------------------------
    subroutine init_mapping_maxnen_fem(self)
        implicit none
        type(FEM) :: self
        ! LOCAL
        INTEGER :: li_id
        INTEGER :: li_maxnen
#ifdef _TRACE
        CALL printlog("init_mapping_maxnen_fem : Start", ai_dtllevel = 1)
#endif

        if (self % oi_nmappings > 0) then
            DO li_id = 0, self % oi_nmappings - 1
                CALL get_maxnen_geos(self % opo_mappings(li_id) % oo_mapping, li_maxnen)
                self % opo_mappings(li_id) % oi_maxnen = li_maxnen
            END DO
        end if

#ifdef _TRACE
        CALL printlog("init_mapping_maxnen_fem : End", ai_dtllevel = 1)
#endif

    end subroutine init_mapping_maxnen_fem
    !---------------------------------------------------------------
    subroutine init_dim_Rd_fem(self)
        implicit none
        TYPE(FEM) :: self
        ! LOCAL
        integer :: li_grids
#ifdef _TRACE
        CALL printlog("init_dim_Rd_fem : Start", ai_dtllevel = 1)
#endif

        DO li_grids = 0, self % oi_nGrids - 1
            self % opi_dim(li_grids) = self % opi_infoGrids(li_grids, INFOGRIDS_DIM)
            self % opi_Rd(li_grids) = self % opi_infoGrids(li_grids, INFOGRIDS_RD)
        END DO
#ifdef _TRACE
        CALL printlog("init_dim_Rd_fem : End", ai_dtllevel = 1)
#endif

    end subroutine init_dim_Rd_fem
    !---------------------------------------------------------------
    subroutine init_maxp_spaces_fem(self)
        implicit none
        type(FEM) :: self
        ! LOCAL
        integer :: li_id
        integer :: li_geo
        integer :: li_maxp
#ifdef _TRACE
        CALL printlog("init_maxp_spaces_fem : Start", ai_dtllevel = 1)
#endif

        DO li_id = 0, self % oi_nspaces - 1

            li_maxp = 0

            DO li_geo = 1, self % opo_spaces(li_id) % oo_mapping % oi_ngeo

                IF (li_maxp < MAXVAL(self % opo_spaces(li_id) % oo_mapping % opo_geo(li_geo) % opi_P)) THEN
                    li_maxp = MAXVAL(self % opo_spaces(li_id) % oo_mapping % opo_geo(li_geo) % opi_P)
                END IF

            END DO

            self % opo_spaces(li_id) % oi_maxp = li_maxp

        END DO
#ifdef _TRACE
        CALL printlog("init_maxp_spaces_fem : End", ai_dtllevel = 1)
#endif

    end subroutine init_maxp_spaces_fem
    !---------------------------------------------------------------
    subroutine init_maxp_mappings_fem(self)
        implicit none
        type(FEM) :: self
        ! LOCAL
        integer :: li_id
        integer :: li_geo
        integer :: li_maxp
#ifdef _TRACE
        CALL printlog("init_maxp_mappings_fem : Start", ai_dtllevel = 1)
#endif

!        print*, 'self % oi_nmappings=', self % oi_nmappings
        DO li_id = 0, self % oi_nmappings - 1
!            print*, 'li_id=', li_id
            li_maxp = 0

            DO li_geo = 1, self % opo_mappings(li_id) % oo_mapping % oi_ngeo
!                print*, 'li_geo=', li_geo

                IF (li_maxp < MAXVAL(self % opo_mappings(li_id) % oo_mapping % opo_geo(li_geo) % opi_P)) THEN
                    li_maxp = MAXVAL(self % opo_mappings(li_id) % oo_mapping % opo_geo(li_geo) % opi_P)
                END IF

            END DO

            self % opo_mappings(li_id) % oi_maxp = li_maxp

        END DO
#ifdef _TRACE
        CALL printlog("init_maxp_mappings_fem : End", ai_dtllevel = 1)
#endif

    end subroutine init_maxp_mappings_fem
    !---------------------------------------------------------------
    subroutine create_fem_partone(self)
        implicit none
        TYPE(FEM) :: self
        ! LOCAL
#ifdef _TRACE
        CALL printlog("create_fem_partone : Start", ai_dtllevel = 1)
#endif

        ALLOCATE ( self % opi_dim(0:self % oi_nGrids - 1))
        ALLOCATE ( self % opi_Rd(0:self % oi_nGrids - 1))

        ALLOCATE ( self % opi_InfoSolver(0:self % oi_nsolvers - 1, NPARAM_INFOSOLVER))
        ALLOCATE ( self % opi_InfoMatrix(0:self % oi_nmatrices - 1, NPARAM_INFOMATRIX))
        ALLOCATE ( self % opi_InfoGraph(0:self % oi_ngraphs - 1, NPARAM_INFOGRAPH))
        ALLOCATE ( self % opi_InfoOperator(0:self % oi_nOperators - 1, NPARAM_INFOOPERATOR))
        ALLOCATE ( self % opi_InfoField(0:self % oi_nFields - 1, NPARAM_INFOFIELD))
        ALLOCATE ( self % opi_InfoSpace(0:self % oi_nspaces - 1, NPARAM_INFOSPACE))
        ALLOCATE ( self % opi_InfoMapping(0:self % oi_nmappings - 1, NPARAM_INFOMAPPING))
        ALLOCATE ( self % opi_InfoGrids(0:self % oi_nGrids - 1, NPARAM_INFOGRIDS))
        ALLOCATE ( self % opi_InfoPatch(0:self % oi_nGrids - 1, 0:self % oi_maxnpatchs - 1, NPARAM_INFOPATCH))
        ALLOCATE ( self % opi_InfoNorm(0:self % oi_nnorms - 1, NPARAM_INFONORM))

        ALLOCATE ( self % opo_mappings(0:self % oi_nmappings - 1))
        ALLOCATE ( self % opo_metrics(0:self % oi_nmetrics - 1))
        ALLOCATE ( self % opo_spaces(0:self % oi_nspaces - 1))
        ALLOCATE ( self % opo_Op(0:self % oi_nOperators - 1))

#ifdef _TRACE
        CALL printlog("create_fem_partone : End", ai_dtllevel = 1)
#endif

    end subroutine create_fem_partone
    !---------------------------------------------------------------
    subroutine create_fem_parttwo(self)
        implicit none
        TYPE(FEM) :: self
        ! LOCAL

#ifdef _TRACE
        CALL printlog("create_fem_parttwo : Start", ai_dtllevel = 1)
#endif

        CALL init_dim_Rd_fem(self)

        CALL allocate_spaces(self)
        
#ifdef _TRACE
        CALL printlog("create_fem_parttwo : End", ai_dtllevel = 1)
#endif

    end subroutine create_fem_parttwo
    !---------------------------------------------------------------
    subroutine create_fem_partthree(self)
        implicit none
        TYPE(FEM) :: self
        ! LOCAL
#ifdef _TRACE
        CALL printlog("create_fem_partthree : Start", ai_dtllevel = 1)
#endif

!        print *, 'allocate_all_grids'
        CALL allocate_all_grids(self)
!        print *, 'done.'

!        print *, 'allocate_all_metrics'
        CALL allocate_all_metrics(self)
!        print *, 'done.'

!        print *, 'allocate_fields'
        CALL allocate_fields(self)
!        print *, 'done.'
        
!        print *, 'allocate_graphs'
        CALL allocate_graphs(self)
!        print *, 'done.'

!        print *, 'allocate_operators'
        CALL allocate_operators(self)
!        print *, 'done.'

!        print *, 'allocate_matrices'
        CALL allocate_matrices(self)
!        print *, 'done.'

!        print *, 'add_info_spaces'
        CALL add_info_spaces(self)
!        print *, 'done.'

!        print *, 'allocate_norms'
        CALL allocate_norms(self)
!        print *, 'done.'

!        print *, 'allocate_solvers'
        CALL allocate_solvers(self)
!        print *, 'done.'

!        STOP

        self % ol_initialized = .TRUE.
        
#ifdef _TRACE
        CALL printlog("create_fem_partthree : End", ai_dtllevel = 1)
#endif

    end subroutine create_fem_partthree
    !---------------------------------------------------------------
    subroutine deallocate_geometries(self)
        implicit none
        TYPE(FEM) :: self
        ! LOCAL
        INTEGER :: li_id
        INTEGER :: li_composed


        DO li_id = 0, self % oi_nspaces - 1
            ! if it is a composed space, we do not need to defin the geometry
!            li_composed = self % opi_InfoSpace(li_id, INFOSPACE_COMPOSED)
!            IF (li_composed == 1) THEN
!                CYCLE
!            END IF
            CALL free_geometries(self % opo_spaces(li_id) % oo_mapping)
        END DO

    end subroutine deallocate_geometries
    !---------------------------------------------------------------
    subroutine deallocate_connectivities(self)
        implicit none
        TYPE(FEM) :: self
        ! LOCAL
        INTEGER :: li_id
        INTEGER :: li_composed

        DO li_id = 0, self % oi_nspaces - 1
            ! if it is a composed space, we do not need to defin the geometry
!            li_composed = self % opi_InfoSpace(li_id, INFOSPACE_COMPOSED)
!            IF (li_composed == 1) THEN
!                CYCLE
!            END IF
            CALL free_connectivity(self % opo_spaces(li_id) % oo_con)
        END DO
        
    end subroutine deallocate_connectivities
    !---------------------------------------------------------------
    subroutine deallocate_mapping_geometries(self)
        implicit none
        TYPE(FEM) :: self
        ! LOCAL
        INTEGER :: li_id
        INTEGER :: li_composed

        DO li_id = 0, self % oi_nmappings - 1
            CALL free_geometries(self % opo_mappings(li_id) % oo_mapping)
        END DO

    end subroutine deallocate_mapping_geometries
    !---------------------------------------------------------------
    subroutine free_fem(self)
        implicit none
        TYPE(FEM) :: self
        ! LOCAL
        IF (self % ol_initialized) THEN
                CALL deallocate_graphs(self)
                CALL deallocate_operators(self)
                CALL deallocate_matrices(self)

                CALL deallocate_solvers(self)

                CALL deallocate_norms(self)
                CALL deallocate_fields(self)

                CALL deallocate_geometries(self)
                CALL deallocate_connectivities(self)
                CALL deallocate_mapping_geometries(self)

                CALL deallocate_all_grids(self)

                CALL deallocate_all_metrics(self)

                DEALLOCATE ( self % opi_InfoSolver)
                DEALLOCATE ( self % opi_InfoMatrix)
                DEALLOCATE ( self % opi_InfoGraph)
                DEALLOCATE ( self % opi_InfoOperator)
                DEALLOCATE ( self % opi_InfoField)
                DEALLOCATE ( self % opi_InfoSpace)
                DEALLOCATE ( self % opi_InfoMapping)
                DEALLOCATE ( self % opi_InfoPatch)
                DEALLOCATE ( self % opi_InfoGrids)
                DEALLOCATE ( self % opi_InfoNorm)

                DEALLOCATE ( self % opo_mappings)
                DEALLOCATE ( self % opo_metrics)
                DEALLOCATE ( self % opo_spaces)
                DEALLOCATE ( self % opo_Op)

                DEALLOCATE ( self % opi_dim)
                DEALLOCATE ( self % opi_Rd)
                
        END IF

    end subroutine free_fem
!----------------------------------------------------------------------------------------------
    subroutine getarrayparamcsrmatrix_fem ( self, ai_ref, apr_a, api_ia, api_ja, ai_nR, ai_nel )
        implicit none
        TYPE(FEM) :: self
        integer, intent(in)  :: ai_ref  !THE MATRIX FROM WICH WE WANT TO ADD CONTRIBUTION
        integer, intent(in)  :: ai_nR   !NUMBER OF ROWS
        integer*8, intent(in)  :: ai_nel  !NUMBER OF NON ZERO ELTS
        real*8 , dimension(ai_nel), intent(out)  :: apr_a
        integer, dimension(ai_nR+1), intent(out)  :: api_ia
        integer*8, dimension(ai_nel), intent(out)  :: api_ja
        ! LOCAL VARIABLES
        integer  :: ierr

#ifdef _TRACE
        CALL printlog("getarrayparamcsrmatrix_fem : Begin", ai_dtllevel = 1)
#endif

        CALL SPM_GetINDPTR  (ai_ref, ai_nR+1, api_ia, ierr)
        CALL SPM_GetINDICES (ai_ref, ai_nel , api_ja, ierr)
        CALL SPM_GetDATA    (ai_ref, ai_nel , apr_a , ierr)

        ! for 0 based indices
        api_ia = api_ia - 1
        api_ja = api_ja - 1

#ifdef _TRACE
        CALL printlog("getarrayparamcsrmatrix_fem : End", ai_dtllevel = 1)
#endif

    end subroutine getarrayparamcsrmatrix_fem
!----------------------------------------------------------------------------------------------
    subroutine createcsrmatrix_fem ( self, ai_ref, api_ia, api_ja, ai_nR, ai_nC, ai_nel )
        implicit none
        TYPE(FEM) :: self
        integer, intent(in)  :: ai_ref  !THE MATRIX FROM WICH WE WANT TO ADD CONTRIBUTION
        integer, intent(in)  :: ai_nR   !NUMBER OF ROWS
        integer, intent(in)  :: ai_nC   !NUMBER OF COLUMNS
        integer*8, intent(in)  :: ai_nel  !NUMBER OF NON ZERO ELTS
        integer, dimension(ai_nR+1), intent(in)  :: api_ia
        integer*8, dimension(ai_nel), intent(in)  :: api_ja
        ! LOCAL VARIABLES
        integer  :: ierr
        integer  :: root 
        integer, dimension(ai_nR+1) :: lpi_ia
        integer*8, dimension(ai_nel)  :: lpi_ja

#ifdef _TRACE
        CALL printlog("createcsrmatrix_fem : Begin", ai_dtllevel = 1)
#endif

        root = -1
        ! for 0 based indices
        lpi_ia = api_ia + 1
        lpi_ja = api_ja + 1

        CALL SPM_GraphGlobalCSR(ai_ref, ai_nC, lpi_ja, lpi_ia, root, ierr)

#ifdef _TRACE
        CALL printlog("createcsrmatrix_fem : End", ai_dtllevel = 1)
#endif

    end subroutine createcsrmatrix_fem
!----------------------------------------------------------------------------------------------
    subroutine setcsrmatrix_fem ( self, ai_ref, apr_a, ai_nel )
        implicit none
        TYPE(FEM) :: self
        integer, intent(in)  :: ai_ref  !THE MATRIX FROM WICH WE WANT TO ADD CONTRIBUTION
        integer*8, intent(in)  :: ai_nel  !NUMBER OF NON ZERO ELTS
        real*8 , dimension(ai_nel), intent(in)  :: apr_a
        ! LOCAL VARIABLES
        integer  :: ierr
        integer  :: root 

#ifdef _TRACE
        CALL printlog("setcsrmatrix_fem : Begin", ai_dtllevel = 1)
#endif

        root = -1
        CALL SPM_SetDATA(ai_ref, ai_nel, apr_a, ierr)

#ifdef _TRACE
        CALL printlog("setcsrmatrix_fem : End", ai_dtllevel = 1)
#endif

    end subroutine setcsrmatrix_fem    
end module fem_module
!**************************************************
