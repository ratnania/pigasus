!     
! File:   assembly.F90
! Author: root
!
! Created on January 2, 2012, 2:25 PM
!

module assembly_module
    use tracelog_module
    use assembly_def
    use grids_def
    use geometries_def
    use bbox_def
    use bbox_module
    use connectivities_def
    use fem_def
    use fem_module
    use operators_module
    use projectors_module
    use norms_module
    use fields_operators_module
    use assembly_tools
    implicit none

#ifdef _DEBUG
    integer, parameter, private  :: mi_dtllevel_base = 0
#else
    integer, parameter, private  :: mi_dtllevel_base = 2
#endif

contains

    !----------------------------------------------------------------------------------------------
    subroutine Assembly_initialize(self,ao_FEM)
        implicit none
        TYPE(ASSEMBLY) :: self
        TYPE(FEM) :: ao_FEM
        ! LOCAL
        INTEGER :: li_id
        INTEGER :: li_npts
        INTEGER :: li_dirmaxnpts
        INTEGER :: li_maxp
        INTEGER :: li_max_nderiv
        INTEGER :: li_space
        INTEGER :: li_mapping
        INTEGER :: li_space_maxnen
        INTEGER :: li_mapping_maxnen
        INTEGER :: li_fieldh_maxnderiv
        INTEGER :: li_space_maxnnz
        LOGICAL :: ll_tensor

        CALL printlog("Assembly_initialize : Start", ai_dtllevel = mi_dtllevel_base  + 1)

    !***********************************************
    !           setting parameters, needed
    !           for memory allocations
    !***********************************************
        ! computing the maximum of ndof
        CALL get_maxndof_fem( ao_FEM, self % oi_maxndof )
        ! computing the max npts per element
        CALL get_maxnpoints_fem( ao_FEM, self % oi_maxnpts )
        ! computing the max nbasis fct per element
        CALL init_space_maxnen_fem   ( ao_FEM )
        CALL init_mapping_maxnen_fem ( ao_FEM )

        ! compute for each space, the maximum polynomial degree
        CALL init_maxp_spaces_fem ( ao_FEM )
        ! compute for each mapping, the maximum polynomial degree
        CALL init_maxp_mappings_fem ( ao_FEM )

        !> compute the highet n-derivatives
        CALL set_operators_nderiv( self, ao_FEM )
        CALL get_mapping_maxnderiv_fem ( ao_FEM, self % oi_mapping_nderiv )
        CALL get_operator_maxnderiv_fem ( ao_FEM, li_max_nderiv )
        CALL get_maxnelts_fem(ao_FEM, ao_FEM % oi_maxnelts)

        ! ********************************************************
        !                   SPACES TREATMENT
        ! COMPUTING THE MAXIMUM OF ORDER DERIVATIVES OF ALL BASIS FUNCTIONS
        ! ********************************************************
        do li_space = 0, ao_FEM % oi_nspaces-1

            !> if we don't use an exterior mapping, then we must use the transformation defined
            !> by the current geometry
!            if ( ao_FEM % opi_InfoSpace ( li_space , INFOSPACE_EXTMAPPING ) == 0 ) then
                li_max_nderiv = max (li_max_nderiv,1)
!            end if

        end do
        self % oi_matrix_nderiv = li_max_nderiv
        li_fieldh_maxnderiv = MAXVAL (ao_FEM % opi_dim(:))**2  ! if we use SECOND_ORDER_FIELD

        ! ********************************************************
        li_space_maxnen = 0
        do li_space = 0, ao_FEM % oi_nspaces-1
         
            if (li_space_maxnen < ao_FEM % opo_spaces(li_space) % oi_maxnen ) then
                li_space_maxnen = ao_FEM % opo_spaces(li_space) % oi_maxnen
            end if
        end do
        ! computing mapping_maxnen
        li_mapping_maxnen = 0
        do li_mapping = 0, ao_FEM % oi_nmappings-1
            if (li_mapping_maxnen < ao_FEM % opo_mappings(li_mapping) % oi_maxnen ) then
                li_mapping_maxnen = ao_FEM % opo_mappings(li_mapping) % oi_maxnen
            end if
        end do

       
    !***********************************************

!    print *, ao_FEM % oi_maxnpatchs
!    print *, ao_FEM % oi_maxnelts
!    print *, ao_FEM % oi_nMatrices
!    print *, ao_FEM % oi_nFields
!    print *, ao_FEM % oi_nnorms
!    print *, ao_FEM % oi_nmappings
!    print *, ao_FEM % oi_nspaces
!    print *, ao_FEM % oi_ngrids
!    print *, li_space_maxnen
!    print *, li_fieldh_maxnderiv
!    print *, self % oi_maxndof
!    print *, self % oi_maxnpts
!    print *, ao_FEM % 
    !***********************************************
    !           MEMORY ALLOCATIONS
    !***********************************************
        !> \todo WE HAVE TO ALLOCATE ONLY IF THE CORRESPONDING FLAG LL_ASSEMBLY IS ACTIVATED

        ! in self % opi_matrices_toassembly(0) we store the number of matrices to assembly
        ALLOCATE (self % opi_patchs_toassembly(0:ao_FEM % oi_maxnpatchs))
        ALLOCATE (self % opi_patchs_toassembly_tmp(0:ao_FEM % oi_maxnpatchs))

        ALLOCATE (self % opi_elts_toassembly(0:ao_FEM % oi_maxnelts))
        ALLOCATE (self % opi_elts_toassembly_tmp(0:ao_FEM % oi_maxnelts))

        ALLOCATE (self % opi_matrices_toassembly(0:ao_FEM % oi_nMatrices))
        ALLOCATE (self % opi_matrices_toassembly_tmp(0:ao_FEM % oi_nMatrices))

        ALLOCATE (self % opi_operators_toassembly(0:ao_FEM % oi_nOperators))
        ALLOCATE (self % opi_operators_toassembly_tmp(0:ao_FEM % oi_nOperators))

        ALLOCATE (self % opi_fields_toassembly(0:ao_FEM % oi_nFields))
        ALLOCATE (self % opi_fields_toassembly_tmp(0:ao_FEM % oi_nFields))

        ALLOCATE (self % opi_norms_toassembly(0:ao_FEM % oi_nnorms))
        ALLOCATE (self % opi_norms_toassembly_tmp(0:ao_FEM % oi_nnorms))

        ALLOCATE (self % opi_spaces_toassembly(0:ao_FEM % oi_nspaces))
        ALLOCATE (self % opi_spaces_toassembly_tmp(0:ao_FEM % oi_nspaces))

        ALLOCATE (self % opo_bbox_mp(0:ao_FEM % oi_nmappings-1))
        ALLOCATE (self % opo_info_mp(0:ao_FEM % oi_nmappings-1))

        ALLOCATE (self % opo_bbox_sp(0:ao_FEM % oi_nspaces-1))
        ALLOCATE (self % opo_info_sp(0:ao_FEM % oi_nspaces-1))

        ALLOCATE (self % opo_info_gr(0:ao_FEM % oi_ngrids-1))

        ALLOCATE (self % opo_pBasis(0:ao_FEM % oi_nspaces-1))

        CALL create_spaces_bbox     (self, ao_FEM)
        CALL create_mappings_bbox   (self, ao_FEM)

        CALL create_spaces_info     (self, ao_FEM)
        CALL create_mappings_info   (self, ao_FEM)

        CALL create_physical_basis  (self, ao_FEM)

        CALL create_grids_info   (self, ao_FEM)

        CALL create_spaces_stored_data( self, ao_FEM )
        CALL create_mappings_stored_data( self, ao_FEM )

        !***********
        ! level 1
        !***********
        ALLOCATE (self % opr_Matrix_elt(0:ao_FEM % oi_nOperators-1, li_space_maxnen, li_space_maxnen))
        ALLOCATE (self % opr_Projection_elt(0:ao_FEM % oi_nFields-1, li_space_maxnen))
        ALLOCATE (self % opr_Fieldh_elt(0:ao_FEM % oi_nFields-1, self % oi_maxndof, li_fieldh_maxnderiv, self % oi_maxnpts))
        ALLOCATE (self % opr_Norm_elt(0:ao_FEM % oi_nNorms-1))
        !***********

        !***********
        ! level 2
        !***********
        li_space_maxnnz = 0
!        do li_space = 0, ao_FEM % oi_nspaces-1
!            IF (ao_FEM % opo_spaces(li_space) % oi_tensorlevel == 2 ) THEN 
!                IF (li_space_maxnnz < assl_get_arbnnz(self % opo_bbox_sp ( li_space ) % oo_assl) ) THEN
!                    li_space_maxnnz = assl_get_arbnnz(self % opo_bbox_sp ( li_space ) % oo_assl)
!                END IF
!            END IF
!        end do
        ALLOCATE (self % opr_Matrix_eltLine(li_space_maxnnz))
        !***********
    !***********************************************

    !***********************************************
    !               USING MURGE
    !***********************************************
        CALL Assembly_initialize_MURGE ( self, ao_FEM)
    !***********************************************

        self % ol_initialized = .TRUE.

        CALL printlog("Assembly_initialize : End", ai_dtllevel = mi_dtllevel_base  + 1)

    end subroutine Assembly_initialize
    !----------------------------------------------------------------------------------------------
    subroutine Assembly_finalize(self,ao_FEM)
        implicit none
        TYPE(ASSEMBLY) :: self
        TYPE(FEM) :: ao_FEM
        ! LOCAL
        INTEGER :: li_id
        INTEGER :: ierr
#ifdef _MURGE
        INTEGER(KIND = MURGE_INTS_KIND) :: ierr
        INTEGER(KIND = MURGE_INTS_KIND) :: id
#endif

        CALL printlog("Assembly_finalize : Start", ai_dtllevel = mi_dtllevel_base  + 1)

        IF (self % ol_initialized) THEN

                CALL free_mappings_stored_data( self, ao_FEM )
                CALL free_mappings_bbox (self, ao_FEM)
                CALL free_mappings_info (self, ao_FEM)
                CALL free_spaces_info (self, ao_FEM)
                CALL free_spaces_bbox (self, ao_FEM)
                CALL free_physical_basis (self, ao_FEM)
                CALL free_grids_info (self, ao_FEM)
                CALL free_spaces_stored_data( self, ao_FEM )

                DEALLOCATE (self % opr_Matrix_elt)
                DEALLOCATE (self % opr_Projection_elt)
                DEALLOCATE (self % opr_Fieldh_elt)
                DEALLOCATE (self % opr_Norm_elt)

                DEALLOCATE (self % opr_Matrix_eltLine)

                DEALLOCATE (self % opi_patchs_toassembly)
                DEALLOCATE (self % opi_patchs_toassembly_tmp)
                DEALLOCATE (self % opi_elts_toassembly)
                DEALLOCATE (self % opi_elts_toassembly_tmp)
                DEALLOCATE (self % opi_operators_toassembly)
                DEALLOCATE (self % opi_operators_toassembly_tmp)
                DEALLOCATE (self % opi_matrices_toassembly)
                DEALLOCATE (self % opi_matrices_toassembly_tmp)
                DEALLOCATE (self % opi_fields_toassembly)
                DEALLOCATE (self % opi_fields_toassembly_tmp)
                DEALLOCATE (self % opi_norms_toassembly)
                DEALLOCATE (self % opi_norms_toassembly_tmp)
                DEALLOCATE (self % opi_spaces_toassembly)
                DEALLOCATE (self % opi_spaces_toassembly_tmp)
                DEALLOCATE (self % opo_bbox_mp)
                DEALLOCATE (self % opo_bbox_sp)
                DEALLOCATE (self % opo_info_mp)
                DEALLOCATE (self % opo_info_sp)
                DEALLOCATE (self % opo_info_gr)
                DEALLOCATE (self % opo_pBasis)

#ifdef _MURGE
                ! I'm Free
                DO id = 0, ao_FEM % oi_nMatrices - 1
                    CALL MURGE_CLEAN(id, ierr)
                    if (ierr /= MURGE_SUCCESS) call abort()
                END DO
                CALL MURGE_FINALIZE(ierr)
                CALL MPI_FINALIZE(ierr)
#endif

        END IF

        CALL printlog("Assembly_finalize : End", ai_dtllevel = mi_dtllevel_base  + 1)

    end subroutine Assembly_finalize
    !----------------------------------------------------------------------------------------------
    subroutine Assembly_Matrix(self, ao_FEM, ai_grids_id)
        !> Assembling process over the grid of id ai_grids_id
        !> it is supposed that all field, matrices, norms are defined on this grid
        implicit none
        TYPE(ASSEMBLY) :: self
        TYPE(FEM) :: ao_FEM
        INTEGER :: ai_grids_id
        ! LOCAL VARIABLES
        integer :: li_elt
        integer :: li_idof
        integer :: li_err
        integer :: li_flag
        integer :: li_index
        integer :: li_quotient_elt
        integer :: li_field
        integer :: li_field1
        integer :: li_field2
        integer :: li_fieldprime
        integer :: li_matrix
        integer :: li_source
        integer :: li_ndof
        integer :: li_ndofprime
        integer :: li_type_BB
        integer :: li_type_vect
        integer :: li_b1
        integer :: li_b2
        integer :: li_ii
        integer :: li_id
        integer :: li_map
        integer :: li_space
        INTEGER :: li_ref_space
        INTEGER :: li_nspaces
        INTEGER :: li_ref_elt
        integer :: li_nel
        INTEGER :: li_npts
        INTEGER :: li_nmatrices
        INTEGER :: li_noperators
        INTEGER :: li_operator
        INTEGER :: li_nfields
        INTEGER :: li_ref
        INTEGER :: li_ref_patch
        INTEGER :: li_npatchs
        INTEGER :: li_nnorms
        INTEGER :: li_norm
        INTEGER :: li_position_sp
        INTEGER :: li_single_sp
        INTEGER :: li_single_grids_id
        INTEGER :: li_realelt
        INTEGER :: li_nfields_operators
        INTEGER :: li_nfields_l2projections
        real(wp) :: lr_u
        logical :: ll_computeU
        logical :: ll_valuesOnGrid
        type(GRID_DATA), pointer :: lp_grid
        type(GRID_DATA), pointer :: lp_single_grid
        type(BBOX), pointer :: lp_bbox
        type(MAPPING), pointer :: lp_mapping
        type(SPACE), pointer :: lp_space
        INTEGER, DIMENSION(:), POINTER :: lpi_fields_operators
        INTEGER, DIMENSION(:), POINTER :: lpi_fields_l2projections
        INTEGER*8 :: ncount
        INTEGER :: ierr
#ifdef _MURGE
        INTEGER(KIND = MURGE_INTS_KIND) :: ierr
        INTEGER(KIND = MURGE_INTS_KIND) :: id
#endif
        call printlog("Assembly_Matrix : Begin", ai_dtllevel = mi_dtllevel_base  + 1)
!        print *, "============== Assembly_Matrix : Begin ================"

        CALL set_spaces_toassembly(self,ao_FEM, ai_grids_id)
        li_nspaces = self % opi_spaces_toassembly(0)
!        print *, 'li_nspaces=', li_nspaces

        CALL set_patchs_toassembly(self,ao_FEM, ai_grids_id)
        li_npatchs = self % opi_patchs_toassembly(0)
!        print *, 'li_npatchs=', li_npatchs

        CALL set_operators_toassembly(self,ao_FEM)
        li_noperators = self % opi_operators_toassembly(0)
!        print *, 'li_noperators=', li_noperators
!        print *, "operators to assembly :", self % opi_operators_toassembly(1:li_noperators)

        CALL set_matrices_toassembly(self,ao_FEM)
        li_nmatrices = self % opi_matrices_toassembly(0)
!        print *, 'li_nmatrices=', li_nmatrices
!        print *, 'matrices to assembly: ', self % opi_matrices_toassembly(1:li_nmatrices)

        CALL set_fields_toassembly(self,ao_FEM)
        li_nfields = self % opi_fields_toassembly(0)
!        print *, 'li_nfields=', li_nfields

        CALL set_norms_toassembly(self,ao_FEM)
        li_nnorms = self % opi_norms_toassembly(0)
!        print *, 'li_nfields=', li_nfields

        ALLOCATE(lpi_fields_operators(0:li_nfields))
        ALLOCATE(lpi_fields_l2projections(0:li_nfields))

        CALL select_fields(self, ao_FEM, FIELD_OPERATOR, lpi_fields_operators)
        CALL select_fields(self, ao_FEM, PROJECTION_L2, lpi_fields_l2projections)

        CALL set_flag_points_basis_assembly(self, ao_FEM)

        DO li_ref = 1, li_nmatrices
            li_matrix = self % opi_matrices_toassembly(li_ref)
            CALL SPM_MATRIXRESET(li_matrix, ierr)
            ncount = 0 

            CALL SPM_ASSEMBLYBEGIN(li_matrix, ncount, SPM_ASSEMBLY_ADD, SPM_ASSEMBLY_ADD, &
            SPM_ASSEMBLY_FOOL, SPM_BOOLEAN_FALSE, ierr)
        END DO
        
        ! USED TO PRINT THE EVOLUTION OF THE COMPUTATION
        li_quotient_elt = 0
        li_err = 1

!        li_npatchs = ao_FEM % opi_infoGrids (ai_grids_id, INFOGRIDS_NPATCHS)

        li_index = 0
        ! LOOP THROUG PATCHS
!        DO li_id = 0, li_npatchs-1
        li_npatchs = self % opi_patchs_toassembly(0)

#ifdef _DEBUG
        call concatmsg("current GRID is = : ", ai_dtllevel = mi_dtllevel_base + 1)
        call concatmsg(ai_grids_id, ai_dtllevel = mi_dtllevel_base + 1)
        call printmsg(ai_dtllevel = mi_dtllevel_base + 1)
        call concatmsg("number of patchs is = : ", ai_dtllevel = mi_dtllevel_base + 1)
        call concatmsg(li_npatchs, ai_dtllevel = mi_dtllevel_base + 1)
        call printmsg(ai_dtllevel = mi_dtllevel_base + 1)
#endif
!        print *, '-- all nel            : ', ao_FEM % opi_InfoPatch ( :, 0, INFOPATCH_NEL)
!        print *, '-- Current Grid       : ', ai_grids_id
!        print *, '-- npatchs            : ', li_npatchs
        do li_ref_patch = 1, li_npatchs

            li_id = self % opi_patchs_toassembly(li_ref_patch)
!            print *, '---- Current Patch    : ', li_id

            lp_grid => ao_FEM % opo_grids(ai_grids_id) % opo_grid ( li_id )

            ! LOOP THROUG ELEMENTS
            li_nel = self % opi_elts_toassembly(0)
            IF ( (li_npatchs > 1) .OR. (li_nel == -1) ) THEN
                li_nel = ao_FEM % opi_InfoPatch ( ai_grids_id, li_id, INFOPATCH_NEL)
                DO li_elt = 1, li_nel
                    self % opi_elts_toassembly(li_elt) = li_elt
                END DO
                self % opi_elts_toassembly(0) = li_nel
!                print *, "opi_elts_toassembly = ", self % opi_elts_toassembly
            END IF
!            print *, 'li_nel = ', li_nel

#ifdef _DEBUG
            call concatmsg("number of elements is : ", ai_dtllevel = 0)
            call concatmsg(li_nel, ai_dtllevel = 0)
            call printmsg(ai_dtllevel = 0)
#endif

            do li_ref_elt = 1, li_nel

                li_elt = self % opi_elts_toassembly(li_ref_elt)
!                print *, '------ Current elt    : ', li_elt

#ifdef _DEBUG
        call concatmsg("current Element is : ", ai_dtllevel = mi_dtllevel_base + 1)
        call concatmsg(li_elt, ai_dtllevel = mi_dtllevel_base + 1)
        call printmsg(ai_dtllevel = mi_dtllevel_base + 1)
#endif

                !> \todo WE HAVE TO INITIALIZE ONLY IF THE CORRESPONDING FLAG LL_ASSEMBLY IS ACTIVATED
                self % opr_Matrix_elt(:,:,:) = 0.0_wp
                self % opr_Projection_elt(:,:) = 0.0_wp
                self % opr_fieldh_elt = 0.0_wp
                
                self % opr_Matrix_eltLine = 0.0_wp

                li_npts = lp_grid % opo_elts ( li_elt ) % oi_npts

                DO li_ref_space = 1, li_nspaces 
                    li_space = self % opi_spaces_toassembly(li_ref_space)

                    ! ********************************************************
                    !                   METRIC TREATMENT
                    ! ********************************************************
                    CALL update_info_metric( self, ao_FEM, li_space, li_id, li_elt, li_npts )
                    ! ********************************************************

                    ! ********************************************************
                    !                   SPACES TREATMENT
                    ! COMPUTING ALL DERIVATIVES OF ALL BASIS FUNCTIONS
                    ! ********************************************************
                    CALL update_info_space( self, ao_FEM, li_space, li_id, li_elt, li_npts )
                    ! ********************************************************

                    ! ********************************************************
                    !               ASSEMBLING PHYSICAL BASIS
                    ! ********************************************************
                    CALL assembly_physical_basis( self, ao_FEM, li_space, li_id, li_elt )
                    ! ********************************************************

                END DO

            ! ********************************************************
            !   ASSEMBLING LOCAL FIELD OPERATOR
            ! ********************************************************
                li_nfields_operators = lpi_fields_operators(0)
                do li_ref = 1, li_nfields_operators

                    li_field = lpi_fields_operators(li_ref)

                    li_space = ao_FEM % opi_InfoField  (li_field , INFOFIELD_SPACE)

                    CALL build_field_operators_local(self, ao_FEM, ai_grids_id, li_id, li_elt, li_field)
                    CALL Assembly_Field_Operator_from_elt(self, ao_FEM, li_id+1, li_elt, li_field)

!#ifdef _DEBUG
!                call concatmsg("self % opr_fieldh_elt= : ", ai_dtllevel = 0)
!                call concatmsg(self % opr_fieldh_elt, ai_dtllevel = 0)
!                call printmsg(ai_dtllevel = 0)

                end do
            ! ********************************************************

            ! ********************************************************
            !   ASSEMBLING LOCAL L2-PROJECTION
            ! ********************************************************
                li_nfields_l2projections = lpi_fields_l2projections(0)
                do li_ref = 1, li_nfields_l2projections

                    li_field = lpi_fields_l2projections(li_ref)

                    li_space = ao_FEM % opi_InfoField  (li_field , INFOFIELD_SPACE)

                    CALL build_field_projectors_Local(self, ao_FEM, ai_grids_id, li_id, li_elt, li_field)
                    CALL Assembly_Projections_from_elt(self, ao_FEM, li_id+1, li_elt, li_field)

                END DO
            ! ********************************************************

            ! ********************************************************
            !   ASSEMBLING LOCAL OPERATORS 
            ! ********************************************************
                do li_ref = 1, li_noperators

                    li_operator = self % opi_operators_toassembly(li_ref)

                    li_space = ao_FEM % opi_Infooperator (li_operator, INFOOPERATOR_SPACE_1)

                    CALL build_matrix_operators_Local(self, ao_FEM, ai_grids_id, li_id, li_elt, li_operator)
                    CALL Assembly_Matrices_from_elt   (self, ao_FEM, ai_grids_id, li_id, li_elt, li_operator)

                end do
            ! ********************************************************

#ifdef _DEBUG
            call concatmsg("self % opr_Matrix_elt= : ", ai_dtllevel = mi_dtllevel_base + 2)
            call concatmsg(self % opr_Matrix_elt, ai_dtllevel = mi_dtllevel_base + 2)
            call printmsg(ai_dtllevel = mi_dtllevel_base + 2)
#endif

            ! ********************************************************
            !   ASSEMBLING LOCAL NORMS
            ! ********************************************************
                do li_ref = 1, li_nnorms

                    li_norm = self % opi_norms_toassembly(li_ref)

                    li_field = ao_FEM % opi_InfoNorm   (li_norm  , INFONORM_FIELD)
                    li_space = ao_FEM % opi_InfoField  (li_field , INFOFIELD_SPACE)

                    CALL build_norm_operators_Local(self, ao_FEM, ai_grids_id, li_id, li_elt, li_norm)
                    CALL Assembly_Norms_from_elt(self, ao_FEM, li_id, li_elt, li_norm)

                end do

            end do

        END DO

        DEALLOCATE(lpi_fields_operators)
        DEALLOCATE(lpi_fields_l2projections)

    DO li_ref = 1, li_nmatrices
        li_matrix = self % opi_matrices_toassembly(li_ref)
        CALL SPM_ASSEMBLYEND(li_matrix, ierr)
    END DO

#ifdef _MURGE
    DO id = 0, ao_FEM % oi_nMatrices - 1
        CALL MURGE_ASSEMBLYEND(id, ierr)
    END DO
#endif
    
!print *, "============== Assembly_Matrix : End    ================"

        call printlog("Assembly_Matrix : End", ai_dtllevel = mi_dtllevel_base + 1)

        call printcputime()

    end subroutine Assembly_Matrix
!----------------------------------------------------------------------------------------------
!    subroutine Assembly_Composed_Matrix(self, ao_FEM, ai_matrix, api_id)
!        implicit none
!        type(ASSEMBLY) :: self
!        type(FEM) :: ao_FEM
!        integer :: ai_matrix
!        integer, dimension(:,:) :: api_id
!        ! LOCAL VARIABLES
!        integer :: li_dim1
!        integer :: li_dim2
!        integer, dimension(:,:), pointer :: lpi_id
!
!        call printlog("Assembly_Composed_Matrix : Start", ai_dtllevel = mi_dtllevel_base + 1)
!
!        li_dim1 = ao_FEM % opi_InfoComposedMatrix ( ai_matrix , INFOCOMPOSEDMATRIX_DIM1 )
!        li_dim2 = ao_FEM % opi_InfoComposedMatrix ( ai_matrix , INFOCOMPOSEDMATRIX_DIM2 )
!
!        ALLOCATE(lpi_id(li_dim1,li_dim2))
!
!        lpi_id (1,1) = ao_FEM % opi_InfoComposedMatrix ( ai_matrix , INFOCOMPOSEDMATRIX_MATID1 )
!        lpi_id (1,2) = ao_FEM % opi_InfoComposedMatrix ( ai_matrix , INFOCOMPOSEDMATRIX_MATID2 )
!        lpi_id (2,1) = ao_FEM % opi_InfoComposedMatrix ( ai_matrix , INFOCOMPOSEDMATRIX_MATID3 )
!        lpi_id (2,2) = ao_FEM % opi_InfoComposedMatrix ( ai_matrix , INFOCOMPOSEDMATRIX_MATID4 )
!
!        CALL copy_csr(ao_FEM % opo_CCR_Matrix ( ai_matrix ) &
!            , ao_FEM % opo_CCR_Matrix ( lpi_id (1,1) )      &
!            , ao_FEM % opo_CCR_Matrix ( lpi_id (1,2) )      &
!            , ao_FEM % opo_CCR_Matrix ( lpi_id (2,1) )      &
!            , ao_FEM % opo_CCR_Matrix ( lpi_id (2,2) )      &
!            , api_id )
!
!!        print*, 'Mat = ', ao_FEM % opo_CCR_Matrix ( ai_matrix ) % opr_a (:)
!
!        DEALLOCATE(lpi_id)
!  
!        call printlog("Assembly_Composed_Matrix : End", ai_dtllevel = mi_dtllevel_base + 1)
!
!    end subroutine Assembly_Composed_Matrix
end module assembly_module
