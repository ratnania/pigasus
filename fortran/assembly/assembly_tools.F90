!     
! File:   assembly_tools.F90
! Author: root
!
! Created on January 30, 2012, 4:56 PM
!

module assembly_tools
    USE tracelog_module
    USE bbox_def
    USE bbox_module

    USE grids_def
    USE grids, ONLY : assembly_logical_basis, assembly_points_elements, create_stored_data, free_stored_data
    USE fem_def
    USE OPERATOR_MODULE
!    use fem_module
    USE assembly_def
    IMPLICIT NONE

!#ifdef _MURGE
!    INCLUDE "mpif.h"
!    INCLUDE "murge.inc"
!#endif

#ifdef _DEBUG
    INTEGER, PARAMETER, PRIVATE  :: mi_dtllevel_base = 0
#else
    INTEGER, PARAMETER, PRIVATE  :: mi_dtllevel_base = 2
#endif

    INTEGER, PARAMETER, PRIVATE  :: DIMENSION_1D = 1
    INTEGER, PARAMETER, PRIVATE  :: DIMENSION_2D = 2
    INTEGER, PARAMETER, PRIVATE  :: DIMENSION_3D = 3
contains
    !---------------------------------------------------------------
    subroutine save_terms_toassembly(self)
        IMPLICIT NONE
        TYPE(ASSEMBLY) :: self
        TYPE(FEM) :: ao_FEM
        ! LOCAL
        INTEGER :: li_n
        INTEGER :: li_i
        INTEGER :: li_ref

        CALL printlog("save_terms_toassembly : Start", ai_dtllevel = mi_dtllevel_base  + 2)

        self % opi_patchs_toassembly_tmp        = self % opi_patchs_toassembly
        self % opi_operators_toassembly_tmp     = self % opi_operators_toassembly
        self % opi_matrices_toassembly_tmp      = self % opi_matrices_toassembly
        self % opi_fields_toassembly_tmp        = self % opi_fields_toassembly
        self % opi_norms_toassembly_tmp         = self % opi_norms_toassembly
        self % opi_spaces_toassembly_tmp        = self % opi_spaces_toassembly

        li_n = self % opi_operators_toassembly(0)
        DO li_i=1, li_n
        li_ref = self % opi_operators_toassembly(li_i)
        CALL save_operator_toassembly(ao_FEM, li_ref)
        END DO 

        CALL printlog("save_terms_toassembly : End", ai_dtllevel = mi_dtllevel_base  + 2)

    end subroutine save_terms_toassembly
    !---------------------------------------------------------------
    subroutine load_terms_toassembly(self)
        IMPLICIT NONE
        TYPE(ASSEMBLY) :: self
        TYPE(FEM) :: ao_FEM
        ! LOCAL
        INTEGER :: li_n
        INTEGER :: li_i
        INTEGER :: li_ref

        CALL printlog("load_terms_toassembly : Start", ai_dtllevel = mi_dtllevel_base  + 2)

        self % opi_patchs_toassembly    = self % opi_patchs_toassembly_tmp
        self % opi_operators_toassembly = self % opi_operators_toassembly_tmp
        self % opi_matrices_toassembly  = self % opi_matrices_toassembly_tmp
        self % opi_fields_toassembly    = self % opi_fields_toassembly_tmp
        self % opi_norms_toassembly     = self % opi_norms_toassembly_tmp
        self % opi_spaces_toassembly    = self % opi_spaces_toassembly_tmp

        li_n = self % opi_operators_toassembly(0)
        DO li_i=1, li_n
        li_ref = self % opi_operators_toassembly(li_i)
        CALL load_operator_toassembly(ao_FEM, li_ref)
        END DO 

        CALL printlog("load_terms_toassembly : End", ai_dtllevel = mi_dtllevel_base  + 2)

    end subroutine load_terms_toassembly
    !---------------------------------------------------------------
    subroutine set_spaces_toassembly(self, ao_FEM, ai_grids)
        IMPLICIT NONE
        TYPE(ASSEMBLY) :: self
        TYPE(FEM) :: ao_FEM
        INTEGER :: ai_grids
        ! LOCAL
        INTEGER :: li_nspaces
        INTEGER :: li_id

        CALL printlog("set_spaces_toassembly : Start", ai_dtllevel = mi_dtllevel_base  + 2)

        li_nspaces = 0        
        self % opi_spaces_toassembly = 0
        DO li_id = 0, ao_FEM % oi_nSpaces-1
            if ( ao_FEM % opi_InfoSpace ( li_id, INFOSPACE_TOASSEMBLY) == 1 ) then
                li_nspaces = li_nspaces + 1
                self % opi_spaces_toassembly(li_nspaces) = li_id
            end if
        END DO
        self % opi_spaces_toassembly(0) = li_nspaces

        CALL printlog("set_spaces_toassembly : End", ai_dtllevel = mi_dtllevel_base  + 2)

    end subroutine set_spaces_toassembly
    !---------------------------------------------------------------
    subroutine set_patchs_toassembly(self, ao_FEM, ai_grids)
        IMPLICIT NONE
        TYPE(ASSEMBLY) :: self
        TYPE(FEM) :: ao_FEM
        INTEGER :: ai_grids
        ! LOCAL
        INTEGER :: li_npatchs
        INTEGER :: li_id

        CALL printlog("set_patchs_toassembly : Start", ai_dtllevel = mi_dtllevel_base  + 2)

        li_npatchs = 0        
        self % opi_patchs_toassembly = 0
        DO li_id = 0, ao_FEM % opi_infoGrids (ai_grids, INFOGRIDS_NPATCHS)-1
            if ( ao_FEM % opi_Infopatch ( ai_grids, li_id, INFOPATCH_TOASSEMBLY) == 1 ) then
                li_npatchs = li_npatchs + 1
                self % opi_patchs_toassembly(li_npatchs) = li_id
            end if
        END DO
        self % opi_patchs_toassembly(0) = li_npatchs

        CALL printlog("set_patchs_toassembly : End", ai_dtllevel = mi_dtllevel_base  + 2)

    end subroutine set_patchs_toassembly
    !---------------------------------------------------------------
    subroutine set_elts_toassembly(self, ao_FEM, ai_grids_id, api_values, ai_size )
        IMPLICIT NONE
        TYPE(ASSEMBLY) :: self
        TYPE(FEM) :: ao_FEM
        INTEGER, INTENT(IN)  :: ai_grids_id
        INTEGER, INTENT(IN)  :: ai_size
        INTEGER , DIMENSION(ai_size), INTENT(IN)  :: api_values
        ! LOCAL
        INTEGER :: li_nelts
        INTEGER :: li_patch_id
        INTEGER :: li_i

        CALL printlog("set_elts_toassembly : Start", ai_dtllevel = mi_dtllevel_base  + 2)

        self % opi_elts_toassembly = 0
        self % opi_elts_toassembly(0) = -1
        IF ( ai_size > 0 ) THEN
!            IF ( self % opi_patchs_toassembly(0) /= 1 ) THEN
!                PRINT *, "ERROR set_elts_toassembly: Can be used for one grids and one patch"
!                STOP
!            END IF
            self % opi_elts_toassembly(1:ai_size) = api_values(1:ai_size)
            self % opi_elts_toassembly(0) = ai_size
        END IF

        CALL printlog("set_elts_toassembly : End", ai_dtllevel = mi_dtllevel_base  + 2)

    end subroutine set_elts_toassembly
    !---------------------------------------------------------------
    subroutine set_matrices_toassembly(self, ao_FEM)
        IMPLICIT NONE
        TYPE(ASSEMBLY) :: self
        TYPE(FEM) :: ao_FEM
        ! LOCAL
        INTEGER :: li_n
        INTEGER :: li_id

        CALL printlog("set_matrices_toassembly: Start", ai_dtllevel = mi_dtllevel_base  + 2)

        li_n = 0
        self % opi_matrices_toassembly = 0
        DO li_id = 0, ao_FEM % oi_nmatrices-1
            if ( ao_FEM % opi_InfoMatrix ( li_id, INFOMATRIX_TOASSEMBLY) == 1 ) then
                li_n = li_n + 1
                self % opi_matrices_toassembly(li_n) = li_id
            end if
        END DO
        self % opi_matrices_toassembly(0) = li_n

        CALL printlog("set_matrices_toassembly: End", ai_dtllevel = mi_dtllevel_base  + 2)

    end subroutine set_matrices_toassembly    
    !---------------------------------------------------------------
    subroutine set_operators_toassembly(self, ao_FEM)
        IMPLICIT NONE
        TYPE(ASSEMBLY) :: self
        TYPE(FEM) :: ao_FEM
        ! LOCAL
        INTEGER :: li_noperators
        INTEGER :: li_id

        CALL printlog("set_operators_toassembly : Start", ai_dtllevel = mi_dtllevel_base  + 2)

        li_noperators = 0
        self % opi_operators_toassembly = 0
        DO li_id = 0, ao_FEM % oi_noperators-1
            if ( ao_FEM % opi_InfoOperator ( li_id, INFOOPERATOR_TOASSEMBLY) == 1 ) then
                li_noperators = li_noperators + 1
                self % opi_operators_toassembly(li_noperators) = li_id
            end if
        END DO
        self % opi_operators_toassembly(0) = li_noperators

        CALL printlog("set_operators_toassembly : End", ai_dtllevel = mi_dtllevel_base  + 2)

    end subroutine set_operators_toassembly
    !---------------------------------------------------------------
    subroutine set_fields_toassembly(self, ao_FEM)
        IMPLICIT NONE
        TYPE(ASSEMBLY) :: self
        TYPE(FEM) :: ao_FEM
        ! LOCAL
        INTEGER :: li_nfields
        INTEGER :: li_id

        CALL printlog("set_fields_toassembly : Start", ai_dtllevel = mi_dtllevel_base  + 2)

        li_nfields = 0
        self % opi_fields_toassembly = 0
        DO li_id = 0, ao_FEM % oi_nFields-1
            if ( ao_FEM % opi_InfoField ( li_id, INFOFIELD_TOASSEMBLY) == 1 ) then
                li_nfields = li_nfields + 1
                self % opi_fields_toassembly(li_nfields) = li_id
            end if
        END DO
        self % opi_fields_toassembly(0) = li_nfields

        CALL printlog("set_fields_toassembly : End", ai_dtllevel = mi_dtllevel_base  + 2)

    end subroutine set_fields_toassembly
    !---------------------------------------------------------------
    subroutine set_norms_toassembly(self, ao_FEM)
        IMPLICIT NONE
        TYPE(ASSEMBLY) :: self
        TYPE(FEM) :: ao_FEM
        ! LOCAL
        INTEGER :: li_nnorms
        INTEGER :: li_id

        CALL printlog("set_norms_toassembly : Start", ai_dtllevel = mi_dtllevel_base  + 2)

        li_nnorms = 0
        self % opi_norms_toassembly = 0
        DO li_id = 0, ao_FEM % oi_nNorms-1
            if ( ao_FEM % opi_InfoNorm ( li_id, INFONORM_TOASSEMBLY) == 1 ) then
                li_nnorms = li_nnorms + 1
                self % opi_norms_toassembly(li_nnorms) = li_id
                ! REset the value of the norm
                ao_FEM % opo_N (li_id) % opr_values (:,:) = 0.0_wp
            end if
        END DO
        self % opi_norms_toassembly(0) = li_nnorms

        CALL printlog("set_norms_toassembly : End", ai_dtllevel = mi_dtllevel_base  + 2)

    end subroutine set_norms_toassembly
    !---------------------------------------------------------------
    subroutine select_fields(self, ao_FEM, ai_TYPE, api_fields_id)
        !> this routine selects the fields whom TYPE is ai_TYPE
        IMPLICIT NONE
        TYPE(ASSEMBLY) :: self
        TYPE(FEM) :: ao_FEM
        INTEGER :: ai_TYPE
        INTEGER, DIMENSION(0:), INTENT(INOUT) :: api_fields_id
        ! LOCAL
        INTEGER :: li_nfields
        INTEGER :: li_field
        INTEGER :: li_ref
        INTEGER :: li_current

        CALL printlog("select_fields : Start", ai_dtllevel = mi_dtllevel_base  + 2)

        li_nfields = self % opi_fields_toassembly(0)

        li_current = 0
        DO li_ref = 1, li_nfields

            li_field = self % opi_fields_toassembly(li_ref)

            IF ( ao_FEM % opi_InfoField(li_field, INFOFIELD_TYPE) == ai_TYPE ) THEN
                li_current = li_current + 1
                api_fields_id(li_current) = li_field
            END IF

        END DO
        api_fields_id(0) = li_current

        CALL printlog("select_fields : End", ai_dtllevel = mi_dtllevel_base  + 2)

    end subroutine select_fields
    !---------------------------------------------------------------
    subroutine set_flag_points_basis_assembly(self, ao_FEM)
        IMPLICIT NONE
        TYPE(ASSEMBLY) :: self
        TYPE(FEM) :: ao_FEM
        ! LOCAL
        INTEGER :: li_ref
        INTEGER :: li_nmatrices
        INTEGER :: li_matrix
        INTEGER :: li_nfields
        INTEGER :: li_field
        INTEGER :: li_nnorms
        INTEGER :: li_norm
        LOGICAL :: ll_assembly_basis_m
        LOGICAL :: ll_assembly_points_m
        LOGICAL :: ll_assembly_basis_f
        LOGICAL :: ll_assembly_points_f
        LOGICAL :: ll_assembly_basis_n
        LOGICAL :: ll_assembly_points_n

        CALL printlog("set_flag_points_basis_assembly : Start", ai_dtllevel = mi_dtllevel_base  + 2)

        ll_assembly_basis_m  = .FALSE.
        ll_assembly_points_m = .FALSE.
        li_nmatrices = self % opi_matrices_toassembly(0)
        DO li_ref = 1, li_nmatrices

            li_matrix = self % opi_matrices_toassembly(li_ref)

            ll_assembly_basis_m  = .TRUE.
            ll_assembly_points_m = .TRUE.

        END DO

        ll_assembly_basis_f  = .FALSE.
        ll_assembly_points_f = .FALSE.
        li_nfields = self % opi_fields_toassembly(0)
        DO li_ref = 1, li_nfields

            li_field = self % opi_fields_toassembly(li_ref)
            
            ll_assembly_basis_f  = .TRUE.
            ll_assembly_points_f = .TRUE.

        END DO

        ll_assembly_basis_n  = .FALSE.
        ll_assembly_points_n = .FALSE.
        li_nnorms = self % opi_norms_toassembly(0)
        DO li_ref = 1, li_nnorms

            li_norm = self % opi_norms_toassembly(li_ref)

            ll_assembly_basis_n  = .TRUE.
            ll_assembly_points_n = .TRUE.

        END DO

        self % ol_assembly_basis  = ll_assembly_basis_m  .OR. ll_assembly_basis_f  .OR. ll_assembly_basis_n
        self % ol_assembly_points = ll_assembly_points_m .OR. ll_assembly_points_f .OR. ll_assembly_points_n

        CALL printlog("set_flag_points_basis_assembly : End", ai_dtllevel = mi_dtllevel_base  + 2)

    end subroutine set_flag_points_basis_assembly
    !---------------------------------------------------------------
    subroutine set_operators_nderiv( self, ao_FEM )
    ! TODO to be done in python
        IMPLICIT NONE
        TYPE(ASSEMBLY) :: self
        TYPE(FEM) :: ao_FEM
        ! LOCAL
        INTEGER :: li_id

        CALL printlog("set_operators_nderiv : Start", ai_dtllevel = mi_dtllevel_base  + 2)

        DO li_id = 0, ao_FEM % oi_nOperators-1

            select case ( ao_FEM % opi_InfoOperator(li_id, INFOOPERATOR_TYPE))
                case ( MASS )
                    ao_FEM % opi_Infooperator ( li_id, INFOOPERATOR_NDERIV) = 0
                case ( STIFFNESS )
                    ao_FEM % opi_Infooperator ( li_id, INFOOPERATOR_NDERIV) = 1
                case ( ADVECTION )
                    ao_FEM % opi_InfoOperator ( li_id, INFOOPERATOR_NDERIV) = 1
                case ( SECOND_DERIV )
                    ao_FEM % opi_InfoOperator ( li_id, INFOOPERATOR_NDERIV) = 2
                case Default
                    call printlog("The following Operator Not Yet implemented", ai_dtllevel = 0)
                    call concatmsg("id ", ai_dtllevel = 0)
                    call concatmsg(li_id, ai_dtllevel = 0)
                    call printmsg(ai_dtllevel = 0)
            end select

        END DO

        CALL printlog("set_operators_nderiv : End", ai_dtllevel = mi_dtllevel_base  + 2)

    end subroutine set_operators_nderiv
    !---------------------------------------------------------------
    subroutine create_spaces_info( self, ao_FEM )
        IMPLICIT NONE
        TYPE(ASSEMBLY) :: self
        TYPE(FEM) :: ao_FEM
        ! LOCAL
        INTEGER :: li_id
        INTEGER :: li_nderiv
        INTEGER :: li_nderiv_code
        INTEGER :: li_grids
        INTEGER :: li_composed
        INTEGER :: li_dim
        INTEGER :: li_nnz
        INTEGER :: li_npts
        INTEGER :: li_tensorlevel

        CALL printlog("create_spaces_info : Start", ai_dtllevel = mi_dtllevel_base  + 2)

        li_nderiv = MAX(self % oi_matrix_nderiv, self % oi_mapping_nderiv, 2)

!        print*,'self % oi_matrix_nderiv =', self % oi_matrix_nderiv

        DO li_id = 0, ao_FEM % oi_nspaces-1
            
            li_grids = ao_FEM % opi_InfoSpace(li_id, INFOSPACE_GRIDS )

            li_dim = ao_FEM % opi_dim(li_grids)

            li_tensorlevel = ao_FEM % opi_InfoPatch ( li_grids, 0, INFOPATCH_TENSOR)

            ! we must compute the number of derivatives in order to allocate opr_points
            if ( ao_FEM % opi_dim(li_grids) == 1 ) then

                select case ( li_nderiv )
                    case ( 0 )
                        li_nderiv_code = 0
                    case ( 1 )
                        li_nderiv_code = 1
    !                    return
                    case ( 2 )
                        li_nderiv_code = 2
    !                    return
                    case Default
                        print*,"create_spaces_info: Type code Not Yet implemented"
                        return
                end select

            end if
            if ( ao_FEM % opi_dim(li_grids) == 2 ) then

                select case ( li_nderiv )
                    case ( 0 )
                        li_nderiv_code = 0
                    case ( 1 )
                        li_nderiv_code = 2
    !                    return
                    case ( 2 )
                        li_nderiv_code = 5
    !                    return
                    case Default
                        print*,"create_spaces_info: Type code Not Yet implemented"
                        return
                end select

            end if
            if ( ao_FEM % opi_dim(li_grids) == 3 ) then

                select case ( li_nderiv )
                    case ( 0 )
                        li_nderiv_code = 0
                    case ( 1 )
                        li_nderiv_code = 3
    !                    return
                    case ( 2 )
                        li_nderiv_code = 6
    !                    return
                    case Default
                        print*,"create_spaces_info: Type code Not Yet implemented"
                        return
                end select

            end if

!            print*,'li_nderiv =', li_nderiv
            ! c est bien maxnderiv_mapping qu'on utilise, car on en a besoin uniquement pour le mapping
            ! dans le ca sou on utilise le mapping issu de la definition du domaine
!            PRINT *, 'TODO : li_nderiv_code a ete mis a 2'
!            PRINT *, 'opi_dim = ', MAXVAL(ao_FEM % opi_dim(:))
!            PRINT *, 'self % oi_maxnpts = ', self % oi_maxnpts
!            li_nderiv_code = 2

            li_npts = self % oi_maxnpts

!            IF (li_tensorlevel==2) THEN
!                li_npts = assl_get_arbnnz(self % opo_bbox_sp(li_id) % oo_asslx)
!                self % oi_maxnpts = MAX(li_npts, self % oi_maxnpts)
!
!!                IF (self % oi_maxnpts /= li_nnz) THEN
!!                    PRINT *, 'SERIOUS ERROR create_spaces_info: maxnpts must be equal to nnz'
!!                    PRINT *, 'maxnpts = ', self % oi_maxnpts
!!                    PRINT *, 'nnz = ', self % opo_bbox_sp(li_id) % oo_asslx % arbnnz
!!                END IF
!            END IF

            ALLOCATE(self % opo_info_sp (li_id) % opr_points (0:li_nderiv_code, MAXVAL(ao_FEM % opi_Rd(:)), li_npts))
            ALLOCATE(self % opo_info_sp (li_id) % opr_jacobians (li_npts))
            ALLOCATE(self % opo_info_sp (li_id) % opr_invJacobian ( li_dim, li_dim, li_npts))
            ALLOCATE(self % opo_info_sp (li_id) % opr_invJacobian_nn ( li_dim, li_dim, li_npts))
            
        END DO

        CALL printlog("create_spaces_info : End", ai_dtllevel = mi_dtllevel_base  + 2)

    end subroutine create_spaces_info
    !---------------------------------------------------------------
    subroutine free_spaces_info( self, ao_FEM )
        IMPLICIT NONE
        TYPE(ASSEMBLY) :: self
        TYPE(FEM) :: ao_FEM
        ! LOCAL
        INTEGER :: li_id

        CALL printlog("free_spaces_info : Start", ai_dtllevel = mi_dtllevel_base  + 2)

        DO li_id = 0, ao_FEM % oi_nspaces-1

            DEALLOCATE(self % opo_info_sp (li_id) % opr_points)
            DEALLOCATE(self % opo_info_sp (li_id) % opr_jacobians)
            DEALLOCATE(self % opo_info_sp (li_id) % opr_invJacobian)
            DEALLOCATE(self % opo_info_sp (li_id) % opr_invJacobian_nn)
            
        END DO

        CALL printlog("free_spaces_info : End", ai_dtllevel = mi_dtllevel_base  + 2)

    end subroutine free_spaces_info
    !---------------------------------------------------------------
    subroutine create_mappings_info( self, ao_FEM )
        IMPLICIT NONE
        TYPE(ASSEMBLY) :: self
        TYPE(FEM) :: ao_FEM
        ! LOCAL
        INTEGER :: li_id
        INTEGER :: li_nderiv
        INTEGER :: li_nderiv_code
        INTEGER :: li_grids
        INTEGER :: li_space

        CALL printlog("create_mappings_info : Start", ai_dtllevel = mi_dtllevel_base  + 2)

        li_nderiv = MAX(self % oi_matrix_nderiv, self % oi_mapping_nderiv, 2)

        DO li_id = 0, ao_FEM % oi_nmappings-1
            li_space   = ao_FEM % opi_InfoMapping (li_id, INFOMAPPING_SPACE )
            li_grids   = ao_FEM % opi_InfoSpace (li_space, INFOSPACE_GRIDS )

            ! we must compute the number of derivatives in order to allocate opr_points
            if ( ao_FEM % opi_dim(li_grids) == 1 ) then

                select case ( li_nderiv )
                    case ( 0 )
                        li_nderiv_code = 0
                    case ( 1 )
                        li_nderiv_code = 1
    !                    return
                    case ( 2 )
                        li_nderiv_code = 2
    !                    return
                    case Default
                        print*,"create_mappings_info: Type code Not Yet implemented"
                        return
                end select

            end if
            if ( ao_FEM % opi_dim(li_grids) == 2 ) then

                select case ( li_nderiv )
                    case ( 0 )
                        li_nderiv_code = 0
                    case ( 1 )
                        li_nderiv_code = 2
    !                    return
                    case ( 2 )
                        li_nderiv_code = 5
    !                    return
                    case Default
                        print*,"create_mappings_info: Type code Not Yet implemented"
                        return
                end select

            end if
            if ( ao_FEM % opi_dim(li_grids) == 3 ) then

                select case ( li_nderiv )
                    case ( 0 )
                        li_nderiv_code = 0
                    case ( 1 )
                        li_nderiv_code = 3
    !                    return
                    case ( 2 )
                        li_nderiv_code = 6
    !                    return
                    case Default
                        print*,"create_mappings_info: Type code Not Yet implemented"
                        return
                end select

            end if

            ALLOCATE(self % opo_info_mp (li_id) % opr_points (0:li_nderiv_code    &
            , MAXVAL(ao_FEM % opi_Rd(:)), self % oi_maxnpts))
            ALLOCATE(self % opo_info_mp (li_id) % opr_jacobians (self % oi_maxnpts))

            ALLOCATE(self % opo_info_mp (li_id) % opr_invJacobian ( ao_FEM % opi_dim(li_grids)  &
            , ao_FEM % opi_dim(li_grids), self % oi_maxnpts))
            ALLOCATE(self % opo_info_mp (li_id) % opr_invJacobian_nn ( ao_FEM % opi_dim(li_grids)  &
            , ao_FEM % opi_dim(li_grids), self % oi_maxnpts))
        END DO

        CALL printlog("create_mappings_info : End", ai_dtllevel = mi_dtllevel_base  + 2)

    end subroutine create_mappings_info
    !---------------------------------------------------------------
    subroutine free_mappings_info( self, ao_FEM )
        IMPLICIT NONE
        TYPE(ASSEMBLY) :: self
        TYPE(FEM) :: ao_FEM
        ! LOCAL
        INTEGER :: li_id

        CALL printlog("free_mappings_info : Start", ai_dtllevel = mi_dtllevel_base  + 2)

        DO li_id = 0, ao_FEM % oi_nmappings-1
            DEALLOCATE(self % opo_info_mp (li_id) % opr_points)
            DEALLOCATE(self % opo_info_mp (li_id) % opr_jacobians)
            DEALLOCATE(self % opo_info_mp (li_id) % opr_invJacobian_nn)
            DEALLOCATE(self % opo_info_mp (li_id) % opr_invJacobian)
        END DO

        CALL printlog("free_mappings_info : End", ai_dtllevel = mi_dtllevel_base  + 2)

    end subroutine free_mappings_info
    !---------------------------------------------------------------
    subroutine create_grids_info( self, ao_FEM )
        IMPLICIT NONE
        TYPE(ASSEMBLY) :: self
        TYPE(FEM) :: ao_FEM
        ! LOCAL
        INTEGER :: li_id
        INTEGER :: li_nderiv
        INTEGER :: li_grids

        CALL printlog("create_grids_info : Start", ai_dtllevel = mi_dtllevel_base  + 2)        
        
        DO li_id = 0, ao_FEM % oi_ngrids-1
            li_grids   = li_id

            ! we must compute the number of derivatives in order to allocate opr_points
            if ( ao_FEM % opi_dim(li_grids) == 1 ) then

                select case ( MAX(self % oi_matrix_nderiv , self % oi_mapping_nderiv) )
                    case ( 0 )
                        li_nderiv = 0
                    case ( 1 )
                        li_nderiv = 1
    !                    return
                    case ( 2 )
                        li_nderiv = 2
    !                    return
                    case Default
                        print*,"create_mappings_info: Type code Not Yet implemented"
                        return
                end select

            end if
            if ( ao_FEM % opi_dim(li_grids) == 2 ) then

                select case ( MAX(self % oi_matrix_nderiv , self % oi_mapping_nderiv) )
                    case ( 0 )
                        li_nderiv = 0
                    case ( 1 )
                        li_nderiv = 2
    !                    return
                    case ( 2 )
                        li_nderiv = 5
    !                    return
                    case Default
                        print*,"create_mappings_info: Type code Not Yet implemented"
                        return
                end select

            end if
            if ( ao_FEM % opi_dim(li_grids) == 3 ) then

                select case ( MAX(self % oi_matrix_nderiv , self % oi_mapping_nderiv))
                    case ( 0 )
                        li_nderiv = 0
                    case ( 1 )
                        li_nderiv = 3
    !                    return
                    case ( 2 )
                        li_nderiv = 6
    !                    return
                    case Default
                        print*,"create_mappings_info: Type code Not Yet implemented"
                        return
                end select

            end if

!            li_nderiv = ao_FEM % opi_dim(li_grids)

            ALLOCATE(self % opo_info_gr (li_id) % opr_points (0:li_nderiv    &
            , MAXVAL(ao_FEM % opi_Rd(:)), self % oi_maxnpts))
            ALLOCATE(self % opo_info_gr (li_id) % opr_jacobians (self % oi_maxnpts))

            ALLOCATE(self % opo_info_gr (li_id) % opr_invJacobian ( ao_FEM % opi_dim(li_grids)  &
            , ao_FEM % opi_dim(li_grids), self % oi_maxnpts))
            ALLOCATE(self % opo_info_gr (li_id) % opr_invJacobian_nn ( ao_FEM % opi_dim(li_grids)  &
            , ao_FEM % opi_dim(li_grids), self % oi_maxnpts))
        END DO

        CALL printlog("create_grids_info : End", ai_dtllevel = mi_dtllevel_base  + 2)

    end subroutine create_grids_info
    !---------------------------------------------------------------
    subroutine free_grids_info( self, ao_FEM )
        IMPLICIT NONE
        TYPE(ASSEMBLY) :: self
        TYPE(FEM) :: ao_FEM
        ! LOCAL
        INTEGER :: li_id

        CALL printlog("free_grids_info : Start", ai_dtllevel = mi_dtllevel_base  + 2)

        DO li_id = 0, ao_FEM % oi_ngrids-1
            DEALLOCATE(self % opo_info_gr (li_id) % opr_points)
            DEALLOCATE(self % opo_info_gr (li_id) % opr_jacobians)
            DEALLOCATE(self % opo_info_gr (li_id) % opr_invJacobian_nn)
            DEALLOCATE(self % opo_info_gr (li_id) % opr_invJacobian)
        END DO

        CALL printlog("free_grids_info : End", ai_dtllevel = mi_dtllevel_base  + 2)

    end subroutine free_grids_info
    !---------------------------------------------------------------
    subroutine create_physical_basis( self, ao_FEM )
        IMPLICIT NONE
        TYPE(ASSEMBLY) :: self
        TYPE(FEM) :: ao_FEM
        ! LOCAL
        INTEGER :: li_id
        INTEGER :: li_nderiv
        INTEGER :: li_grids
        INTEGER :: li_dim
        INTEGER :: li_maxnen
        INTEGER :: li_maxnpts
        INTEGER :: li_nH

        CALL printlog("create_physical_basis: Start", ai_dtllevel = mi_dtllevel_base  + 2)

        li_nderiv = MAX(self % oi_matrix_nderiv, self % oi_mapping_nderiv, 2)

        DO li_id = 0, ao_FEM % oi_nspaces-1
            
            li_grids = ao_FEM % opi_InfoSpace(li_id, INFOSPACE_GRIDS )

            li_dim = ao_FEM % opi_dim(li_grids)
            li_maxnpts = self % oi_maxnpts
            li_maxnen = ao_FEM % opo_spaces (li_id) % oi_maxnen
            li_nH = li_dim * (li_dim + 1) / 2
            
            ALLOCATE(self % opo_pBasis (li_id) % opr_B        (          1:li_maxnen, 1:li_maxnpts))
            ALLOCATE(self % opo_pBasis (li_id) % opr_gradB    (1:li_dim, 1:li_maxnen, 1:li_maxnpts))
            ALLOCATE(self % opo_pBasis (li_id) % opr_curlB    (1:li_dim, 1:li_maxnen, 1:li_maxnpts))
            ALLOCATE(self % opo_pBasis (li_id) % opr_HessianB (1:li_nH , 1:li_maxnen, 1:li_maxnpts))

        END DO

        CALL printlog("create_physical_basis : End", ai_dtllevel = mi_dtllevel_base  + 2)

    end subroutine create_physical_basis
    !---------------------------------------------------------------
    subroutine free_physical_basis( self, ao_FEM )
        IMPLICIT NONE
        TYPE(ASSEMBLY) :: self
        TYPE(FEM) :: ao_FEM
        ! LOCAL
        INTEGER :: li_id

        CALL printlog("free_physical_basis : Start", ai_dtllevel = mi_dtllevel_base  + 2)

        DO li_id = 0, ao_FEM % oi_nspaces-1

            DEALLOCATE(self % opo_pBasis (li_id) % opr_B)
            DEALLOCATE(self % opo_pBasis (li_id) % opr_gradB)
            DEALLOCATE(self % opo_pBasis (li_id) % opr_curlB)
            DEALLOCATE(self % opo_pBasis (li_id) % opr_HessianB)

        END DO

        CALL printlog("free_physical_basis: End", ai_dtllevel = mi_dtllevel_base  + 2)

    end subroutine free_physical_basis

    !---------------------------------------------------------------
    subroutine create_spaces_stored_data( self, ao_FEM )
        IMPLICIT NONE
        TYPE(ASSEMBLY) :: self
        TYPE(FEM) :: ao_FEM
        ! LOCAL
        INTEGER :: li_space
        INTEGER :: li_grids
        INTEGER :: li_npatchs
        INTEGER :: li_patch
        INTEGER :: li_dim
        INTEGER :: li_composed
        INTEGER, DIMENSION(:), POINTER :: lpi_nderiv

        CALL printlog("create_spaces_stored_data : Start", ai_dtllevel = mi_dtllevel_base  + 2)

        DO li_space = 0, ao_FEM % oi_nspaces-1
            ! if it is a composed space, we do not need to defin the geometry
!            li_composed = ao_FEM % opi_InfoSpace(li_space, INFOSPACE_COMPOSED)
!            IF (li_composed == 1) THEN
!                CYCLE
!            END IF

            if ( ao_FEM % opi_InfoSpace (li_space, INFOSPACE_STOREDDATA ) == 0 ) then
                cycle
            end if

            li_grids   = ao_FEM % opi_InfoSpace (li_space, INFOSPACE_GRIDS )
            li_npatchs = ao_FEM % opi_infoGrids (li_grids, INFOGRIDS_NPATCHS)
            li_dim = ao_FEM % opi_dim(li_grids)
            ALLOCATE(lpi_nderiv(li_dim))
            lpi_nderiv = self % oi_matrix_nderiv
            DO li_patch = 0,li_npatchs-1
                CALL create_stored_data(ao_FEM % opo_grids (li_grids) % opo_grid(li_patch) &
                , ao_FEM % opo_spaces ( li_space ) % oo_mapping, li_patch+1, lpi_nderiv)
            END DO
            DEALLOCATE(lpi_nderiv)

        END DO

        CALL printlog("create_spaces_stored_data : End", ai_dtllevel = mi_dtllevel_base  + 2)

    end subroutine create_spaces_stored_data
    !---------------------------------------------------------------
    subroutine free_spaces_stored_data( self, ao_FEM )
        IMPLICIT NONE
        TYPE(ASSEMBLY) :: self
        TYPE(FEM) :: ao_FEM
        ! LOCAL
        INTEGER :: li_space
        INTEGER :: li_grids
        INTEGER :: li_npatchs
        INTEGER :: li_patch
        INTEGER :: li_composed

        CALL printlog("free_spaces_stored_data : Start", ai_dtllevel = mi_dtllevel_base  + 2)

        DO li_space = 0, ao_FEM % oi_nspaces-1

            if ( ao_FEM % opi_InfoSpace (li_space, INFOSPACE_STOREDDATA ) == 0 ) then
                cycle
            end if

            li_grids   = ao_FEM % opi_InfoSpace (li_space, INFOSPACE_GRIDS )
            li_npatchs = ao_FEM % opi_infoGrids (li_grids, INFOGRIDS_NPATCHS)

            DO li_patch = 0,li_npatchs-1
                CALL free_stored_data(ao_FEM % opo_grids (li_grids) % opo_grid(li_patch))
            END DO

        END DO

        CALL printlog("free_spaces_stored_data : End", ai_dtllevel = mi_dtllevel_base  + 2)

    end subroutine free_spaces_stored_data
    !---------------------------------------------------------------
    subroutine create_mappings_stored_data( self, ao_FEM )
        IMPLICIT NONE
        TYPE(ASSEMBLY) :: self
        TYPE(FEM) :: ao_FEM
        ! LOCAL
        INTEGER :: li_mapping
        INTEGER :: li_grids
        INTEGER :: li_npatchs
        INTEGER :: li_space
        INTEGER :: li_patch
        INTEGER :: li_dim
        INTEGER, DIMENSION(:), POINTER :: lpi_nderiv

        CALL printlog("create_mappings_stored_data : Start", ai_dtllevel = mi_dtllevel_base  + 2)

        DO li_mapping = 0, ao_FEM % oi_nmappings-1
            if ( ao_FEM % opi_InfoMapping (li_mapping, INFOMAPPING_STOREDDATA ) == 0 ) then
                cycle
            end if

            li_space   = ao_FEM % opi_InfoMapping (li_mapping, INFOMAPPING_SPACE )
            li_grids   = ao_FEM % opi_InfoSpace (li_space, INFOSPACE_GRIDS )
            li_npatchs = ao_FEM % opi_InfoGrids(li_grids, INFOGRIDS_NPATCHS)

            li_dim = ao_FEM % opi_dim(li_grids)
            ALLOCATE(lpi_nderiv(li_dim))
            lpi_nderiv = self % oi_mapping_nderiv
            DO li_patch = 0,li_npatchs-1
                CALL create_stored_data(ao_FEM % opo_grids (li_grids) % opo_grid(li_patch) &
                , ao_FEM % opo_mappings ( li_mapping ) % oo_mapping, li_patch+1, lpi_nderiv)
            END DO
            DEALLOCATE(lpi_nderiv)

        END DO

        CALL printlog("create_mappings_stored_data : End", ai_dtllevel = mi_dtllevel_base  + 2)

    end subroutine create_mappings_stored_data
    !---------------------------------------------------------------
    subroutine free_mappings_stored_data( self, ao_FEM )
        IMPLICIT NONE
        TYPE(ASSEMBLY) :: self
        TYPE(FEM) :: ao_FEM
        ! LOCAL
        INTEGER :: li_mapping
        INTEGER :: li_grids
        INTEGER :: li_patch
        INTEGER :: li_npatchs
        INTEGER :: li_space

        CALL printlog("free_mappings_stored_data : Start", ai_dtllevel = mi_dtllevel_base  + 2)

        DO li_mapping = 0, ao_FEM % oi_nmappings-1
            if ( ao_FEM % opi_InfoMapping (li_mapping, INFOMAPPING_STOREDDATA ) == 0 ) then
                cycle
            end if

            li_space   = ao_FEM % opi_InfoMapping (li_mapping, INFOMAPPING_SPACE )
            li_grids   = ao_FEM % opi_InfoSpace (li_space, INFOSPACE_GRIDS )
            li_npatchs = ao_FEM % opi_InfoGrids(li_grids, INFOGRIDS_NPATCHS)

            DO li_patch = 0,li_npatchs-1
                CALL free_stored_data(ao_FEM % opo_grids (li_grids) % opo_grid(li_patch))
            END DO

        END DO

        CALL printlog("free_mappings_stored_data : End", ai_dtllevel = mi_dtllevel_base  + 2)

    end subroutine free_mappings_stored_data
    !---------------------------------------------------------------
    subroutine create_spaces_bbox( self, ao_FEM )
        IMPLICIT NONE
        TYPE(ASSEMBLY) :: self
        TYPE(FEM) :: ao_FEM
        ! LOCAL
        INTEGER :: li_id
        INTEGER :: li_maxp
        INTEGER :: li_maxnen
        INTEGER :: li_dirmaxnpts
        INTEGER :: li_npatchs
        INTEGER :: li_grids
        INTEGER :: li_composed
        INTEGER :: li_nderiv
        INTEGER :: li_tensorlevel
        INTEGER :: li_dim
        INTEGER :: li_n
        INTEGER :: li_nel
        INTEGER :: li_d
        INTEGER :: li_patch_id
        TYPE(GEOMETRIES), POINTER :: lp_geos
        INTEGER, DIMENSION(:), POINTER :: lpi_k, lpi_g

        CALL printlog("create_spaces_bbox : Start", ai_dtllevel = mi_dtllevel_base  + 1)

        DO li_id = 0, ao_FEM % oi_nspaces-1

            lp_geos => ao_FEM % opo_spaces (li_id) % oo_mapping

            li_maxp = ao_FEM % opo_spaces (li_id) % oi_maxp
            li_maxnen = ao_FEM % opo_spaces (li_id) % oi_maxnen

            li_tensorlevel = ao_fem % opi_InfoSpace ( li_id, INFOSPACE_TENSOR ) 

            li_grids   = ao_FEM % opi_InfoSpace(li_id, INFOSPACE_GRIDS )
            li_npatchs = ao_FEM % opi_infoGrids (li_grids, INFOGRIDS_NPATCHS)
            li_dirmaxnpts = MAXVAL (ao_FEM % opi_InfoPatch ( li_grids, 0:li_npatchs-1, INFOPATCH_DIRMAXNPTS) )

            li_dim = ao_FEM % opi_dim(li_grids)

            !> \todo il faut faire un REALloc si le nderiv change
            li_nderiv = MAX(self % oi_matrix_nderiv, self % oi_mapping_nderiv)
!            li_nderiv = MAX(self % oi_matrix_nderiv, self % oi_mapping_nderiv, 2)

        ! *************************************************
        !                 TENSOR-LEVEL   1
        ! *************************************************
            IF (li_tensorlevel==1) THEN
                CALL create_bbox(self % opo_bbox_sp(li_id)  &
                , li_dim    &
                , li_nderiv       &
                , li_maxnen                     &
                , self % oi_maxnpts             &
    !            , ai_maxnpts = li_dirmaxnpts    &
                , li_dirmaxnpts    &
    !            , ai_maxp = li_maxp             &
                , li_maxp )
            END IF
        ! *************************************************

        ! *************************************************
        !                 TENSOR-LEVEL   2
        ! *************************************************
!            IF (li_tensorlevel==2) THEN
!                ALLOCATE(lpi_k(li_dim))
!                ALLOCATE(lpi_g(li_dim))
!                !> \todo ne marchera pas si on a plusieurs patchs, a modifier
!                li_patch_id = 0
!                li_n = ao_FEM % opo_spaces (li_id) % oo_mapping % opo_geo (li_patch_id+1) % opi_ndiffu(1) - 1
!
!                li_nel = 1
!                DO li_d = 2, li_dim
!                    li_nel = li_nel * ( ao_FEM % opo_spaces (li_id) % oo_mapping % opo_geo (li_patch_id+1) % opi_ndiffu(li_d) - 1 )
!                END DO
!#ifdef _MYDEBUG
!                print *, 'create bbox with ', li_n, li_nel
!#endif
!
!#ifdef _TRACE
!        CALL printlog("TODO : a mettre les bonnes valeurs pour  li_n (=n1) et li_nel (n2*n3)" &
!        , ai_dtllevel = mi_dtllevel_base + 1)
!#endif
!                lpi_k = lp_geos % opo_geo(li_patch_id+1) % opi_P + 1
!                !> \todo a changer avec le vrai nombre de pts de quad
!!                lpi_g = lp_geos % opo_geo(li_patch_id+1) % opi_P + 1
!
!                lpi_g = ao_FEM % opo_grids(li_grids) % opo_grid (li_patch_id) % opo_elts (1) % opi_npts (:)
!                lpi_g(1) = lpi_g(1) / li_n
!                
!#ifdef _MYDEBUG
!                print *, 'li_n, lpi_k, lpi_g =', li_n, lpi_k, lpi_g
!#endif
!                
!                CALL create_bbox(self % opo_bbox_sp(li_id) &
!                , ao_FEM % opi_dim(li_grids) &
!                , li_nderiv &
!                , li_n &
!                , li_nel &
!                , lpi_k &
!                , lpi_g)
!                DEALLOCATE(lpi_k, lpi_g)
!            END IF
        ! *************************************************

!            print*,'li_maxp =', li_maxp
!            print*,'li_maxnen =', li_maxnen
!            print*,'self % oi_maxnpts =', self % oi_maxnpts
!            print*,'li_dirmaxnpts =', li_dirmaxnpts
!            print*,'li_nderiv =', li_nderiv
!            print*,'ao_FEM % opi_dim(li_grids) =', ao_FEM % opi_dim(li_grids)

        END DO

        CALL printlog("create_spaces_bbox : End", ai_dtllevel = mi_dtllevel_base  + 1)

    end subroutine create_spaces_bbox
    !---------------------------------------------------------------
    subroutine free_spaces_bbox( self, ao_FEM )
        IMPLICIT NONE
        TYPE(ASSEMBLY) :: self
        TYPE(FEM) :: ao_FEM
        ! LOCAL
        INTEGER :: li_id
        INTEGER :: li_composed

        CALL printlog("free_spaces_bbox : Start", ai_dtllevel = mi_dtllevel_base  + 1)

        DO li_id = 0, ao_FEM % oi_nspaces-1
            CALL free_bbox(self % opo_bbox_sp(li_id))
        END DO

        CALL printlog("free_spaces_bbox : End", ai_dtllevel = mi_dtllevel_base  + 1)

    end subroutine free_spaces_bbox
    !---------------------------------------------------------------
    subroutine create_mappings_bbox( self, ao_FEM )
        IMPLICIT NONE
        TYPE(ASSEMBLY) :: self
        TYPE(FEM) :: ao_FEM
        ! LOCAL
        INTEGER :: li_mapping
        INTEGER :: li_maxp
        INTEGER :: li_maxnen
        INTEGER :: li_dirmaxnpts
        INTEGER :: li_npatchs
        INTEGER :: li_grids
        INTEGER :: li_space
        INTEGER :: li_tensorlevel
        INTEGER :: li_dim
        INTEGER :: li_n
        INTEGER, DIMENSION(:), POINTER :: lpi_k, lpi_g
        LOGICAL :: ll_continue

        CALL printlog("create_mappings_bbox : Start", ai_dtllevel = mi_dtllevel_base  + 1)

        DO li_mapping = 0, ao_FEM % oi_nmappings-1

            li_maxp = ao_FEM % opo_mappings (li_mapping) % oi_maxp
            li_maxnen = ao_FEM % opo_mappings (li_mapping) % oi_maxnen

            li_tensorlevel =  ao_fem % opi_InfoMapping ( li_mapping, INFOMAPPING_TENSOR )

            li_space   = ao_FEM % opi_InfoMapping (li_mapping, INFOMAPPING_SPACE )
            li_grids   = ao_FEM % opi_InfoSpace (li_space, INFOSPACE_GRIDS )
            li_npatchs = ao_FEM % opi_InfoGrids(li_grids, INFOGRIDS_NPATCHS)
            
            li_dirmaxnpts = MAXVAL (ao_FEM % opi_InfoPatch ( li_grids, 0:li_npatchs-1, INFOPATCH_DIRMAXNPTS) )

            li_dim = ao_FEM % opi_dim(li_grids)

#ifdef _DEBUG
            print*,'current mapping id',li_mapping
            print*,'li_grids=',li_grids
            print*,'li_npatchs=',li_npatchs
            print*,'li_dirmaxnpts=',li_dirmaxnpts
            print*,'ao_FEM % opi_dim(li_grids)=',ao_FEM % opi_dim(li_grids)
            print*,'li_maxnen=',li_maxnen
            print*,'self % oi_maxnpts=',self % oi_maxnpts
            print*,'self % oi_mapping_nderiv=',self % oi_mapping_nderiv
            print*,'li_maxp=',li_maxp
            print*,'li_tensorlevel=',li_tensorlevel
#endif
            IF (li_tensorlevel==1) THEN
                CALL create_bbox(self % opo_bbox_mp(li_mapping)  &
                , li_dim    &
                , self % oi_mapping_nderiv      &
                , li_maxnen                     &
                , self % oi_maxnpts             &
    !            , ai_maxnpts = li_dirmaxnpts    &
                , li_dirmaxnpts    &
    !            , ai_maxp = li_maxp             &
                , li_maxp   )
            END IF

            IF (li_tensorlevel==2) THEN
                print *, 'par ici'
                ALLOCATE(lpi_k(li_dim))
                ALLOCATE(lpi_g(li_dim))
                li_n = 5
                lpi_k = 2 + 1
                lpi_g = 2 + 1
!                CALL create_bbox_t2(self % opo_bbox_sp(li_id) &
!                , ao_FEM % opi_dim(li_grids) &
!                , li_nderiv &
!                , lpi_n &
!                , lpi_k &
!                , lpi_g)

                DEALLOCATE(lpi_k, lpi_g)
            END IF

        END DO

        CALL printlog("create_mappings_bbox : End", ai_dtllevel = mi_dtllevel_base  + 1)

    end subroutine create_mappings_bbox
    !---------------------------------------------------------------
    subroutine free_mappings_bbox( self, ao_FEM )
        IMPLICIT NONE
        TYPE(ASSEMBLY) :: self
        TYPE(FEM) :: ao_FEM
        ! LOCAL
        INTEGER :: li_id

        CALL printlog("free_mappings_bbox : Start", ai_dtllevel = mi_dtllevel_base  + 1)

        DO li_id = 0, ao_FEM % oi_nmappings-1
            CALL free_bbox(self % opo_bbox_mp(li_id))
        END DO

        CALL printlog("free_mappings_bbox : End", ai_dtllevel = mi_dtllevel_base  + 1)

    end subroutine free_mappings_bbox
    !---------------------------------------------------------------
    subroutine determinant_1D( apr_points, ai_npts, apr_Jacobian, apr_invJacobian, apr_invJacobian_nn )
        !> evaluate xi_1 
        IMPLICIT NONE
        REAL(wp), DIMENSION(0:,:,:) :: apr_points
        INTEGER :: ai_npts
        REAL(wp), DIMENSION(:), intent(inout) :: apr_Jacobian
        REAL(wp), DIMENSION(:,:,:), intent(inout) :: apr_invJacobian
        REAL(wp), DIMENSION(:,:,:), intent(inout) :: apr_invJacobian_nn
        ! LOCAL
        INTEGER, PARAMETER :: li_X = 1
        INTEGER, PARAMETER :: li_ksi = 1
        INTEGER :: li_i
        INTEGER :: li_Rd = 3

        apr_jacobian(1:ai_npts) = 0.0_wp
        DO li_i = 1, li_Rd
!        DO li_i = 1, 2 
!        DO li_i = 1, 1 
            apr_jacobian(1:ai_npts)= apr_jacobian(1:ai_npts) + apr_points(li_ksi,li_i,1:ai_npts)**2
        END DO
!        print *, 'jacobian-1D ', apr_jacobian
        apr_jacobian(1:ai_npts) = SQRT(apr_jacobian(1:ai_npts))
        apr_invJacobian(li_ksi,li_X,1:ai_npts)= 1.0_wp / apr_jacobian(1:ai_npts)
        apr_invJacobian_nn(li_ksi,li_X,1:ai_npts)= 1.0_wp

    end subroutine determinant_1D    
    !---------------------------------------------------------------
    subroutine determinant_2D( apr_points, ai_npts, apr_Jacobian, apr_invJacobian, apr_invJacobian_nn )
        !> evaluate xi_1 xj_2 - xi_2 * xj_1
        IMPLICIT NONE
        REAL(wp), DIMENSION(0:,:,:) :: apr_points
        INTEGER :: ai_npts
        REAL(wp), DIMENSION(:), intent(inout) :: apr_Jacobian
        REAL(wp), DIMENSION(:,:,:), intent(inout) :: apr_invJacobian
        REAL(wp), DIMENSION(:,:,:), intent(inout) :: apr_invJacobian_nn
        ! LOCAL
        INTEGER, PARAMETER :: li_X = 1
        INTEGER, PARAMETER :: li_Y = 2
        INTEGER, PARAMETER :: li_ksi = 1
        INTEGER, PARAMETER :: li_eta = 2
        INTEGER :: li_pt

!        print *, 'apr_points(0,1,1:ai_npts)=', apr_points(0,1,1:ai_npts)
!        print *, 'apr_points(1,1,1:ai_npts)=', apr_points(1,1,1:ai_npts)
!        print *, 'apr_points(2,1,1:ai_npts)=', apr_points(2,1,1:ai_npts)

!        apr_result (1:ai_npts) = apr_points(1,li_i,1:ai_npts) * apr_points(2,li_j,1:ai_npts) &
!        - apr_points(2,li_i,1:ai_npts) * apr_points(1,li_j,1:ai_npts)

        ! x_eta * y_theta - x_theta * y_eta
        apr_Jacobian(1:ai_npts) = 0.0_wp    &
        + apr_points(li_ksi,li_X,1:ai_npts) * apr_points(li_eta,li_Y,1:ai_npts) &
        - apr_points(li_eta,li_X,1:ai_npts) * apr_points(li_ksi,li_Y,1:ai_npts)
!        print *, 'jacobian-2D ', apr_jacobian

        apr_invJacobian (1,1,1:ai_npts) =   apr_points(li_eta,li_Y,1:ai_npts)
        apr_invJacobian (1,2,1:ai_npts) = - apr_points(li_ksi,li_Y,1:ai_npts)
        apr_invJacobian (2,1,1:ai_npts) = - apr_points(li_eta,li_X,1:ai_npts)
        apr_invJacobian (2,2,1:ai_npts) =   apr_points(li_ksi,li_X,1:ai_npts)

        apr_invJacobian_nn (:,:,:) = apr_invJacobian (:,:,:)
        DO li_pt = 1, ai_npts
            apr_invJacobian (:,:,li_pt) = apr_invJacobian (:,:,li_pt) / apr_Jacobian(li_pt)
        END DO

    end subroutine determinant_2D
    !---------------------------------------------------------------
    subroutine determinant_3D( apr_points, ai_npts, apr_Jacobian, apr_invJacobian, apr_invJacobian_nn )
        !> evaluate xi_1 xj_2 - xi_2 * xj_1
        IMPLICIT NONE
        REAL(wp), DIMENSION(0:,:,:) :: apr_points
        INTEGER :: ai_npts
        REAL(wp), DIMENSION(:), intent(inout) :: apr_Jacobian
        REAL(wp), DIMENSION(:,:,:), intent(inout) :: apr_invJacobian
        REAL(wp), DIMENSION(:,:,:), intent(inout) :: apr_invJacobian_nn
        ! LOCAL
        INTEGER, PARAMETER :: li_X = 1
        INTEGER, PARAMETER :: li_Y = 2
        INTEGER, PARAMETER :: li_Z = 3
        INTEGER, PARAMETER :: li_ksi = 1
        INTEGER, PARAMETER :: li_eta = 2
        INTEGER, PARAMETER :: li_theta = 3
        INTEGER :: li_pt

!        print *, 'apr_points(0,1,1:ai_npts)=', apr_points(0,1,1:ai_npts)
!        print *, 'apr_points(1,1,1:ai_npts)=', apr_points(1,1,1:ai_npts)
!        print *, 'apr_points(2,1,1:ai_npts)=', apr_points(2,1,1:ai_npts)

        apr_Jacobian(1:ai_npts) = 0.0_wp    &
        ! + x_ksi * (y_eta * z_theta - y_theta * z_eta)
        + apr_points(li_ksi,li_X,1:ai_npts)   &
        * ( &
        apr_points(li_eta,li_Y,1:ai_npts) * apr_points(li_theta,li_Z,1:ai_npts) &
        - apr_points(li_theta,li_Y,1:ai_npts) * apr_points(li_eta,li_Z,1:ai_npts) &
        )   &
        ! - y_ksi * (x_eta * z_theta - x_theta * z_eta)
        - apr_points(li_ksi,li_Y,1:ai_npts) &
        * ( &
        apr_points(li_eta,li_X,1:ai_npts) * apr_points(li_theta,li_Z,1:ai_npts) &
        - apr_points(li_theta,li_X,1:ai_npts) * apr_points(li_eta,li_Z,1:ai_npts) &
        )   &
        ! + z_ksi * (x_eta * y_theta - x_theta * y_eta)
        + apr_points(li_ksi,li_Z,1:ai_npts) &
        * ( &
        apr_points(li_eta,li_X,1:ai_npts) * apr_points(li_theta,li_Y,1:ai_npts) &
        - apr_points(li_theta,li_X,1:ai_npts) * apr_points(li_eta,li_Y,1:ai_npts) &
        )

        ! (y_eta * z_theta - y_theta * z_eta)
        apr_invJacobian (li_ksi,li_X,1:ai_npts) =   &
              apr_points(li_eta,li_Y,1:ai_npts) * apr_points(li_theta,li_Z,1:ai_npts) &
            - apr_points(li_theta,li_Y,1:ai_npts) * apr_points(li_eta,li_Z,1:ai_npts)
        ! -(y_ksi * z_theta - y_theta * z_ksi)
        apr_invJacobian (li_eta,li_X,1:ai_npts) =   &
            - apr_points(li_ksi,li_Y,1:ai_npts) * apr_points(li_theta,li_Z,1:ai_npts) &
            + apr_points(li_theta,li_Y,1:ai_npts) * apr_points(li_ksi,li_Z,1:ai_npts)
        ! (y_ksi * z_eta - y_eta * z_ksi)
        apr_invJacobian (li_theta,li_X,1:ai_npts) =   &
              apr_points(li_ksi,li_Y,1:ai_npts) * apr_points(li_eta,li_Z,1:ai_npts) &
            - apr_points(li_eta,li_Y,1:ai_npts) * apr_points(li_ksi,li_Z,1:ai_npts)

        ! -(x_eta * z_theta - x_theta * z_eta)
        apr_invJacobian (li_ksi,li_Y,1:ai_npts) =   &
            - apr_points(li_eta,li_X,1:ai_npts) * apr_points(li_theta,li_Z,1:ai_npts) &
            + apr_points(li_theta,li_X,1:ai_npts) * apr_points(li_eta,li_Z,1:ai_npts)
        ! (x_ksi * z_theta - x_theta * z_ksi)
        apr_invJacobian (li_eta,li_Y,1:ai_npts) =   &
            apr_points(li_ksi,li_X,1:ai_npts) * apr_points(li_theta,li_Z,1:ai_npts) &
            - apr_points(li_theta,li_X,1:ai_npts) * apr_points(li_ksi,li_Z,1:ai_npts)
        ! - (x_ksi * z_eta - x_eta * z_ksi)
        apr_invJacobian (li_theta,li_Y,1:ai_npts) =   &
            - apr_points(li_ksi,li_X,1:ai_npts) * apr_points(li_eta,li_Z,1:ai_npts) &
            + apr_points(li_eta,li_X,1:ai_npts) * apr_points(li_ksi,li_Z,1:ai_npts)

        ! (x_eta * y_theta - x_theta * y_eta)
        apr_invJacobian (li_ksi,li_Z,1:ai_npts) =   &
            apr_points(li_eta,li_X,1:ai_npts) * apr_points(li_theta,li_Y,1:ai_npts) &
            - apr_points(li_theta,li_X,1:ai_npts) * apr_points(li_eta,li_Y,1:ai_npts)
        ! -(x_ksi * y_theta - x_theta * y_ksi)
        apr_invJacobian (li_eta,li_Z,1:ai_npts) =   &
            apr_points(li_ksi,li_X,1:ai_npts) * apr_points(li_theta,li_Y,1:ai_npts) &
            - apr_points(li_theta,li_X,1:ai_npts) * apr_points(li_ksi,li_Y,1:ai_npts)
        ! (x_ksi * y_eta - x_eta * y_ksi)
        apr_invJacobian (li_theta,li_Z,1:ai_npts) =   &
            apr_points(li_ksi,li_X,1:ai_npts) * apr_points(li_eta,li_Y,1:ai_npts) &
            - apr_points(li_eta,li_X,1:ai_npts) * apr_points(li_ksi,li_Y,1:ai_npts)

        apr_invJacobian_nn (:,:,:) = apr_invJacobian (:,:,:)
        
        DO li_pt = 1, ai_npts
            apr_invJacobian (:,:,li_pt) = apr_invJacobian (:,:,li_pt) / apr_Jacobian(li_pt)
        END DO

    end subroutine determinant_3D
    !---------------------------------------------------------------
    subroutine Assembly_space_jacobian( self, ao_FEM, ai_id, ai_npts )
        IMPLICIT NONE
        TYPE(ASSEMBLY) :: self
        TYPE(FEM) :: ao_FEM
        INTEGER :: ai_id
        INTEGER :: ai_npts
        ! LOCAL
        INTEGER :: li_grids
        TYPE(METRIC_INFO), POINTER :: lp_info

        CALL printlog("Assembly_space_jacobian : Start", ai_dtllevel = mi_dtllevel_base  + 1)

        lp_info => self % opo_info_sp (ai_id)
        li_grids   = ao_FEM % opi_InfoSpace(ai_id, INFOSPACE_GRIDS )
!        print *, 'ai_id, li_grids, ao_FEM % opi_dim(li_grids) =', ai_id, li_grids,ao_FEM % opi_dim(li_grids)
!        print *, 'dpoints =', lp_info % opr_points(1,1,1:ai_npts)
        select case ( ao_FEM % opi_dim(li_grids) )

            case ( DIMENSION_1D )
!                lp_info % opr_jacobians(1:ai_npts)= lp_info % opr_points(1,1,1:ai_npts)
!                lp_info % opr_invJacobian(1,1,1:ai_npts)= 1.0_wp / lp_info % opr_points(1,1,1:ai_npts)
!                lp_info % opr_invJacobian_nn(1,1,1:ai_npts)= 1.0_wp
                CALL determinant_1D( lp_info % opr_points, ai_npts  &
                , lp_info % opr_jacobians   &
                , lp_info % opr_invJacobian &
                , lp_info % opr_invJacobian_nn )

            case ( DIMENSION_2D )
                CALL determinant_2D( lp_info % opr_points, ai_npts  &
                , lp_info % opr_jacobians   &
                , lp_info % opr_invJacobian &
                , lp_info % opr_invJacobian_nn )

            case ( DIMENSION_3D )
                CALL determinant_3D( lp_info % opr_points, ai_npts  &
                , lp_info % opr_jacobians   &
                , lp_info % opr_invJacobian &
                , lp_info % opr_invJacobian_nn )

            case Default
                call printlog("Dimension Not Yet implemented", ai_dtllevel = mi_dtllevel_base +  0)

        end select

#ifdef _DEBUG
call concatmsg("lp_info % opr_jacobians = ", ai_dtllevel = mi_dtllevel_base + 1)
call concatmsg(lp_info % opr_jacobians, ai_dtllevel = mi_dtllevel_base + 1)
call printmsg(ai_dtllevel = mi_dtllevel_base + 1)
#endif
        
        CALL printlog("Assembly_space_jacobian : End", ai_dtllevel = mi_dtllevel_base  + 1)

    end subroutine Assembly_space_jacobian
    !---------------------------------------------------------------
    subroutine Assembly_mapping_jacobian( self, ao_FEM, ai_id, ai_npts )
        IMPLICIT NONE
        TYPE(ASSEMBLY) :: self
        TYPE(FEM) :: ao_FEM
        INTEGER :: ai_id
        INTEGER :: ai_npts
        ! LOCAL
        INTEGER :: li_grids
        INTEGER :: li_space
        TYPE(METRIC_INFO), POINTER :: lp_info

        CALL printlog("Assembly_mapping_jacobian : Start", ai_dtllevel = mi_dtllevel_base  + 1)

        lp_info => self % opo_info_mp (ai_id)
        li_space   = ao_FEM % opi_InfoMapping (ai_id, INFOMAPPING_SPACE )
        li_grids   = ao_FEM % opi_InfoSpace (li_space, INFOSPACE_GRIDS )
        select case ( ao_FEM % opi_dim(li_grids) )

            case ( DIMENSION_1D )
!                lp_info % opr_jacobians(1:ai_npts)= lp_info % opr_points(1,1,1:ai_npts)
!                lp_info % opr_invJacobian(1,1,1:ai_npts)= 1.0_wp / lp_info % opr_points(1,1,1:ai_npts)
!                lp_info % opr_invJacobian_nn(1,1,1:ai_npts)= 1.0_wp
                CALL determinant_1D( lp_info % opr_points, ai_npts  &
                , lp_info % opr_jacobians   &
                , lp_info % opr_invJacobian &
                , lp_info % opr_invJacobian_nn )

            case ( DIMENSION_2D )
                CALL determinant_2D( lp_info % opr_points, ai_npts  &
                , lp_info % opr_jacobians   &
                , lp_info % opr_invJacobian &
                , lp_info % opr_invJacobian_nn )

            case ( DIMENSION_3D )
                CALL determinant_3D( lp_info % opr_points, ai_npts  &
                , lp_info % opr_jacobians   &
                , lp_info % opr_invJacobian &
                , lp_info % opr_invJacobian_nn )

            case Default
                call printlog("Dimension Not Yet implemented", ai_dtllevel = mi_dtllevel_base +  0)

        end select

!        print *, 'lp_info % opr_jacobians=', lp_info % opr_jacobians

        CALL printlog("Assembly_mapping_jacobian : End", ai_dtllevel = mi_dtllevel_base  + 1)

    end subroutine Assembly_mapping_jacobian
    !---------------------------------------------------------------
    subroutine Assembly_metric_jacobian( self, ao_FEM, ai_grids, ai_id, ai_npts )
        IMPLICIT NONE
        TYPE(ASSEMBLY) :: self
        TYPE(FEM) :: ao_FEM
        INTEGER :: ai_grids
        INTEGER :: ai_id
        INTEGER :: ai_npts
        ! LOCAL
        INTEGER :: li_space
        TYPE(METRIC_INFO), POINTER :: lp_info

        CALL printlog("Assembly_metric_jacobian : Start", ai_dtllevel = mi_dtllevel_base  + 1)

        lp_info => self % opo_info_gr (ai_id)
        select case ( ao_FEM % opi_dim(ai_grids) )

            case ( DIMENSION_1D )
!                lp_info % opr_jacobians(1:ai_npts)= lp_info % opr_points(1,1,1:ai_npts)
!                lp_info % opr_invJacobian(1,1,1:ai_npts)= 1.0_wp / lp_info % opr_points(1,1,1:ai_npts)
!                lp_info % opr_invJacobian_nn(1,1,1:ai_npts)= 1.0_wp
                CALL determinant_1D( lp_info % opr_points, ai_npts  &
                , lp_info % opr_jacobians   &
                , lp_info % opr_invJacobian &
                , lp_info % opr_invJacobian_nn )

            case ( DIMENSION_2D )
                CALL determinant_2D( lp_info % opr_points, ai_npts  &
                , lp_info % opr_jacobians   &
                , lp_info % opr_invJacobian &
                , lp_info % opr_invJacobian_nn )

            case ( DIMENSION_3D )
                CALL determinant_3D( lp_info % opr_points, ai_npts  &
                , lp_info % opr_jacobians   &
                , lp_info % opr_invJacobian &
                , lp_info % opr_invJacobian_nn )

            case Default
                call printlog("Dimension Not Yet implemented", ai_dtllevel = mi_dtllevel_base +  0)

        end select

        CALL printlog("Assembly_metric_jacobian : End", ai_dtllevel = mi_dtllevel_base  + 1)

    end subroutine Assembly_metric_jacobian
    !---------------------------------------------------------------
    subroutine copy_info_metric ( ao_in, ao_out )
        !> remember that points may comptain second derivatives in the case of bilaplacian operator
        !> we need to cap the shape with output one
        IMPLICIT NONE
        TYPE(METRIC_INFO) :: ao_in
        TYPE(METRIC_INFO) :: ao_out
        ! LOCAL
        INTEGER  :: li_nderiv
        INTEGER  :: li_Rd
        INTEGER  :: li_npts

        li_nderiv = MIN(SIZE(ao_in % opr_points, 1) - 1, SIZE(ao_out % opr_points, 1) - 1)
        li_Rd = MIN(SIZE(ao_in % opr_points, 2), SIZE(ao_out % opr_points, 2))
        li_npts = MIN(SIZE(ao_in % opr_points, 3), SIZE(ao_out % opr_points, 3))

        ao_out % opr_points (0:li_nderiv, 1:li_Rd, 1:li_npts) = ao_in % opr_points (0:li_nderiv, 1:li_Rd, 1:li_npts)
        ao_out % opr_jacobians = ao_in % opr_jacobians
        ao_out % opr_invJacobian = ao_in % opr_invJacobian
        ao_out % opr_invJacobian_nn = ao_in % opr_invJacobian_nn

    end subroutine copy_info_metric
    !---------------------------------------------------------------
    subroutine update_info_metric( self, ao_FEM, ai_space, ai_patch_id, ai_elt, ai_npts )
        IMPLICIT NONE
        TYPE(ASSEMBLY) :: self
        TYPE(FEM) :: ao_FEM
        INTEGER :: ai_space
        INTEGER :: ai_patch_id
        INTEGER :: ai_elt
        INTEGER :: ai_npts
        ! LOCAL
        INTEGER :: li_grids
        INTEGER :: li_extmap
        INTEGER :: li_usemetric

        CALL printlog("update_info_metric : Start", ai_dtllevel = mi_dtllevel_base  + 1)
!        print*, 'update_info_metric : Start'

!        IF (  self % ol_assembly_basis .OR. self % ol_assembly_points) THEN

            li_grids       = ao_FEM % opi_InfoSpace (ai_space, INFOSPACE_GRIDS )
            li_extmap      = ao_FEM % opi_InfoSpace (ai_space, INFOSPACE_EXTMAPPING )
            li_usemetric  = ao_FEM % opi_InfoGrids (li_grids, INFOGRIDS_USEMETRIC )

!            print *, 'li_usemetric, li_extmap = ', li_usemetric, li_extmap

            IF (li_usemetric == 1) THEN
                IF (li_extmap == 1) THEN
                    CALL update_info_mapping( self, ao_FEM, ai_space, ai_patch_id, ai_elt, ai_npts )
                ELSE
                    CALL update_points_metric( self, ao_FEM, ai_space, ai_patch_id, ai_elt, ai_npts )
                END IF
            END IF

!        END IF
!        print*, 'update_info_metric : End'

        CALL printlog("update_info_metric : End", ai_dtllevel = mi_dtllevel_base  + 1)

    end subroutine update_info_metric
    !---------------------------------------------------------------
    subroutine update_points_metric( self, ao_FEM, ai_space, ai_patch_id, ai_elt, ai_npts )
        IMPLICIT NONE
        TYPE(ASSEMBLY) :: self
        TYPE(FEM) :: ao_FEM
        INTEGER :: ai_space
        INTEGER :: ai_patch_id
        INTEGER :: ai_elt
        INTEGER :: ai_npts
        ! LOCAL
        INTEGER :: li_grids
        INTEGER :: li_Rd
        INTEGER :: li_d
        INTEGER :: li_deriv
        INTEGER :: li_nderiv
        INTEGER :: li_usemetric
        INTEGER :: li_metric_id
        TYPE(METRIC_INFO), POINTER :: lp_info
        TYPE(METRIC), POINTER :: lp_metric

        CALL printlog("update_points_metric : Start", ai_dtllevel = mi_dtllevel_base  + 1)
!        print*, 'update_points_metric: Start'


        li_grids       = ao_FEM % opi_InfoSpace (ai_space, INFOSPACE_GRIDS )
        li_metric_id  = ao_FEM % opi_InfoGrids (li_grids, INFOGRIDS_METRIC_ID )
        li_Rd           = ao_FEM % opi_InfoGrids (li_grids, INFOGRIDS_RD )
        li_nderiv      = self % oi_mapping_nderiv        
        
#ifdef _DEBUG
        call concatmsg("Treatment of the metric of id : ", ai_dtllevel = 0)
        call concatmsg(li_metric_id, ai_dtllevel = 0)
        call printmsg(ai_dtllevel = 0)
#endif

        lp_info => self % opo_info_gr (li_grids)
        lp_metric => ao_FEM % opo_metrics (li_metric_id)

!        print *, ao_FEM % opo_metrics (li_metric_id) % opr_points(ai_patch_id, 1, 1, 0, 1:li_Rd)
!        print *, ao_FEM % opo_metrics (li_metric_id) % opr_points(ai_patch_id, 1, 2, 0, 1:li_Rd)
!        print *, ao_FEM % opo_metrics (li_metric_id) % opr_points(ai_patch_id, 1, 3, 0, 1:li_Rd)
!        print *, ao_FEM % opo_metrics (li_metric_id) % opr_points(ai_patch_id, 1, 4, 0, 1:li_Rd)

        DO li_deriv = 0, li_nderiv
            DO li_d = 1, li_Rd
                lp_info % opr_points (li_deriv,li_d,1:ai_npts) =  &
                lp_metric % opr_points(ai_patch_id, ai_elt,1:ai_npts, li_deriv, li_d)
            END DO

!            if (ai_elt==1) then
!                print *, '====',lp_info % opr_points (0,:,1), '===='
!                print *, '====',lp_info % opr_points (0,:,2), '===='
!                print *, '====',lp_info % opr_points (0,:,3), '===='
!                print *, '====',lp_info % opr_points (0,:,4), '===='
!            end if
        END DO

        !> evaluation of the jacobian on all points of the grid, inside the current element
        !> \todo
        CALL Assembly_metric_jacobian( self, ao_FEM, li_grids, ai_patch_id, ai_npts )

#ifdef _DEBUG
        call printlog("done.", ai_dtllevel = 0)
#endif
!        print*, 'update_points_metric: End'

        CALL printlog("update_points_metric : End", ai_dtllevel = mi_dtllevel_base  + 1)

    end subroutine update_points_metric
    !---------------------------------------------------------------
    subroutine update_info_mapping( self, ao_FEM, ai_space, ai_patch_id, ai_elt, ai_npts )
        IMPLICIT NONE
        TYPE(ASSEMBLY) :: self
        TYPE(FEM) :: ao_FEM
        INTEGER :: ai_space
        INTEGER :: ai_patch_id
        INTEGER :: ai_elt
        INTEGER :: ai_npts
        ! LOCAL
        INTEGER :: li_grids
        INTEGER :: li_REALelt
        INTEGER :: li_map
        INTEGER :: li_extmap
        INTEGER :: li_usemetric
        INTEGER :: li_metric_id
        INTEGER :: li_ptw_evaluation
        TYPE(GRID_DATA), POINTER :: lp_grid
        TYPE(BBOX), POINTER :: lp_bbox
        TYPE(MAPPING), POINTER :: lp_mapping
        TYPE(SPACE), POINTER :: lp_space
        TYPE(METRIC_INFO), POINTER :: lp_info

        CALL printlog("update_info_mapping : Start", ai_dtllevel = mi_dtllevel_base  + 1)

        li_ptw_evaluation = 1

        li_grids       = ao_FEM % opi_InfoSpace (ai_space, INFOSPACE_GRIDS )
        li_map          = ao_FEM % opi_InfoSpace (ai_space, INFOSPACE_MAPPING )
        li_metric_id  = ao_FEM % opi_InfoGrids (li_grids, INFOGRIDS_METRIC_ID )

#ifdef _DEBUG
        call concatmsg("Treatment of the mapping of id : ", ai_dtllevel = 0)
        call concatmsg(li_map, ai_dtllevel = 0)
        call printmsg(ai_dtllevel = 0)
#endif

        lp_info => self % opo_info_gr (li_grids)
        lp_bbox => self % opo_bbox_mp ( li_map )
        lp_mapping => ao_FEM % opo_mappings ( li_map )
        lp_grid => ao_FEM % opo_grids(li_grids) % opo_grid ( ai_patch_id )

        li_REALelt = ao_FEM % opo_spaces (ai_space) % oo_con % opi_REAL_elts (ai_patch_id,ai_elt)

        !> \todo a changer: ai_patch_id+1 en ai_patch_id, ds la lib grids
        CALL assembly_logical_basis(lp_grid, lp_mapping % oo_mapping &
        , self % oi_mapping_nderiv, ai_patch_id+1, ai_elt, li_REALelt &
        , self % opo_bbox_mp ( li_map ) &
!        , ai_nen = lp_mapping % oi_maxnen &
        , lp_mapping % oi_maxnen &
!        , ai_ptw_evaluation = li_ptw_evaluation)
        , li_ptw_evaluation)

        CALL assembly_points_elements(lp_grid, lp_mapping % oo_mapping    &
        , self % opo_bbox_mp ( li_map )   &
        , ai_patch_id+1, ai_elt, li_REALelt, ai_npts  &
        , self % opo_info_mp (li_map) % opr_points &
!        , ai_ptw_evaluation = li_ptw_evaluation)
        , li_ptw_evaluation)

!#ifdef _MYDEBUG
!#endif

        !> evaluation of the jacobian on all points of the grid, inside the current element
        !> \todo
        CALL Assembly_mapping_jacobian( self, ao_FEM, li_map, ai_npts )

        CALL copy_info_metric( self % opo_info_mp (li_map), self % opo_info_gr (li_metric_id))

#ifdef _DEBUG
        call printlog("done.", ai_dtllevel = 0)
#endif

        CALL printlog("update_info_mapping : End", ai_dtllevel = mi_dtllevel_base  + 1)

    end subroutine update_info_mapping
    !---------------------------------------------------------------
    subroutine update_info_space( self, ao_FEM, ai_space, ai_patch_id, ai_elt, ai_npts )
        IMPLICIT NONE
        TYPE(ASSEMBLY) :: self
        TYPE(FEM) :: ao_FEM
        INTEGER :: ai_space
        INTEGER :: ai_patch_id
        INTEGER :: ai_elt
        INTEGER :: ai_npts
        ! LOCAL
        INTEGER :: li_grids
        INTEGER :: li_REALelt
        INTEGER :: li_map
        INTEGER :: li_extmap
        INTEGER :: li_usemetric
        INTEGER :: li_metric_id
        INTEGER :: li_composed
        INTEGER :: li_composed_nsp
        INTEGER :: li_position_sp
        INTEGER :: li_single_sp
        INTEGER :: li_single_grids_id
        TYPE(GRID_DATA), POINTER :: lp_grid
        TYPE(GRID_DATA), POINTER :: lp_single_grid
        TYPE(BBOX), POINTER :: lp_bbox
        TYPE(SPACE), POINTER :: lp_space
        TYPE(METRIC_INFO), POINTER :: lp_info

        CALL printlog("update_info_space : Start", ai_dtllevel = mi_dtllevel_base  + 1)
!        print*, "update_info_space : Start"

        li_grids   = ao_FEM % opi_InfoSpace (ai_space, INFOSPACE_GRIDS )
        li_extmap      = ao_FEM % opi_InfoSpace (ai_space, INFOSPACE_EXTMAPPING )
        li_usemetric  = ao_FEM % opi_InfoGrids (li_grids, INFOGRIDS_USEMETRIC )
!        li_metric_id  = ao_FEM % opi_InfoGrids (li_grids, INFOGRIDS_METRIC_ID )
        
#ifdef _DEBUG
        call concatmsg("Treatment of the space of id : ", ai_dtllevel = 0)
        call concatmsg(ai_space, ai_dtllevel = 0)
        call printmsg(ai_dtllevel = 0)
#endif

        lp_bbox => self % opo_bbox_sp ( ai_space )
        lp_space => ao_FEM % opo_spaces ( ai_space )
        lp_grid => ao_FEM % opo_grids(li_grids) % opo_grid ( ai_patch_id )

        li_REALelt = ao_FEM % opo_spaces (ai_space) % oo_con % opi_REAL_elts (ai_patch_id,ai_elt)

        !> \todo a changer: ai_patch_id+1 en ai_patch_id, ds la lib grids
        CALL assembly_logical_basis(lp_grid, lp_space % oo_mapping &
        , self % oi_matrix_nderiv, ai_patch_id+1, ai_elt, li_REALelt  &
        , self % opo_bbox_sp ( ai_space ) &
!        , ai_nen = lp_space % oi_maxnen &
        , lp_space % oi_maxnen &
!        , ai_ptw_evaluation = 0)
        , 0)

        !> \todo a traiter le cas des pullback
        !> \todo
        !> \todo a changer: ai_patch_id+1 en ai_patch_id, ds la lib grids

        IF (li_usemetric == 0) THEN
            CALL assembly_points_elements(lp_grid, lp_space % oo_mapping    &
            , self % opo_bbox_sp ( ai_space )   &
            , ai_patch_id+1, ai_elt, li_REALelt, ai_npts  &
            , self % opo_info_sp (ai_space) % opr_points &
!            , ai_ptw_evaluation = 0)
            , 0)
            
            !> evaluation of the jacobian on all points of the grid, inside the current element
            CALL Assembly_space_jacobian( self, ao_FEM, ai_space, ai_npts )

            CALL copy_info_metric( self % opo_info_sp (ai_space), self % opo_info_gr (li_grids))

!            print *, '==========================================================='
!            print *, 'x  =', self % opo_info_sp (ai_space) % opr_points(0, 1, 1:ai_npts)
!            print *, 'y  =', self % opo_info_sp (ai_space) % opr_points(0, 2, 1:ai_npts)
!!            print *, 'z  =', self % opo_info_sp (ai_space) % opr_points(0, 3, 1:ai_npts)
!            print *, '==========================================================='
!            print *, '==========================================================='
!            print *, 'x1 =', self % opo_info_sp (ai_space) % opr_points(1, 1, 1:ai_npts)
!            print *, 'y1 =', self % opo_info_sp (ai_space) % opr_points(1, 2, 1:ai_npts)
!!            print *, 'z1 =', self % opo_info_sp (ai_space) % opr_points(1, 3, 1:ai_npts)
!            print *, '==========================================================='
!            print *, '==========================================================='
!            print *, 'x2 =', self % opo_info_sp (ai_space) % opr_points(2, 1, 1:ai_npts)
!            print *, 'y2 =', self % opo_info_sp (ai_space) % opr_points(2, 2, 1:ai_npts)
!!            print *, 'z2 =', self % opo_info_sp (ai_space) % opr_points(2, 3, 1:ai_npts)
!            print *, '==========================================================='
!            print *, '==========================================================='
!            print *, 'jacobian =', self % opo_info_sp (ai_space) % opr_jacobians(1:ai_npts)
!            print *, '==========================================================='

        END IF

#ifdef _DEBUG
        call printlog("done.", ai_dtllevel = 0)
#endif
!        print*, "update_info_space : End"
        CALL printlog("update_info_space : End", ai_dtllevel = mi_dtllevel_base  + 1)

    end subroutine update_info_space
    !---------------------------------------------------------------
    subroutine assembly_physical_basis( self, ao_FEM, ai_space, ai_patch, ai_elt )
        IMPLICIT NONE
        TYPE(ASSEMBLY) :: self
        TYPE(FEM) :: ao_FEM
        INTEGER, INTENT(IN)  :: ai_space
        INTEGER, INTENT(IN)  :: ai_patch
        INTEGER, INTENT(IN)  :: ai_elt
        ! LOCAL VARIABLES
        INTEGER  :: li_dim
        INTEGER  :: li_grids
        INTEGER  :: li_npts
        TYPE(GRID_DATA), POINTER :: lp_grid

        CALL printlog("assembly_physical_basis: Begin", ai_dtllevel = 1)

        ! ...
        li_grids        = ao_FEM % opi_InfoSpace ( ai_space, INFOSPACE_GRIDS )
        li_dim          = ao_FEM % opi_dim(li_grids)
        lp_grid         => ao_FEM % opo_grids(li_grids) % opo_grid ( ai_patch )
        li_npts         = lp_grid % opo_elts ( ai_elt ) % oi_npts
        ! ...

        IF (li_dim==1) THEN
                CALL assembly_physical_basis_1D( self, ao_FEM, ai_space, ai_patch, ai_elt, li_npts)
        END IF

        IF (li_dim==2) THEN
                CALL assembly_physical_basis_2D( self, ao_FEM, ai_space, ai_patch, ai_elt, li_npts)
        END IF

        IF (li_dim==3) THEN
                CALL assembly_physical_basis_3D( self, ao_FEM, ai_space, ai_patch, ai_elt, li_npts)
        END IF

        CALL printlog("assembly_physical_basis: End", ai_dtllevel = 1)

    end subroutine assembly_physical_basis    
    !---------------------------------------------------------------
    subroutine assembly_physical_basis_1D( self, ao_FEM, ai_space, ai_patch, ai_elt, ai_npts)
        IMPLICIT NONE
        TYPE(ASSEMBLY) :: self
        TYPE(FEM) :: ao_FEM
        INTEGER, INTENT(IN)  :: ai_space
        INTEGER, INTENT(IN)  :: ai_patch
        INTEGER, INTENT(IN)  :: ai_elt
        INTEGER, INTENT(IN)  :: ai_npts
        ! LOCAL VARIABLES
        INTEGER, PARAMETER :: XI        = 1
        INTEGER, PARAMETER :: D_XI      = 1 
        INTEGER, PARAMETER :: D_XIXI    = 2
        INTEGER  :: li_dim
        INTEGER  :: li_elt
        INTEGER  :: li_grids
        INTEGER  :: li_npts
        INTEGER  :: li_pt
        INTEGER  :: li_deriv
        INTEGER  :: li_nen
        INTEGER  :: li_b
        INTEGER  :: li_nderiv
        TYPE(GRID_DATA), POINTER :: lp_grid
        TYPE(BBOX), POINTER :: lp_bbox
        TYPE(SPACE), POINTER :: lp_space
        REAL(wp), DIMENSION(:,:,:), POINTER :: lpr_invJacobian
        REAL(wp), DIMENSION(:,:), POINTER :: lpr_B
        REAL(wp), DIMENSION(:,:,:), POINTER :: lpr_gradB
        REAL(wp), DIMENSION(:,:,:), POINTER :: lpr_curlB
        REAL(wp), DIMENSION (:), POINTER :: lpr_XI_X
        REAL(wp), DIMENSION (:), POINTER :: lpr_XI_XX
        REAL(wp), DIMENSION (:,:), POINTER :: lpr_Phi_XI
        REAL(wp), DIMENSION (:,:), POINTER :: lpr_Phi_XIXI        

        CALL printlog("assembly_physical_basis_1D: Begin", ai_dtllevel = 1)
        
        li_grids = ao_FEM % opi_InfoSpace ( ai_space, INFOSPACE_GRIDS )

        lp_bbox         => self % opo_bbox_sp ( ai_space )
        lp_space        => ao_FEM % opo_spaces ( ai_space )
        lp_grid         => ao_FEM % opo_grids(li_grids) % opo_grid ( ai_patch )
        lpr_invJacobian => self % opo_info_gr (li_grids) % opr_invJacobian
        lpr_B           => self % opo_pBasis (ai_space) % opr_B
        lpr_gradB       => self % opo_pBasis (ai_space) % opr_gradB
        lpr_curlB       => self % opo_pBasis (ai_space) % opr_curlB

        li_dim          = ao_FEM % opi_dim(li_grids)
        li_npts         = lp_grid % opo_elts ( ai_elt ) % oi_npts
        li_nen          = lp_space % oo_con % opi_nen(ai_patch+1)

        li_nderiv       = MAX(self % oi_matrix_nderiv, self % oi_mapping_nderiv, 2)

        lpr_XI_X        => self % opo_info_gr (li_grids) % opr_points (D_XI , XI, 1:li_npts)
        lpr_Phi_XI      => lp_bbox % opr_dBasis (D_XI  , 0:li_nen-1, 1:li_npts)

        ! ...
        ! assembling B
        ! ...
        self % opo_pBasis (ai_space) % opr_B (1:li_nen, 1:li_npts) = lp_bbox % opr_dBasis (0,0:li_nen-1,1:li_npts)
        ! ...

        ! ...
        ! assembling grad B = ( dB/dx, dB/dy, dB/dz ) 
        ! ...
        do li_b =1, li_nen
        do li_pt=1, li_npts
            lpr_gradB (1:li_dim, li_b, li_pt) = MATMUL(lpr_invJacobian(1:li_dim, 1:li_dim, li_pt) &
            , lp_bbox % opr_dBasis (1:li_dim, li_b-1, li_pt))
        end do
        end do
        ! ...

        ! ...
        ! assembling curl B = 
        ! ...

        ! ...

        IF (li_nderiv < 2) THEN
            return
        END IF

        ! ...
        ! assembling Hessian B = ( dB/dxx, dB/dyy, dB/dzz, dB/dxy, dB/dyz, dB/dzx )
        ! ...
        lpr_XI_XX       => self % opo_info_gr (li_grids) % opr_points (D_XIXI, XI, 1:li_npts)
        lpr_Phi_XIXI    => lp_bbox % opr_dBasis (D_XIXI, 0:li_nen-1, 1:li_npts)
        
        ! phi_{x,x} = xi_{x,x} phi_xi + (xi_x)^2 phi_{xi,xi}
        DO li_pt = 1, li_npts
                self % opo_pBasis (ai_space) % opr_HessianB (1, 1:li_nen, li_pt) = &
                lpr_XI_XX(li_pt) * lpr_Phi_XI (1:li_nen, li_pt) &
                + lpr_XI_X(li_pt)**2 * lpr_Phi_XIXI (1:li_nen, li_pt)               
        END DO
        ! ...

        CALL printlog("assembly_physical_basis_1D: End", ai_dtllevel = 1)

    end subroutine assembly_physical_basis_1D    
    !---------------------------------------------------------------
    subroutine assembly_physical_basis_2D( self, ao_FEM, ai_space, ai_patch, ai_elt, ai_npts)
        IMPLICIT NONE
        TYPE(ASSEMBLY) :: self
        TYPE(FEM) :: ao_FEM
        INTEGER, INTENT(IN)  :: ai_space
        INTEGER, INTENT(IN)  :: ai_patch
        INTEGER, INTENT(IN)  :: ai_elt
        INTEGER, INTENT(IN)  :: ai_npts
        ! LOCAL VARIABLES
        INTEGER, PARAMETER :: li_MAXDIM = 3
        INTEGER, PARAMETER :: li_MAXDIM_invH = 3
        INTEGER, PARAMETER :: li_NPARAM = 3
        INTEGER, PARAMETER :: X         = 1
        INTEGER, PARAMETER :: Y         = 2
        INTEGER, PARAMETER :: D_XI      = 1
        INTEGER, PARAMETER :: D_ETA     = 2 
        INTEGER, PARAMETER :: D_XIXI    = 3
        INTEGER, PARAMETER :: D_XIETA   = 4
        INTEGER, PARAMETER :: D_ETAETA  = 5 
        INTEGER, PARAMETER :: D_XX      = 1
        INTEGER, PARAMETER :: D_XY      = 2
        INTEGER, PARAMETER :: D_YY      = 3         
        INTEGER  :: li_dim
        INTEGER  :: li_elt
        INTEGER  :: li_nen
        INTEGER  :: li_grids
        INTEGER  :: li_npts
        INTEGER  :: li_pt
        INTEGER  :: li_deriv
        INTEGER  :: li_nenA
        INTEGER  :: li_b
        INTEGER  :: li_d
        INTEGER  :: li_nderiv
        INTEGER  :: li_nH
        TYPE(GRID_DATA), POINTER :: lp_grid
        TYPE(BBOX), POINTER :: lp_bbox
        TYPE(SPACE), POINTER :: lp_space
        REAL(wp), DIMENSION (:,:,:), POINTER :: lpr_invJacobian
        REAL(wp), DIMENSION (:,:), POINTER :: lpr_B
        REAL(wp), DIMENSION (:,:,:), POINTER :: lpr_gradB
        REAL(wp), DIMENSION (:,:,:), POINTER :: lpr_curlB
        REAL(wp), DIMENSION (:,:,:), POINTER :: lpr_HessianB
        REAL(wp), DIMENSION (:,:), POINTER :: lpr_Phi
        REAL(wp), DIMENSION (:,:), POINTER :: lpr_Phi_X
        REAL(wp), DIMENSION (:,:), POINTER :: lpr_Phi_Y
        REAL(wp), DIMENSION (:,:), POINTER :: lpr_Phi_XI
        REAL(wp), DIMENSION (:,:), POINTER :: lpr_Phi_ETA
        REAL(wp), DIMENSION (:,:), POINTER :: lpr_Phi_XIXI
        REAL(wp), DIMENSION (:,:), POINTER :: lpr_Phi_XIETA
        REAL(wp), DIMENSION (:,:), POINTER :: lpr_Phi_ETAETA
        REAL(wp), DIMENSION (:), POINTER :: lpr_X
        REAL(wp), DIMENSION (:), POINTER :: lpr_Y        
        REAL(wp), DIMENSION (:), POINTER :: lpr_X_XI
        REAL(wp), DIMENSION (:), POINTER :: lpr_X_ETA
        REAL(wp), DIMENSION (:), POINTER :: lpr_Y_XI
        REAL(wp), DIMENSION (:), POINTER :: lpr_Y_ETA
        REAL(wp), DIMENSION (:), POINTER :: lpr_X_XIXI
        REAL(wp), DIMENSION (:), POINTER :: lpr_X_XIETA
        REAL(wp), DIMENSION (:), POINTER :: lpr_X_ETAETA
        REAL(wp), DIMENSION (:), POINTER :: lpr_Y_XIXI
        REAL(wp), DIMENSION (:), POINTER :: lpr_Y_XIETA
        REAL(wp), DIMENSION (:), POINTER :: lpr_Y_ETAETA        
        REAL(wp), DIMENSION (ai_npts) :: lpr_XI_X     ! d XI / d x
        REAL(wp), DIMENSION (ai_npts) :: lpr_ETA_X    ! d ETA / d x
        REAL(wp), DIMENSION (ai_npts) :: lpr_XI_Y     ! d XI / d y
        REAL(wp), DIMENSION (ai_npts) :: lpr_ETA_Y    ! d ETA / d y
        REAL(wp), DIMENSION (ai_npts) :: lpr_J        ! jacobian
        REAL(wp), DIMENSION (ai_npts) :: lpr_uJ_X      ! ( d jacobian d x ) / jacobian
        REAL(wp), DIMENSION (ai_npts) :: lpr_uJ_Y      ! ( d jacobian d y ) / jacobian
        REAL(wp), DIMENSION (ai_npts) :: lpr_XI_XX     ! d2 XI / d xx
        REAL(wp), DIMENSION (ai_npts) :: lpr_ETA_XX    ! d2 ETA / d xx
        REAL(wp), DIMENSION (ai_npts) :: lpr_XI_XY     ! d2 XI / d xy
        REAL(wp), DIMENSION (ai_npts) :: lpr_ETA_XY    ! d2 ETA / d xy        
        REAL(wp), DIMENSION (ai_npts) :: lpr_XI_YY     ! d2 XI / d yy
        REAL(wp), DIMENSION (ai_npts) :: lpr_ETA_YY    ! d2 ETA / d yy        
        REAL(wp), DIMENSION (ao_FEM % opo_spaces (ai_space) % oi_maxnen, ai_npts) :: lpr_C1
        REAL(wp), DIMENSION (ao_FEM % opo_spaces (ai_space) % oi_maxnen, ai_npts) :: lpr_C2
        REAL(wp), DIMENSION (ao_FEM % opo_spaces (ai_space) % oi_maxnen, ai_npts) :: lpr_C3

        CALL printlog("assembly_physical_basis_2D: Begin", ai_dtllevel = 1)
!        print*,"assembly_physical_basis_2D: Begin"
        
        li_grids = ao_FEM % opi_InfoSpace ( ai_space, INFOSPACE_GRIDS )

        ! ...
        lp_bbox         => self % opo_bbox_sp ( ai_space )
        lp_space        => ao_FEM % opo_spaces ( ai_space )
        lp_grid         => ao_FEM % opo_grids(li_grids) % opo_grid ( ai_patch )

        li_dim          = ao_FEM % opi_dim(li_grids)
        li_npts         = lp_grid % opo_elts ( ai_elt ) % oi_npts
        li_nen          = lp_space % oo_con % opi_nen(ai_patch+1)

        li_nderiv       = MAX(self % oi_matrix_nderiv, self % oi_mapping_nderiv, 2)
        li_nH           = li_dim * (li_dim + 1) / 2

        lpr_invJacobian => self % opo_info_gr (li_grids) % opr_invJacobian (1:li_dim, 1:li_dim, 1:li_npts)
        lpr_B           => self % opo_pBasis (ai_space) % opr_B (1:li_nen, 1:li_npts)
        lpr_gradB       => self % opo_pBasis (ai_space) % opr_gradB (1:li_dim, 1:li_nen, 1:li_npts)
        lpr_curlB       => self % opo_pBasis (ai_space) % opr_curlB (1:li_dim, 1:li_nen, 1:li_npts)
        lpr_HessianB    => self % opo_pBasis (ai_space) % opr_HessianB (1:li_nH, 1:li_nen, 1:li_npts)
        ! ...

        ! ...
        lpr_X           => self % opo_info_gr (li_grids) % opr_points(0 , X, 1:li_npts)
        lpr_Y           => self % opo_info_gr (li_grids) % opr_points(0 , Y, 1:li_npts)

        lpr_X_XI        => self % opo_info_gr (li_grids) % opr_points(D_XI , X, 1:li_npts)
        lpr_X_ETA       => self % opo_info_gr (li_grids) % opr_points(D_ETA, X, 1:li_npts)
        lpr_Y_XI        => self % opo_info_gr (li_grids) % opr_points(D_XI , Y, 1:li_npts)
        lpr_Y_ETA       => self % opo_info_gr (li_grids) % opr_points(D_ETA, Y, 1:li_npts)

        lpr_J           = lpr_X_XI * lpr_Y_ETA - lpr_X_ETA * lpr_Y_XI

        lpr_XI_X        =   lpr_Y_ETA / lpr_J
        lpr_XI_Y        = - lpr_X_ETA / lpr_J
        lpr_ETA_X       = - lpr_Y_XI  / lpr_J
        lpr_ETA_Y       =   lpr_X_XI  / lpr_J

        lpr_Phi         => lp_bbox % opr_dBasis (0    , 0:li_nen-1, 1:li_npts)
        lpr_Phi_XI      => lp_bbox % opr_dBasis (D_XI , 0:li_nen-1, 1:li_npts)
        lpr_Phi_ETA     => lp_bbox % opr_dBasis (D_ETA, 0:li_nen-1, 1:li_npts)
        ! ...        

        ! ...
        ! assembling B
        ! ...
        self % opo_pBasis (ai_space) % opr_B (1:li_nen, 1:li_npts) = lpr_Phi (1:li_nen,1:li_npts)
!        self % opo_pBasis (ai_space) % opr_B (1:li_nen, 1:li_npts) = lp_bbox % opr_dBasis (0,0:li_nen-1,1:li_npts)
        ! ...

        ! ...
        ! assembling grad B = ( dB/dx, dB/dy, dB/dz ) 
        ! ...
        DO li_b =1, li_nen
        DO li_pt=1, li_npts
            lpr_gradB (1:li_dim, li_b, li_pt) = MATMUL(lpr_invJacobian(1:li_dim, 1:li_dim, li_pt) &
            , lp_bbox % opr_dBasis (1:li_dim, li_b-1, li_pt))
        END DO
        END DO
        ! ...

        lpr_Phi_X      => lpr_gradB (X, :, :)
        lpr_Phi_Y      => lpr_gradB (Y, :, :)

        ! ...
        ! assembling curl B = 
        ! ...
        lpr_curlB (1, :, :) =   lpr_gradB (Y, :, :)
        lpr_curlB (2, :, :) = - lpr_gradB (X, :, :)
        ! ...

        IF (li_nderiv < 2) THEN
                RETURN
        END IF

        ! ...
        ! 
        ! ...
        lpr_X_XIXI      => self % opo_info_gr (li_grids) % opr_points(D_XIXI  , X, 1:li_npts)
        lpr_X_XIETA     => self % opo_info_gr (li_grids) % opr_points(D_XIETA , X, 1:li_npts)
        lpr_X_ETAETA    => self % opo_info_gr (li_grids) % opr_points(D_ETAETA, X, 1:li_npts)
        lpr_Y_XIXI      => self % opo_info_gr (li_grids) % opr_points(D_XIXI  , Y, 1:li_npts)
        lpr_Y_XIETA     => self % opo_info_gr (li_grids) % opr_points(D_XIETA , Y, 1:li_npts)
        lpr_Y_ETAETA    => self % opo_info_gr (li_grids) % opr_points(D_ETAETA, Y, 1:li_npts)
        ! ...   

        ! ...
        lpr_Phi_XIXI    => lp_bbox % opr_dBasis (D_XIXI,0:li_nen-1,1:li_npts)
        lpr_Phi_XIETA   => lp_bbox % opr_dBasis (D_XIETA,0:li_nen-1,1:li_npts)
        lpr_Phi_ETAETA  => lp_bbox % opr_dBasis (D_ETAETA,0:li_nen-1,1:li_npts)        
        ! ...

        ! ...
        DO li_b =1, li_nen
        DO li_pt=1, li_npts        
        lpr_C1(li_b, li_pt) = lpr_PHI_XIXI(li_b, li_pt)   &
        - lpr_X_XIXI(li_pt) * lpr_Phi_X(li_b, li_pt) &
        - lpr_Y_XIXI(li_pt) * lpr_PHI_Y(li_b, li_pt)

        lpr_C2(li_b, li_pt) = lpr_PHI_XIETA(li_b, li_pt)  &
        - lpr_X_XIETA(li_pt)  * lpr_Phi_X(li_b, li_pt) &
        - lpr_Y_XIETA(li_pt)  * lpr_PHI_Y(li_b, li_pt)

        lpr_C3(li_b, li_pt) = lpr_PHI_ETAETA(li_b, li_pt) &
        - lpr_X_ETAETA(li_pt) * lpr_Phi_X(li_b, li_pt) &
        - lpr_Y_ETAETA(li_pt) * lpr_PHI_Y(li_b, li_pt)
        END DO
        END DO
        ! ...

        ! ...
        DO li_b =1, li_nen
        lpr_HessianB (D_XX, li_b, :) = lpr_Y_ETA**2 * lpr_C1(li_b, :) &
        - 2 * lpr_Y_XI * lpr_Y_ETA * lpr_C2(li_b, :) &
        + lpr_Y_XI**2 * lpr_C3 (li_b, :)

        lpr_HessianB (D_XY, li_b, :) = -lpr_X_ETA * lpr_Y_ETA * lpr_C1(li_b, :) &
        + (lpr_X_XI * lpr_Y_ETA + lpr_X_ETA * lpr_Y_XI) * lpr_C2(li_b, :) &
        - lpr_X_XI * lpr_Y_XI * lpr_C3 (li_b, :)

        lpr_HessianB (D_YY, li_b, :) = lpr_X_ETA**2 * lpr_C1(li_b, :) &
        - 2 * lpr_X_XI * lpr_X_ETA * lpr_C2(li_b, :) &
        + lpr_X_XI**2 * lpr_C3 (li_b, :)

        DO li_d = 1, li_NPARAM
        lpr_HessianB(li_d, li_b, :) = lpr_HessianB(li_d, li_b, :) / lpr_J**2
        END DO
        END DO
        ! ...
!        print*,"assembly_physical_basis_2D: End"
        CALL printlog("assembly_physical_basis_2D: End", ai_dtllevel = 1)

    end subroutine assembly_physical_basis_2D    
    !---------------------------------------------------------------
    subroutine assembly_physical_basis_3D( self, ao_FEM, ai_space, ai_patch, ai_elt, ai_npts)
        IMPLICIT NONE
        TYPE(ASSEMBLY) :: self
        TYPE(FEM) :: ao_FEM
        INTEGER, INTENT(IN)  :: ai_space
        INTEGER, INTENT(IN)  :: ai_patch
        INTEGER, INTENT(IN)  :: ai_elt
        INTEGER, INTENT(IN)  :: ai_npts
        ! LOCAL VARIABLES

        CALL printlog("assembly_physical_basis_3D: Begin", ai_dtllevel = 1)

        PRINT *, 'assembly_physical_basis_3D: Not yet implemented'
        STOP

        CALL printlog("assembly_physical_basis_3D: End", ai_dtllevel = 1)

    end subroutine assembly_physical_basis_3D    
    !---------------------------------------------------------------
    subroutine assembly_space_points_advanced ( self, ao_FEM, ai_space, ai_patch, ai_nderiv, apr_points )
        IMPLICIT NONE
        TYPE(ASSEMBLY) :: self
        TYPE(FEM) :: ao_FEM
        INTEGER, INTENT(IN)  :: ai_space
        INTEGER, INTENT(IN)  :: ai_patch
        INTEGER, INTENT(IN)  :: ai_nderiv
        REAL*8 , DIMENSION(:,:,:,0:), INTENT(OUT)  :: apr_points
        ! LOCAL VARIABLES
        INTEGER  :: li_dim
        INTEGER  :: li_elt
        INTEGER  :: li_nel
        INTEGER  :: li_map
        INTEGER  :: li_grids
        INTEGER  :: li_npts
        INTEGER  :: li_REALelt
        INTEGER  :: li_pt
        INTEGER  :: li_deriv
        TYPE(GRID_DATA), POINTER :: lp_grid
        TYPE(BBOX), POINTER :: lp_bbox
        TYPE(MAPPING), POINTER :: lp_mapping
        TYPE(SPACE), POINTER :: lp_space

        CALL printlog("assembly_space_points_advanced : Begin", ai_dtllevel = 1)

!        PRINT *, "assembly_space_points_advanced : Begin"
        
        li_grids = ao_FEM % opi_InfoSpace ( ai_space, INFOSPACE_GRIDS )

        li_dim = ao_FEM % opi_dim(li_grids)

        lp_bbox => self % opo_bbox_sp ( ai_space )
        lp_space => ao_FEM % opo_spaces ( ai_space )

        li_nel  = ao_FEM % opi_InfoPatch ( li_grids, ai_patch, INFOPATCH_NEL)

        lp_grid => ao_FEM % opo_grids(li_grids) % opo_grid ( ai_patch )

        DO li_elt = 1, li_nel
!        PRINT *, ">> ELT ", li_elt
            li_npts = lp_grid % opo_elts ( li_elt ) % oi_npts

            CALL update_info_metric( self, ao_FEM, ai_space, ai_patch, li_elt, li_npts )
            CALL update_info_space ( self, ao_FEM, ai_space, ai_patch, li_elt, li_npts )

            DO li_deriv = 0, ai_nderiv
                apr_points (li_elt,1:li_npts,1:li_dim, li_deriv) = TRANSPOSE( &
                self % opo_info_gr (li_grids) % opr_points (li_deriv,1:li_dim,1:li_npts) )
            END DO
        END DO

!        PRINT *, "assembly_space_points_advanced : End"

        CALL printlog("assembly_space_points_advanced : End", ai_dtllevel = 1)

    end subroutine assembly_space_points_advanced
!----------------------------------------------------------------------------------------------
    subroutine set_metric_points_ass ( self, ao_FEM, ai_ref, ai_patch, apr_F, apr_DF, ai_nel, ai_npts, ai_Rd, ai_Rd2 )
        IMPLICIT NONE
        TYPE(ASSEMBLY) :: self
        TYPE(FEM) :: ao_FEM
        INTEGER, INTENT(IN)  :: ai_ref  !THE METRIC ID
        INTEGER, INTENT(IN)  :: ai_patch  !THE ai_patch ID
        INTEGER, INTENT(IN)  :: ai_nel   !NUMBER OF ELEMENTS
        INTEGER, INTENT(IN)  :: ai_npts   !NUMBER OF POINTS
        INTEGER, INTENT(IN)  :: ai_Rd   !THE DOMAIN DIMENSION
        INTEGER, INTENT(IN)  :: ai_Rd2  !GENERALLY THIS ai_Rd**2
        REAL*8 , DIMENSION(ai_nel, ai_npts,ai_Rd), INTENT(IN)  :: apr_F
        REAL*8 , DIMENSION(ai_nel, ai_npts,ai_Rd2), INTENT(IN)  :: apr_DF
        ! LOCAL VARIABLES
        INTEGER  :: li_elt
        INTEGER  :: li_i
        INTEGER  :: li_j
        INTEGER  :: li_index
        INTEGER  :: li_der

        CALL printlog("set_metric_points_ass : Begin", ai_dtllevel = 1)
!        print *, 'ai_nel, ai_npts,ai_Rd=', ai_nel, ai_npts,ai_Rd
!        print *, 'ai_nel, ai_npts,ai_Rd2=', ai_nel, ai_npts,ai_Rd2
!        print *, 'SIZE(apr_F,1)=', SIZE(apr_F,1)
!        print *, 'SIZE(apr_F,2)=', SIZE(apr_F,2)
!        print *, 'SIZE(apr_F,3)=', SIZE(apr_F,3)
        DO li_elt = 1, ai_nel
!            print *, 'li_elt = ', li_elt
            DO li_i = 1, ai_npts
!                print *, 'li_i = ', li_i
                ! position
                ao_FEM % opo_metrics (ai_ref) % opr_points(ai_patch, li_elt,li_i, 0, 1:ai_Rd) = apr_F (li_elt, li_i, 1:ai_Rd)
                ! derivatives
                li_index = 1
                DO li_j = 1, ai_Rd
!                    print *, 'li_j = ', li_j
                    DO li_der = 1, ai_Rd
!                        print *, 'li_der = ', li_der
!                        print *, 'li_index = ', li_index
                        ao_FEM % opo_metrics (ai_ref) % opr_points(ai_patch, li_elt,li_i, li_der, li_j) = &
                        apr_DF (li_elt, li_i, li_index)
                        li_index = li_index + 1
                    END DO
                END DO

            END DO
        END DO

        CALL printlog("set_metric_points_ass : End", ai_dtllevel = 1)

    end subroutine set_metric_points_ass
!----------------------------------------------------------------------------------------------
    subroutine set_metric_points_advanced_ass (self, ao_FEM, ai_ref, ai_patch, apr_points, ai_nel, ai_npts, ai_Rd, ai_nderivp1 )
        IMPLICIT NONE
        TYPE(ASSEMBLY) :: self
        TYPE(FEM) :: ao_FEM
        INTEGER, INTENT(IN)  :: ai_ref  !THE METRIC ID
        INTEGER, INTENT(IN)  :: ai_patch  !THE ai_patch ID
        INTEGER, INTENT(IN)  :: ai_nel   !NUMBER OF ELEMENTS
        INTEGER, INTENT(IN)  :: ai_npts   !NUMBER OF POINTS
        INTEGER, INTENT(IN)  :: ai_Rd   !THE DOMAIN DIMENSION
        INTEGER, INTENT(IN)  :: ai_nderivp1  !GENERALLY THIS ai_Rd
        REAL*8 , DIMENSION(ai_nel, ai_npts,ai_Rd, 0:ai_nderivp1-1), INTENT(IN)  :: apr_points
        ! LOCAL VARIABLES
        INTEGER  :: li_elt
        INTEGER  :: li_i
        INTEGER  :: li_nderiv

        CALL printlog("set_metric_points_advanced_ass : Begin", ai_dtllevel = 1)

        li_nderiv = ai_nderivp1 - 1

        DO li_elt = 1, ai_nel
            DO li_i = 1, ai_npts
                ao_FEM % opo_metrics (ai_ref) % opr_points(ai_patch, li_elt, li_i, 0:li_nderiv, 1:ai_Rd) = &
                TRANSPOSE(apr_points (li_elt, li_i, 1:ai_Rd, 0:li_nderiv))
!                PRINT *, 'Points : ', apr_points (li_elt, li_i, 1:ai_Rd, 0:li_nderiv)
            END DO
        END DO

        CALL printlog("set_metric_points_advanced_ass : End", ai_dtllevel = 1)

    end subroutine set_metric_points_advanced_ass
!----------------------------------------------------------------------------------------------
    subroutine get_space_parametric_points ( self, ao_FEM, ai_space, ai_patch, apr_points, ai_nelts, ai_npts, ai_dim )
        IMPLICIT NONE
        TYPE(ASSEMBLY) :: self
        TYPE(FEM) :: ao_FEM
        INTEGER, INTENT(IN)  :: ai_space
        INTEGER, INTENT(IN)  :: ai_patch
        INTEGER, INTENT(IN)  :: ai_nelts
        INTEGER, INTENT(IN)  :: ai_npts
        INTEGER, INTENT(IN)  :: ai_dim
        REAL*8 , DIMENSION(ai_nelts,ai_npts,ai_dim), INTENT(OUT)  :: apr_points
        ! LOCAL VARIABLES
        INTEGER  :: li_elt
        INTEGER  :: li_nel
        INTEGER  :: li_patch
        INTEGER  :: li_map
        INTEGER  :: li_grids
        INTEGER  :: li_npts
        INTEGER  :: li_REALelt
        INTEGER  :: li_pt
        INTEGER  :: li_dim
        INTEGER  :: li_ni
        INTEGER  :: li_nj
        INTEGER  :: li_i
        INTEGER  :: li_j
        REAL*8   :: lr_ksi
        REAL*8   :: lr_eta
        TYPE(GRID_DATA), POINTER :: lp_grid
        TYPE(BBOX), POINTER :: lp_bbox
        TYPE(MAPPING), POINTER :: lp_mapping
        TYPE(SPACE), POINTER :: lp_space

        CALL printlog("get_space_parametric_points : Begin", ai_dtllevel = 1)

        li_grids = ao_FEM % opi_InfoSpace ( ai_space, INFOSPACE_GRIDS )

        lp_bbox => self % opo_bbox_sp ( ai_space )
        lp_space => ao_FEM % opo_spaces ( ai_space )

        li_patch = ai_patch

        li_nel  = ao_FEM % opi_InfoPatch ( li_grids, li_patch, INFOPATCH_NEL)

        lp_grid => ao_FEM % opo_grids(li_grids) % opo_grid ( li_patch )

        li_dim = lp_grid % oi_dim
        IF  (li_dim > 2) THEN
            PRINT *, 'ERROR get_space_parametric_points : Only 1D and 2D case is treated'
            STOP
        END IF

        IF  (li_dim == 1) THEN
        DO li_elt = 1, li_nel

            li_pt = 1

            li_ni = lp_grid % opo_elts (li_elt) % opi_npts (1)

            DO li_i = 1, li_ni

                lr_ksi = lp_grid % opo_elts (li_elt) % opr_pts(li_i, 1)

                apr_points (li_elt, li_pt, 1) = lr_ksi

                li_pt = li_pt + 1

            END DO

        END DO
        END IF

        IF  (li_dim == 2) THEN
        DO li_elt = 1, li_nel

            li_pt = 1

            li_ni = lp_grid % opo_elts (li_elt) % opi_npts (1)
            li_nj = lp_grid % opo_elts (li_elt) % opi_npts (2)

            DO li_j = 1, li_nj

                lr_eta = lp_grid % opo_elts (li_elt) % opr_pts(li_j, 2)

                DO li_i = 1, li_ni

                    lr_ksi = lp_grid % opo_elts (li_elt) % opr_pts(li_i, 1)

                    apr_points (li_elt, li_pt, 1) = lr_ksi
                    apr_points (li_elt, li_pt, 2) = lr_eta

                    li_pt = li_pt + 1

                END DO

            END DO

        END DO
        END IF

        CALL printlog("get_space_parametric_points : End", ai_dtllevel = 1)

    end subroutine get_space_parametric_points
!----------------------------------------------------------------------------------------------
    subroutine Assembly_initialize_MURGE ( self, ao_FEM)
        IMPLICIT NONE
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
        INTEGER :: li_composed
        INTEGER :: li_matrix
        INTEGER :: li_i
        INTEGER :: li_nj
        LOGICAL :: ll_tensor
#ifdef _MURGE
        INTEGER(KIND=MURGE_INTS_KIND)          :: ierr
        INTEGER       :: Me
        INTEGER       :: NTasks
        INTEGER       :: required
        INTEGER       :: provided
        ! MURGE Data
        INTEGER(KIND=MURGE_INTS_KIND)                         :: id
        INTEGER(KIND=MURGE_INTS_KIND)                         :: m
        ! CSC Data
        INTEGER(KIND=MURGE_INTS_KIND)                         :: n
        INTEGER(KIND=MURGE_INTS_KIND)                         :: dof
        INTEGER(KIND=MURGE_INTL_KIND)                         :: nnzeros
        INTEGER(KIND=MURGE_INTL_KIND)                         :: edgenbr
        REAL(KIND=MURGE_COEF_KIND)                             :: val
        ! Local Data
        INTEGER(KIND=MURGE_INTS_KIND)                         :: root
        INTEGER(KIND=MURGE_INTS_KIND)                         :: base
        ! Other data
        REAL(KIND=8)                                                :: prec
        INTEGER(KIND=MURGE_INTS_KIND)                         :: i
        INTEGER(KIND=MURGE_INTS_KIND)                         :: j
        INTEGER(KIND=MURGE_INTS_KIND)                         :: k
        INTEGER(KIND=MURGE_INTS_KIND)                         :: l
        INTEGER(KIND=MURGE_INTS_KIND)                         :: solver
        INTEGER(KIND=MURGE_INTS_KIND)                         :: zero=0
        INTEGER(KIND=MURGE_INTS_KIND)                         :: one=1
        INTEGER(KIND=MURGE_INTS_KIND)                         :: two=2
#endif
        
        CALL printlog("Assembly_initialize_MURGE : Begin", ai_dtllevel = 1)

#ifdef _MURGE

        ! Set Options
        prec = 1e-7
        root = -1
        dof = 1
        base = 1
        
    ! ********************************************************
    !                MPI AND MURGE INTIALIZATION
    ! ********************************************************
        required = MPI_THREAD_MULTIPLE
        CALL MPI_Init_thread(required, provided, ierr)

        CALL MPI_COMM_SIZE(MPI_Comm_world, NTasks, ierr)
        CALL MPI_COMM_RANK(MPI_Comm_world, Me, ierr)

        ! Starting MURGE
        CALL MURGE_INITIALIZE(ao_FEM % oi_nMatrices, ierr)
        
        IF (ierr /= MURGE_SUCCESS) CALL abort()
    ! ********************************************************

    ! ********************************************************
    !                    CHOOSE THE SOLVER
    ! ********************************************************
        solver = MURGE_SOLVER_PASTIX
    ! ********************************************************

        DO id = 0, ao_FEM % oi_nMatrices - 1

        ! ********************************************************
        !                   PASTIX INTIALIZATION
        ! ********************************************************
            IF (solver == MURGE_SOLVER_PASTIX) THEN
                CALL MURGE_SetDefaultOptions(id, zero, ierr)
                CALL MURGE_SetOptionINT(id, IPARM_VERBOSE, API_VERBOSE_NOT, ierr)
                CALL MURGE_SetOptionINT(id, IPARM_MATRIX_VERIFICATION, API_YES, ierr)
            END IF
        ! ********************************************************
            
        ! ********************************************************
        !              MATRIX INITIALIZATION - PARAMS
        ! ********************************************************
            CALL MURGE_SetOptionINT(id, MURGE_IPARAM_DOF, dof, ierr)
            CALL MURGE_SetOptionINT(id, MURGE_IPARAM_SYM, MURGE_BOOLEAN_FALSE, ierr)
            CALL MURGE_SetOptionINT(id, MURGE_IPARAM_BASEVAL, base, ierr)
            
            CALL MURGE_SetOptionREAL(id, MURGE_RPARAM_EPSILON_ERROR, PREC, ierr)
        ! ********************************************************

        ! ********************************************************
        !              MATRIX INITIALIZATION : GRAPH
        ! ********************************************************
!            IF (Me == 0) THEN
!                edgenbr = ao_FEM % opo_CCR_Matrix(id) % oi_nel
!                n         = ao_FEM % opo_CCR_Matrix(id) % oi_nR
!
!                CALL MURGE_GRAPHBEGIN(id, n, edgenbr, ierr)
!                DO k = 1, edgenbr
!                    i = ao_FEM % opo_CCR_Matrix(id) % opi_ind_row(k)
!                    j = ao_FEM % opo_CCR_Matrix(id) % opi_ind_col(k)
!                    CALL MURGE_GRAPHEDGE(id, i, j, ierr)
!                END DO
!            ELSE
!                edgenbr = 0
!                CALL MURGE_GRAPHBEGIN(id, n, edgenbr,  ierr)
!            END IF
!
!            CALL MURGE_GRAPHEND(id, ierr)
!            
!            IF (ierr /= MURGE_SUCCESS) CALL abort()
        ! ********************************************************

        END DO
        
#endif

        CALL printlog("Assembly_initialize_MURGE : End", ai_dtllevel = 1)
        return 

    end subroutine Assembly_initialize_MURGE
end module assembly_tools


