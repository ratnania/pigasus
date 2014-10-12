!     
! File:   mass_t2.F90
! Author: root
!
! Created on January 15, 2013, 5:39 PM
!

!
! File:   mass.F90
! Author: root
!
! Created on January 16, 2012, 9:06 AM
!

module mass_t2_module
    use used_precision
    use tracelog_module
    use grids_def
    use geometries_def
    use geometry_tools
    use bbox_def
    use bbox_module
    use fem_def
    use assembly_def
    implicit none

    private

    public :: build_Mass_Local_t2

#ifdef _DEBUG
    integer, parameter, private  :: mi_dtllevel_base = 0
#else
    integer, parameter, private  :: mi_dtllevel_base = 2
#endif

contains

    !----------------------------------------------------------------------------------------------
    subroutine build_Mass_Local_t2(ao_ASS, ao_FEM, ai_grids_id, ai_id, ai_elt, ai_matrix)
        implicit none
        TYPE(ASSEMBLY) :: ao_ASS
        type(FEM) :: ao_FEM
        !> grids id
        integer :: ai_grids_id
        !> patch id
        integer :: ai_id
        !> current element
        integer :: ai_elt
        !> current matrix-id
        integer :: ai_matrix
        ! LOCAL VARIABLES
        integer :: li_npts
        integer :: li_sp, li_spprime
        integer :: li_bloc, li_bprimeloc
        integer :: li_map
        INTEGER :: li_loc_id
        INTEGER :: li_nparam
        INTEGER :: li_annz
        INTEGER :: li_rnnz
        INTEGER :: li_arbnnz
        INTEGER :: li_blocn
        real(wp) :: lr_contribution
        real(wp), dimension(ao_ASS % oi_maxnpts) :: lpr_jacobian
        real(wp), DIMENSION(ao_ASS % oi_maxnpts) :: lpr_contribution
        TYPE(GRID_DATA), pointer :: lp_grid
        TYPE(BBOX), POINTER :: lp_bbox
        real(wp), dimension(:), POINTER :: lpr_rho
        real(wp), dimension(:), pointer :: lpr_a
        real(wp), dimension(:), pointer :: lpr_b
        real(wp), dimension(:), pointer :: lpr_x

#ifdef _DEBUG
        call printlog("build_Mass_Local_t2 : Begin", ai_dtllevel = mi_dtllevel_base + 1)
#endif
!print *, 'ai_elt=', ai_elt
        print *, 'build_Mass_Local_t2: Not yet implemented'
        STOP
        !\todo prblm du id => id+1
        lp_grid => ao_FEM % opo_grids ( ai_grids_id ) % opo_grid ( ai_id-1 )


        li_npts = lp_grid % opo_elts ( ai_elt ) % oi_npts

        li_sp       = ao_FEM % opi_InfoMatrix (ai_matrix, INFOMATRIX_SPACE_1)
        li_spprime  = ao_FEM % opi_InfoMatrix (ai_matrix, INFOMATRIX_SPACE_2)

!        print*,'current spaces are ', li_sp, li_spprime

        li_loc_id = ao_FEM % opi_InfoMatrix(ai_matrix, INFOMATRIX_LOCID)
        li_nparam = ao_FEM % opi_InfoMatrix(ai_matrix, INFOMATRIX_NPARAM)

        lpr_jacobian = ao_ASS % opo_info_gr (ai_grids_id) % opr_jacobians
        lpr_jacobian = dabs(lpr_jacobian)
!        !> \todo a supprimer apres
!        lpr_jacobian = 1.0_wp

#ifdef _TRACE
        CALL printlog("TODO : mettre oo_assl dans info_gr" &
        , ai_dtllevel = mi_dtllevel_base + 1)
#endif

        lp_bbox => ao_ASS % opo_bbox_sp ( li_sp )
        
        li_annz = assl_get_annz(lp_bbox % oo_assl)
        li_rnnz = assl_get_rnnz(lp_bbox % oo_assl)
        li_arbnnz = assl_get_arbnnz(lp_bbox % oo_assl)
        li_blocn = assl_get_blocn(lp_bbox % oo_assl)

        ! ...
        ! updating the parameter function (the R matrix)
        ! ...
        ALLOCATE(lpr_a(li_annz))
        ALLOCATE(lpr_rho(li_rnnz))
        ALLOCATE(lpr_b(li_blocn))
        ALLOCATE(lpr_x(li_arbnnz))

        lpr_jacobian = 1.0
        lp_grid % opo_elts ( ai_elt ) % opr_w (1:li_npts) = 1.0
        lp_grid % opo_elts (ai_elt) % opr_values_m(li_loc_id, 1, 1 : li_npts) = 1.0
!        print *, 'lpr_jacobian = ', lpr_jacobian
!        print *, 'lpr_w = ', lp_grid % opo_elts ( ai_elt ) % opr_w (1:li_npts)
!        print *, 'lpr_f = ', lp_grid % opo_elts (ai_elt) % opr_values_m(li_loc_id, 1, 1 : li_npts)

        lpr_rho (1:li_rnnz)  = lp_grid % opo_elts ( ai_elt ) % opr_w (1:li_rnnz) &
                                   * lp_grid % opo_elts (ai_elt) % opr_values_m(li_loc_id, 1, 1 : li_rnnz) &
                                   * lpr_jacobian (1:li_rnnz)

        call assl_update_r(lp_bbox % oo_assl, lpr_rho)

        lpr_a(1:li_annz) = lp_bbox % opr_gBasis_t(0, 1:li_annz)
        lpr_b(:) = lp_bbox % opr_lBasis_t(0, 1:li_blocn, ai_elt)
        
        call assl_update_a(lp_bbox % oo_assl, lpr_a)
        call assl_update_bloc(lp_bbox % oo_assl, lpr_b)
        ! ...

        ! ...
        ! assembling the whole line
        ! ...
        call assl_compute_tarb(lp_bbox % oo_assl)
        ! ...

        CALL assl_get_tarb(lp_bbox % oo_assl, lpr_x)

        ao_ASS % opr_Matrix_eltLine (1:li_arbnnz) = lpr_x(1:li_arbnnz)

!        print *, 'Matrix_eltLine (:)= ', ao_ASS % opr_Matrix_eltLine (1:li_arbnnz)

!        STOP

        DEALLOCATE(lpr_a)
        DEALLOCATE(lpr_rho)
        DEALLOCATE(lpr_b)
        DEALLOCATE(lpr_x)

#ifdef _DEBUG
        call printlog("build_Mass_Local_t2 : End", ai_dtllevel = mi_dtllevel_base + 1)
#endif
    end subroutine build_Mass_Local_t2

end module mass_t2_module
