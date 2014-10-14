!     
! File:   l2_projector.F90
! Author: root
!
! Created on January 27, 2012, 8:33 AM
!


module l2_projector_module
    use used_precision
    use tracelog_module
    use grids_def
    use geometries_def
    use geometry_tools
    use bbox_def
    use bbox_module
    use connectivities_def
    use fem_def
    use assembly_def
    implicit none

    private

    public :: build_l2_projector_Local

#ifdef _DEBUG
    integer, parameter, private  :: mi_dtllevel_base = 0
#else
    integer, parameter, private  :: mi_dtllevel_base = 2
#endif

contains

    !----------------------------------------------------------------------------------------------
    subroutine build_l2_projector_Local(ao_ASS, ao_FEM, ai_grids_id, ai_id, ai_elt, ai_field, apr_Projection_elt)
        implicit none
        TYPE(ASSEMBLY) :: ao_ASS
        type(FEM) :: ao_FEM
        !> grids id
        integer :: ai_grids_id
        !> patch id
        integer :: ai_id
        !> current element
        integer :: ai_elt
        !> current field-id
        integer :: ai_field
        real(wp), dimension(0:,:) :: apr_Projection_elt
        ! LOCAL VARIABLES
        integer :: li_npts
        integer :: li_b
        integer :: li_i
        integer :: li_A
        integer :: li_sp
        integer :: li_bloc
        integer :: li_ndof
        integer :: li_map
        INTEGER :: li_loc_id
        INTEGER :: li_nparam
        INTEGER :: li_conelt
        INTEGER :: li_dim
        real(wp) :: lr_contribution
        real(wp), dimension (ao_ASS % oi_maxndof,ao_ASS % oi_maxnpts) :: lpr_B
        real(wp), dimension(ao_ASS % oi_maxnpts) :: lpr_jacobian
        real(wp), DIMENSION(ao_ASS % oi_maxnpts) :: lpr_contribution
        real(wp), dimension(ao_ASS % oi_maxnpts) :: lpr_w
        TYPE(CONNECTIVITY), pointer :: lp_con
        TYPE(CONNECTIVITY), pointer :: lp_conprime
        TYPE(GRID_DATA), pointer :: lp_grid
        TYPE(PHYSICAL_BASIS), pointer :: lp_pBasis

#ifdef _DEBUG
        call printlog("build_l2_projector_Local : Begin", ai_dtllevel = mi_dtllevel_base + 1)
#endif
!print *, 'ai_elt=', ai_elt

        lpr_B = 0.0_wp
        li_dim = ao_FEM % opi_dim(ai_grids_id)

        !\todo prblm du id => id+1
        lp_grid => ao_FEM % opo_grids ( ai_grids_id ) % opo_grid ( ai_id-1 )
        li_npts = lp_grid % opo_elts ( ai_elt ) % oi_npts

        li_sp       = ao_FEM % opi_InfoField(ai_field, INFOFIELD_SPACE )
        li_ndof     = ao_FEM % opi_InfoField(ai_field, INFOFIELD_NDOF )
        li_loc_id   = ao_FEM % opi_InfoField(ai_field, INFOFIELD_LOCID )
!        li_nparam   = ao_FEM % opi_InfoField(ai_field, INFOFIELD_NPARAM )

        lp_con          => ao_FEM % opo_spaces (li_sp) % oo_con
        lp_pBasis       => ao_ASS % opo_pBasis (li_sp)
        
        lpr_jacobian = ao_ASS % opo_info_gr (ai_grids_id) % opr_jacobians
        lpr_jacobian = dabs(lpr_jacobian)

        li_conelt = lp_con % opi_real_elts (ai_id-1,ai_elt)

!        print *, '===='
!        print *, 'w = ',lp_grid % opo_elts ( ai_elt ) % opr_w (1:li_npts)

!print *, 'ai_id, ai_elt=', ai_id, ai_elt
!        print *, 'lpr_jacobian=', lpr_jacobian

!print *,'lp_con % opi_nen=',lp_con % opi_nen
!print *,'nen=',lp_con % opi_nen(ai_id, ai_elt)
!print *,'lp_con % opi_LM=',lp_con % opi_LM

        do li_b = 1, lp_con % opi_nen(ai_id)
!print *,'li_b=',li_b
            li_A = lp_con % opi_LM ( ai_id, li_b, li_conelt )

            if (li_A == 0) then
                cycle
            end if

            lpr_B (1,1:li_npts) = lp_pBasis % opr_B (li_b, 1:li_npts)
!print *, 'Basis=', lpr_B
!print *, 'values_f=', lp_grid % opo_elts (ai_elt) % opr_values_f(li_loc_id, 1, 1 : li_npts)

            !> \todo only work if ndof=1 (scalar unknowns)
            ! multiply by input data
            lpr_contribution (1:li_npts) = lp_grid % opo_elts (ai_elt) % opr_values_f(li_loc_id, 1, 1 : li_npts)   &
            ! multiply by the left basis
            * lpr_B (1,1:li_npts)   &
            ! multiply by quadratures weights
            * lp_grid % opo_elts ( ai_elt ) % opr_w (1:li_npts) &
            ! multiply by jacobians
            * lpr_jacobian (1:li_npts)

            lr_contribution = SUM ( lpr_contribution (1:li_npts) )

            apr_Projection_elt(ai_field, li_b) = apr_Projection_elt(ai_field, li_b) &
            + lr_contribution

        end do
!        print *, '===='
!#ifdef _DEBUG
!print*,'apr_Projection_elt',apr_Projection_elt(ai_field,1:lp_con % opi_nen (ai_id, li_conelt))
!print *,"field exact=", lp_grid % opo_elts (ai_elt) % opr_values_f(li_loc_id, 1:li_dim, 1 : li_npts)
!#endif
#ifdef _DEBUG
        call printlog("build_l2_projector_Local : End", ai_dtllevel = mi_dtllevel_base + 1)
#endif
    end subroutine build_l2_projector_Local

end module l2_projector_module

