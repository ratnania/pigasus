!     
! File:   second_deriv_s.F90
! Author: root
!
! Created on October 11, 2012, 9:05 AM
!


module second_deriv_s_field_module
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

    public :: build_second_deriv_s_field_Local

#ifdef _DEBUG
    integer, parameter, private  :: mi_dtllevel_base = 0
#else
    integer, parameter, private  :: mi_dtllevel_base = 2
#endif
    !> \todo to remove from here, and use an adequate variable
    integer, parameter, private :: mi_MAXDIM = 3
    integer, parameter, private :: mi_MAXDIM_invH = 3
    integer, parameter, private :: mi_NPARAM = 3

contains

    !----------------------------------------------------------------------------------------------
    subroutine build_second_deriv_s_field_Local(ao_ASS, ao_FEM, ai_grids_id, ai_id, ai_elt, ai_field, apr_fieldh)
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
        real*8, dimension(0:,:,:,:), intent(inout) :: apr_fieldh
        ! LOCAL VARIABLES
        integer :: li_npts
        integer :: li_b
        integer :: li_i
        integer :: li_iprime
        integer :: li_A
        integer :: li_sp
        integer :: li_bloc
        integer :: li_ndof
        integer :: li_map
        INTEGER :: li_loc_id
        INTEGER :: li_nparam
        INTEGER :: li_field_id
        INTEGER :: li_dim
        INTEGER :: li_conelt
        INTEGER :: li_pt
        INTEGER :: li_iparam
        real(wp) :: lr_contribution
        INTEGER, PARAMETER :: MAXDIM = 3
        real(wp), dimension(ao_ASS % oi_maxnpts) :: lpr_HessianB
        real(wp), dimension(ao_ASS % oi_maxnpts) :: lpr_jacobian
        real(wp), DIMENSION(ao_ASS % oi_maxnpts) :: lpr_contribution
        REAL(wp), DIMENSION (mi_NPARAM, ao_ASS % oi_maxnpts) :: lpr_A
        TYPE(CONNECTIVITY), pointer :: lp_con
        TYPE(GRID_DATA), pointer :: lp_grid
        TYPE(PHYSICAL_BASIS), pointer :: lp_pBasis

#ifdef _DEBUG
        call printlog("build_second_deriv_s_field_Local : Begin", ai_dtllevel = mi_dtllevel_base + 1)
#endif
!print *, 'ai_elt=', ai_elt

        li_dim = ao_FEM % opi_dim(ai_grids_id)
        lpr_HessianB = 0.0_wp

        !\todo prblm du id => id+1
        lp_grid => ao_FEM % opo_grids ( ai_grids_id ) % opo_grid ( ai_id-1 )
        li_npts = lp_grid % opo_elts ( ai_elt ) % oi_npts

        li_sp       = ao_FEM % opi_InfoField(ai_field, INFOFIELD_SPACE )
        li_ndof     = ao_FEM % opi_InfoField(ai_field, INFOFIELD_NDOF )
        li_loc_id   = ao_FEM % opi_InfoField(ai_field, INFOFIELD_LOCID )
        li_nparam   = ao_FEM % opi_InfoField(ai_field, INFOFIELD_NPARAM )

        li_field_id = ao_FEM % opi_InfoField(ai_field, INFOFIELD_OPERANDE )

        lp_con          => ao_FEM % opo_spaces (li_sp) % oo_con
        lp_pBasis       => ao_ASS % opo_pBasis (li_sp)

        li_conelt = lp_con % opi_real_elts (ai_id-1,ai_elt)

        lpr_jacobian = ao_ASS % opo_info_gr (ai_grids_id) % opr_jacobians
        lpr_jacobian = dabs(lpr_jacobian)

    !...
    ! compute the param matrix
    !...
        lpr_A (1:li_nparam, 1:li_npts) = lp_grid % opo_elts (ai_elt) % opr_values_f(li_loc_id, 1:li_nparam, 1:li_npts)
    !...

        !> \todo attention +1
        do li_b = 1, lp_con % opi_nen(ai_id)

            li_A = lp_con % opi_LM(ai_id, li_b, li_conelt)

            if (li_A == 0) then
                cycle
            end if

            DO li_pt=1, li_npts
            lpr_HessianB (li_pt) = DOT_PRODUCT(lpr_A(1:li_nparam, li_pt) &
            , lp_pBasis % opr_HessianB (1:li_nparam, li_b, li_pt))
            END DO

            apr_fieldh(ai_field,1,1, 1:li_npts) = apr_fieldh(ai_field,1,1, 1:li_npts) &
            + ao_FEM % opo_F(li_field_id) % opr_c(li_A) &
            * lpr_HessianB (1:li_npts)

        end do

#ifdef _DEBUG
        call printlog("build_second_deriv_s_field_Local : End", ai_dtllevel = mi_dtllevel_base + 1)
#endif
    end subroutine build_second_deriv_s_field_Local

end module second_deriv_s_field_module


