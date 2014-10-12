!     
! File:   H1_norm.F90
! Author: root
!
! Created on October 12, 2012, 11:48 AM
!

module H1_norm_module
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

    public :: build_H1_norm_Local

#ifdef _DEBUG
    integer, parameter, private  :: mi_dtllevel_base = 0
#else
    integer, parameter, private  :: mi_dtllevel_base = 2
#endif
    integer, parameter, private  :: mi_MAXDIM = 3

contains

    !----------------------------------------------------------------------------------------------
    subroutine build_H1_norm_Local(ao_ASS, ao_FEM, ai_grids_id, ai_id, ai_elt   &
        , ai_norm, apr_norm_elt)
        implicit none
        TYPE(ASSEMBLY) :: ao_ASS
        type(FEM) :: ao_FEM
        !> grids id
        integer :: ai_grids_id
        !> patch id
        integer :: ai_id
        !> current element
        integer :: ai_elt
        !> current norm-id
        integer :: ai_norm
        real(wp), dimension(0:) :: apr_norm_elt
        ! LOCAL VARIABLES
        integer :: li_npts
        integer :: li_field
        integer :: li_b
        integer :: li_dim
        integer :: li_iparam
        integer :: li_d
        integer :: li_dprime
        integer :: li_i
        integer :: li_A
        integer :: li_sp
        integer :: li_bloc
        integer :: li_ndof
        integer :: li_map
        INTEGER :: li_field_loc_id
        INTEGER :: li_norm_loc_id
        INTEGER :: li_nparam
        INTEGER :: li_conelt
        real(wp), dimension (ao_ASS % oi_maxndof, mi_MAXDIM,ao_ASS % oi_maxnpts) :: lpr_fieldh
        real(wp), dimension (ao_ASS % oi_maxndof, mi_MAXDIM, ao_ASS % oi_maxnpts) :: lpr_gradB
        real(wp), DIMENSION(ao_ASS % oi_maxnpts) :: lpr_contribution
        real(wp), dimension(ao_ASS % oi_maxnpts) :: lpr_jacobian
        real(wp), dimension(ao_ASS % oi_maxnpts) :: lpr_w
        TYPE(CONNECTIVITY), pointer :: lp_con
        TYPE(GRID_DATA), pointer :: lp_grid
        TYPE(PHYSICAL_BASIS), pointer :: lp_pBasis

#ifdef _DEBUG
        call printlog("build_H1_norm_Local : Begin", ai_dtllevel = mi_dtllevel_base + 1)
#endif
!print *, 'ai_elt=', ai_elt

        li_dim = ao_FEM % opi_dim(ai_grids_id)
        lpr_fieldh = 0.0_wp

        !\todo prblm du id => id+1
        lp_grid => ao_FEM % opo_grids ( ai_grids_id ) % opo_grid ( ai_id-1 )
        li_npts = lp_grid % opo_elts ( ai_elt ) % oi_npts

        li_norm_loc_id = ao_FEM % opi_InfoNorm(ai_norm, INFONORM_LOCID)
        li_field = ao_FEM % opi_InfoNorm(ai_norm, INFONORM_FIELD)
        li_sp       = ao_FEM % opi_InfoField(li_field, INFOFIELD_SPACE )
        li_ndof     = ao_FEM % opi_InfoField(li_field, INFOFIELD_NDOF )
        li_field_loc_id   = ao_FEM % opi_InfoField(li_field, INFOFIELD_LOCID )
!        li_nparam   = ao_FEM % opi_InfoField(li_field, INFOFIELD_NPARAM )

        lp_con          => ao_FEM % opo_spaces (li_sp) % oo_con
        lp_pBasis       => ao_ASS % opo_pBasis (li_sp)

        li_conelt = lp_con % opi_real_elts (ai_id-1,ai_elt)

        !> \todo attention +1
        do li_b = 1, lp_con % opi_nen(ai_id)

            li_A = lp_con % opi_LM(ai_id, li_b, li_conelt)

            if (li_A == 0) then
                cycle
            end if

            !> \todo traiter le cas ou ndof > 1
            ! apr_Basis (1:li_ndof,1:li_npts)
            !> \todo on doit enlever 1 a cause des indices
            !> bbox % dBasis commence a 0

            lpr_gradB (1,1:li_dim, 1:li_npts) = lp_pBasis % opr_gradB (1:li_dim, li_b, 1:li_npts)

            DO li_d = 1, li_dim
                lpr_fieldh(1, li_d, 1:li_npts) = lpr_fieldh(1, li_d, 1:li_npts) &
                + ao_FEM % opo_F(li_field) % opr_c(li_A) * lpr_gradB (1,li_d,1:li_npts)
            END DO
            
        end do

        ! if we don't use an exterior mapping, then we must use the transformation defined
        ! by the current geometry
        if ( ao_FEM % opi_InfoSpace ( li_sp , INFOSPACE_EXTMAPPING ) == 1 ) then
            li_map = ao_FEM % opi_InfoSpace ( li_sp , INFOSPACE_MAPPING )
            lpr_jacobian = ao_ASS % opo_info_mp (li_map) % opr_jacobians
        else
            lpr_jacobian = ao_ASS % opo_info_sp (li_sp) % opr_jacobians
        end if
        lpr_jacobian = dabs(lpr_jacobian)

!        lpr_contribution = 0.0_wp
!        DO li_d = 1, li_dim
!            lpr_contribution (1:li_npts) = lpr_contribution (1:li_npts) + (  lpr_fieldh (1, li_d, 1:li_npts)    &
!                                            - lp_grid % opo_elts (ai_elt) % opr_values_f(li_field_loc_id, li_d, 1 : li_npts) )**2
!        END DO

        lpr_contribution = 0.0_wp
        li_iparam = 1
        DO li_dprime = 1, li_dim
            DO li_d = 1, li_dim
                ! multiply by the dot product of the two basis
                lpr_contribution (1:li_npts) = lpr_contribution (1:li_npts) &
                + ( lpr_fieldh (1, li_d, 1:li_npts) &
                - lp_grid % opo_elts (ai_elt) % opr_values_f(li_field_loc_id, li_d, 1 : li_npts) ) &
                * ( lpr_fieldh (1, li_dprime, 1:li_npts) &
                - lp_grid % opo_elts (ai_elt) % opr_values_f(li_field_loc_id, li_dprime, 1 : li_npts) ) &
                * lp_grid % opo_elts (ai_elt) % opr_values_n(li_norm_loc_id, li_iparam, 1 : li_npts)
                li_iparam = li_iparam + 1
            END DO
        END DO

        lpr_contribution (1:li_npts) = lpr_contribution (1:li_npts) &
            ! multiply by quadratures weights
            * lp_grid % opo_elts ( ai_elt ) % opr_w (1:li_npts) &
            ! multiply by jacobians
            * lpr_jacobian (1:li_npts)

        apr_norm_elt(ai_norm) = SUM ( lpr_contribution (1:li_npts) )

!#ifdef _DEBUG
!        print *,"fieldh=", lpr_fieldh
!        print *,"field exact=", lp_grid % opo_elts (ai_elt) % opr_values_f(li_field_loc_id, 1, 1 : li_npts)
!        print *,"ar_norm_elt =", ar_norm_elt
!#endif

#ifdef _DEBUG
        call printlog("build_H1_norm_Local : End", ai_dtllevel = mi_dtllevel_base + 1)
#endif
    end subroutine build_H1_norm_Local

end module H1_norm_module

