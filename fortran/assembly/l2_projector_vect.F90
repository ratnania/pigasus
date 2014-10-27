!     
! File:   l2_projector_vect_vect.F90
! Author: root
!
! Created on March 29, 2012, 6:25 PM
!

module l2_projector_vect_module
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

    public :: build_l2_projector_vect_Local

#ifdef _DEBUG
    integer, parameter, private  :: mi_dtllevel_base = 0
#else
    integer, parameter, private  :: mi_dtllevel_base = 2
#endif

contains

    !----------------------------------------------------------------------------------------------
    subroutine build_l2_projector_vect_Local(ao_ASS, ao_FEM, ai_grids_id, ai_id, ai_elt, ai_field, apr_Projection_elt)
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
        INTEGER :: li_cur_sp
        INTEGER :: li_cur_b
        INTEGER :: li_cur_nen
        integer :: li_cur_sp_id
        integer :: li_conelt
        real(wp) :: lr_contribution
        real(wp), dimension (ao_ASS % oi_maxndof,ao_ASS % oi_maxnpts) :: lpr_Basis
        real(wp), dimension(ao_ASS % oi_maxnpts) :: lpr_jacobian
        real(wp), DIMENSION(ao_ASS % oi_maxnpts) :: lpr_contribution
        real(wp), dimension(ao_ASS % oi_maxnpts) :: lpr_w
        real(wp), dimension(ao_ASS % oi_maxndof,ao_ASS % oi_maxnpts) :: lpr_f
        real(wp), dimension(ao_ASS % oi_maxndof) :: lpr_v ! the unit vector
        TYPE(CONNECTIVITY), pointer :: lp_con
        TYPE(CONNECTIVITY), pointer :: lp_conprime
        TYPE(GRID_DATA), pointer :: lp_grid

#ifdef _DEBUG
        call printlog("build_l2_projector_vect_Local : Begin", ai_dtllevel = mi_dtllevel_base + 1)
#endif
!print *, 'ai_elt=', ai_elt

        lpr_Basis = 0.0_wp

        !\todo prblm du id => id+1
        lp_grid => ao_FEM % opo_grids ( ai_grids_id ) % opo_grid ( ai_id-1 )
        li_npts = lp_grid % opo_elts ( ai_elt ) % oi_npts

        li_sp       = ao_FEM % opi_InfoField(ai_field, INFOFIELD_SPACE )
        li_ndof     = ao_FEM % opi_InfoField(ai_field, INFOFIELD_NDOF )
        li_loc_id   = ao_FEM % opi_InfoField(ai_field, INFOFIELD_LOCID )
!        li_nparam   = ao_FEM % opi_InfoField(ai_field, INFOFIELD_NPARAM )

        lp_con          => ao_FEM % opo_spaces (li_sp) % oo_con

!        print*,'current space is ', li_sp
!        lpr_b (1 : li_nparam, 1 : li_npts) =    &
!        lp_grid % opo_elts (ai_elt) % opr_values_m(li_loc_id, 1 : li_nparam, 1 : li_npts)

        ! if we don't use an exterior mapping, then we must use the transformation defined
        ! by the current geometry
!        if ( ao_FEM % opi_InfoSpace ( li_sp , INFOSPACE_EXTMAPPING ) == 1 ) then
!            li_map = ao_FEM % opi_InfoSpace ( li_sp , INFOSPACE_MAPPING )
!            lpr_jacobian = ao_ASS % opo_info_mp (li_map) % opr_jacobians
!        else
!            lpr_jacobian = ao_ASS % opo_info_sp (li_sp) % opr_jacobians
!        end if
!        lpr_jacobian = dabs(lpr_jacobian)
!        !> \todo a supprimer apres
        lpr_jacobian = 1.0_wp

        li_cur_sp = 1
        li_cur_sp_id = ao_FEM % opi_InfoComposedSpace(li_sp, li_cur_sp)
        li_cur_nen = ao_FEM % opo_spaces (li_cur_sp_id) % oo_con % opi_nen (ai_id, ai_elt)
        li_cur_b = 1
        lpr_v (:) = 0.0_wp ; lpr_v (li_cur_sp) = 1.0_wp

        li_conelt = lp_con % opi_real_elts (ai_id-1,ai_elt)
        do li_b = 1, lp_con % opi_nen (ai_id, li_conelt)
!print *,'li_b=',li_b
            li_A = lp_con % opi_LM ( ai_id, li_b, li_conelt )

            if (li_A == 0) then
                if ((li_cur_b == li_cur_nen) .AND. (li_b /= lp_con % opi_nen (ai_id, ai_elt)) ) then
                    li_cur_b = 0
                    !> \todo this works only for V = V1xV2
                    li_cur_sp = li_cur_sp + 1
                    li_cur_sp_id = ao_FEM % opi_InfoComposedSpace(li_sp, li_cur_sp)
                    li_cur_nen = ao_FEM % opo_spaces (li_cur_sp_id) % oo_con % opi_nen (ai_id, ai_elt)
                    lpr_v (:) = 0.0_wp ; lpr_v (li_cur_sp) = 1.0_wp
                end if
                li_cur_b = li_cur_b + 1
                cycle
            end if

            lpr_f (1:li_ndof, 1 : li_npts) = lp_grid % opo_elts (ai_elt) % opr_values_f(li_loc_id, 1:li_ndof, 1 : li_npts)

            CALL af_Basis(ao_ASS, ao_FEM, li_sp, li_ndof, ai_grids_id, ai_id, ai_elt, li_b  &
            , li_cur_sp_id, li_cur_b, li_cur_nen, lpr_Basis )
!print *, 'lpr_Basis=', lpr_Basis

            !> \todo only work if ndof=1 (scalar unknowns)
            ! multiply by input data
            lpr_contribution (1:li_npts) =  MATMUL (TRANSPOSE(lpr_f), lpr_v)   &
            ! multiply by the left basis
            * lpr_Basis (1,1:li_npts)   &
            ! multiply by quadratures weights
            * lp_grid % opo_elts ( ai_elt ) % opr_w (1:li_npts) &
            ! multiply by jacobians
            * lpr_jacobian (1:li_npts)

            lr_contribution = SUM ( lpr_contribution (1:li_npts) )

            apr_Projection_elt(ai_field, li_b) = apr_Projection_elt(ai_field, li_b) &
            + lr_contribution

            if ((li_cur_b == li_cur_nen) .AND. (li_b /= lp_con % opi_nen (ai_id, ai_elt)) ) then
                li_cur_b = 0
                !> \todo this works only for V = V1xV2
                li_cur_sp = li_cur_sp + 1
                li_cur_sp_id = ao_FEM % opi_InfoComposedSpace(li_sp, li_cur_sp)
                li_cur_nen = ao_FEM % opo_spaces (li_cur_sp_id) % oo_con % opi_nen (ai_id, ai_elt)
                lpr_v (:) = 0.0_wp ; lpr_v (li_cur_sp) = 1.0_wp
            end if
            li_cur_b = li_cur_b + 1
        end do

#ifdef _DEBUG        
        call printlog("build_l2_projector_vect_Local : End", ai_dtllevel = mi_dtllevel_base + 1)
#endif
    end subroutine build_l2_projector_vect_Local
    !---------------------------------------------------------------
    subroutine af_Basis(ao_ASS, ao_FEM, ai_space, ai_ndof, ai_grids_id, ai_id, ai_elt, ai_b &
        , ai_cur_sp, ai_cur_b, ai_cur_nen, apr_Basis)
        implicit none
        TYPE(ASSEMBLY) :: ao_ASS
        TYPE(FEM) :: ao_FEM
        INTEGER :: ai_space
        INTEGER :: ai_ndof
        INTEGER :: ai_grids_id
        INTEGER :: ai_id
        INTEGER :: ai_elt
        INTEGER :: ai_b
        INTEGER :: ai_cur_sp
        INTEGER :: ai_cur_b
        INTEGER :: ai_cur_nen
        real(wp), DIMENSION(:,:) :: apr_Basis
        ! LOCAL
        TYPE(BBOX), POINTER :: lp_bbox
        TYPE(GRID_DATA), pointer :: lp_grid
        TYPE(GEOMETRIES), pointer :: lp_geos
        INTEGER :: li_d
        INTEGER :: li_npts

#ifdef _DEBUG
        call printlog("af_Basis : Begin", ai_dtllevel = mi_dtllevel_base + 3)
#endif        

        !\todo prblm du id => id+1
        lp_grid => ao_FEM % opo_grids (ai_grids_id) % opo_grid ( ai_id-1 )

        li_npts = lp_grid % opo_elts ( ai_elt ) % oi_npts
!        print *, 'current sp :', ai_cur_sp
        lp_bbox => ao_ASS % opo_bbox_sp ( ai_cur_sp )
        lp_geos => ao_FEM % opo_spaces ( ai_space ) % oo_mapping

        !> \todo traiter le cas ou ndof > 1
        ! apr_Basis (1:li_ndof,1:li_npts)
        !> \todo on doit enlever 1 a cause des indices
        !> bbox % dBasis commence a 0
        apr_Basis (1,:) = lp_bbox % opr_dBasis (0,ai_cur_b-1,:)

#ifdef _DEBUG
        call printlog("af_Basis : End", ai_dtllevel = mi_dtllevel_base + 3)
#endif
    end subroutine af_Basis

end module l2_projector_vect_module


