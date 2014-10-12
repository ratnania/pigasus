!     
! File:   l2_projector.F90
! Author: root
!
! Created on January 27, 2012, 8:33 AM
!


MODULE l2_projector_module
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
    SUBROUTINE build_l2_projector_Local(ao_ASS, ao_FEM, ai_grids_id, ai_id, ai_elt, ai_color, apr_Projection_elt)
        implicit none
        TYPE(ASSEMBLY) :: ao_ASS
        type(FEM) :: ao_FEM
        !> grids id
        integer :: ai_grids_id
        !> patch id
        integer :: ai_id
        !> current element
        integer :: ai_elt
        !> current color 
        integer :: ai_color
        real(wp), dimension(0:,:) :: apr_Projection_elt
        ! LOCAL VARIABLES
        integer :: li_f_ref
        integer :: li_nfields
        integer :: li_field
        integer :: li_npts
        integer :: li_b
        integer :: li_i
        integer :: li_A
        integer :: li_sp
        integer :: li_bloc
        integer :: li_ndof
        integer :: li_map
        INTEGER :: li_loc_id
        INTEGER :: li_nen
        INTEGER :: li_dim
        real(wp) :: lr_contribution
        real(wp), dimension (ao_ASS % oi_maxndof,ao_ASS % oi_maxnpts) :: lpr_B
        real(wp), dimension(ao_ASS % oi_maxnpts) :: lpr_jacobian
        real(wp), dimension(ao_ASS % oi_maxnpts) :: lpr_w
        real(wp), DIMENSION(ao_FEM % oi_maxcoloraddto, ao_ASS % oi_maxnpts) :: lpr_contribution
        real(wp), dimension(ao_FEM % oi_maxcoloraddto, ao_ASS % oi_maxnpts) :: lpr_f
        REAL(WP), dimension(ao_FEM % oi_maxcoloraddto, ao_FEM % oi_maxnen) :: lpr_Projection
        TYPE(CONNECTIVITY), pointer :: lp_con
        TYPE(CONNECTIVITY), pointer :: lp_conprime
        TYPE(GRID_DATA), pointer :: lp_grid
        TYPE(PHYSICAL_BASIS), pointer :: lp_pBasis

#ifdef _DEBUG
        call printlog("build_l2_projector_Local : Begin", ai_dtllevel = mi_dtllevel_base + 1)

        call concatmsg("ai_elt : ", ai_dtllevel = mi_dtllevel_base + 1)
        call concatmsg(ai_elt     , ai_dtllevel = mi_dtllevel_base + 1)
        call printmsg(              ai_dtllevel = mi_dtllevel_base + 1)
#endif

        lpr_B           = 0.0_wp
        lpr_Projection  = 0.0_wp
        lpr_f           = 0.0_wp

        li_nfields  = ao_FEM % opo_colors(ai_color) % opi_objects_toassembly (0)

        ! ...
        ! commun to all fields of the same color
        ! ...
        ! get field-id from color
        li_f_ref = 0
        li_field = ao_FEM % opo_colors(ai_color) % opi_objects_toassembly (li_f_ref)

        li_dim = ao_FEM % opi_dim(ai_grids_id)

        !\todo prblm du id => id+1
        lp_grid => ao_FEM % opo_grids ( ai_grids_id ) % opo_grid ( ai_id-1 )
        li_npts = lp_grid % opo_elts ( ai_elt ) % oi_npts

        li_sp       = ao_FEM % opi_InfoField(li_field, INFOFIELD_SPACE )
        li_ndof     = ao_FEM % opi_InfoField(li_field, INFOFIELD_NDOF )

        lp_con          => ao_FEM % opo_spaces (li_sp) % oo_con
        lp_pBasis       => ao_ASS % opo_pBasis (li_sp)
        
        lpr_jacobian = ao_ASS % opo_info_gr (ai_grids_id) % opr_jacobians
        lpr_jacobian = dabs(lpr_jacobian)

        li_nen          = lp_con % opi_nen(ai_id)

        lpr_w (1 : li_npts)        = lp_grid % opo_elts ( ai_elt ) % opr_w (1:li_npts)

        DO li_f_ref = 1, li_nfields
            ! get field-id from color
            li_field    = ao_FEM % opo_colors(ai_color) % opi_objects_toassembly (li_f_ref)
            li_loc_id   = ao_FEM % opi_InfoField(li_field, INFOFIELD_LOCID )

            lpr_f (li_f_ref, 1:li_npts) = lp_grid % opo_elts (ai_elt) % opr_values_f(li_loc_id, 1, 1 : li_npts)
        END DO

#ifdef _PERF        
        CALL CPU_TIME ( time_begin )
#endif

        DO li_b = 1, li_nen

            lpr_B (1,1:li_npts) = lp_pBasis % opr_B (li_b, 1:li_npts)

            lpr_contribution = 0.0_wp
            DO li_f_ref = 1, li_nfields
                ! multiply by input data
                lpr_contribution (li_f_ref,1:li_npts) = lpr_f(li_f_ref,1 : li_npts)   &
                ! multiply by the left basis
                * lpr_B (1,1:li_npts)      
            END DO

            DO li_f_ref = 1, li_nfields
                ! multiply by input data
                lpr_contribution (li_f_ref,1:li_npts) = lpr_contribution (li_f_ref,1:li_npts) &              
                ! multiply by quadratures weights
                * lpr_w (1:li_npts) &
                ! multiply by jacobians
                * lpr_jacobian (1:li_npts)

                lr_contribution = SUM ( lpr_contribution (li_f_ref,1:li_npts) )

                lpr_Projection(li_f_ref,li_b) = lpr_Projection(li_f_ref,li_b) &
                + lr_contribution
            END DO

        END DO

#ifdef _PERF        
        CALL CPU_TIME ( time_end )
        PRINT *, 'Time of operation was ', &
        time_end - time_begin, ' seconds'
#endif

#ifdef _DEBUG
        call concatmsg("li_nfields : ", ai_dtllevel = mi_dtllevel_base + 1)
        call concatmsg(li_nfields     , ai_dtllevel = mi_dtllevel_base + 1)
        call printmsg(                  ai_dtllevel = mi_dtllevel_base + 1)
#endif

        DO li_f_ref = 1, li_nfields
            ! get field-id from color
            li_field = ao_FEM % opo_colors(ai_color) % opi_objects_toassembly (li_f_ref)
            apr_Projection_elt(li_field, 1:li_nen) = lpr_Projection (li_f_ref, 1:li_nen)

#ifdef _DEBUG
        call concatmsg("+++ field : ", ai_dtllevel = mi_dtllevel_base + 1)
        call concatmsg(li_field      , ai_dtllevel = mi_dtllevel_base + 1)
        call printmsg(                 ai_dtllevel = mi_dtllevel_base + 1)
#endif

        END DO
        
#ifdef _DEBUG
        call printlog("build_l2_projector_Local : End", ai_dtllevel = mi_dtllevel_base + 1)
#endif
    END SUBROUTINE build_l2_projector_Local

END MODULE l2_projector_module

