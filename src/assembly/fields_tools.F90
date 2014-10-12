!     
! File:   fields_tools.F90
! Author: abdelkaderratnani
!
! Created on February 4, 2012, 9:38 AM
!

module fields_tools_module
    use used_precision
    use tracelog_module
    use fem_def
    use fem_module
    use assembly_def
    use assembly_module
    use fields_operators_module
    implicit none

#ifdef _DEBUG
    integer, parameter, private  :: mi_dtllevel_base = 0
#else
    integer, parameter, private  :: mi_dtllevel_base = 2
#endif

contains
    !----------------------------------------------------------------------------------------------
    subroutine assembly_field_on_grids(ao_ASS, ao_FEM, ai_field, ai_patch, apr_fieldh &
        , ai_ndof, ai_maxder, ai_nel, ai_maxnpts)
        implicit none
        type(ASSEMBLY) :: ao_ASS
        type(FEM) :: ao_FEM
        integer, intent(in) :: ai_field
        integer, intent(in) :: ai_patch
        integer, intent(in) :: ai_ndof
        integer, intent(in) :: ai_maxder
        integer, intent(in) :: ai_nel
        integer, intent(in) :: ai_maxnpts
        real*8, dimension(ai_ndof, ai_maxder, ai_nel, ai_maxnpts), intent(out) :: apr_fieldh
        ! LOCAL VARIABLES
        integer :: li_b
        integer :: li_A
        integer :: li_space
        integer :: li_elt
        integer :: li_nel
        integer :: li_patch
        integer :: li_id
        integer :: li_map
        integer :: li_grids
        integer :: li_npts
        integer :: li_realelt
        type(GRID_DATA), pointer :: lp_grid
        type(BBOX), pointer :: lp_bbox
        type(MAPPING), pointer :: lp_mapping
        type(SPACE), pointer :: lp_space
        TYPE(CONNECTIVITY), pointer :: lp_con
        INTEGER, PARAMETER :: NFIELD = 2
        real*8, dimension(NFIELD,ai_ndof,ai_maxder, ai_maxnpts) :: lpr_fieldh

        CALL printlog("assembly_field_on_grids : Start", ai_dtllevel = 1)

        print *, "assembly_field_on_grids"
        STOp

        li_space = ao_FEM % opi_InfoField(ai_field, INFOFIELD_SPACE)
        li_grids = ao_FEM % opi_InfoSpace(li_space, INFOSPACE_GRIDS)

        lp_bbox => ao_ASS % opo_bbox_sp(li_space)
        lp_space => ao_FEM % opo_spaces(li_space)
        lp_con => lp_space % oo_con

        ! LOOP THROUG PATCHS
        li_patch = ai_patch
        li_id = li_patch + 1

        li_nel = ao_FEM % opi_InfoPatch(li_grids, li_patch, INFOPATCH_NEL)

        lp_grid => ao_FEM % opo_grids(li_grids) % opo_grid(li_patch)

        apr_fieldh = 0.0_wp

        DO li_elt = 1, li_nel
            
            lpr_fieldh = 0.0_wp
            
            li_npts = lp_grid % opo_elts(li_elt) % oi_npts
            li_realelt = lp_con % opi_real_elts (li_patch,li_elt)

            !> \todo a changer: li_patch+1 en li_patch, ds la lib grids
!            CALL assembly_logical_basis(lp_grid, lp_space % oo_mapping &
!            , ao_ASS % oi_matrix_nderiv, li_id, li_elt, li_realelt &
!            , ao_ASS % opo_bbox_sp(li_space) &
!            , ai_nen = lp_space % oi_maxnen &
!            , ai_ptw_evaluation = 0)

            ! ********************************************************
            !                   METRIC TREATMENT
            ! ********************************************************
            CALL update_info_metric( ao_ASS, ao_FEM, li_space, li_patch, li_elt, li_npts )
            ! ********************************************************

            ! ********************************************************
            !                   SPACES TREATMENT
            ! COMPUTING ALL DERIVATIVES OF ALL BASIS FUNCTIONS
            ! ********************************************************
            CALL update_info_space( ao_ASS, ao_FEM, li_space, li_patch, li_elt, li_npts )
            ! ********************************************************
            
            select case ( ao_FEM % opi_InfoField(ai_field, INFOFIELD_OPERATOR))
                ! ********************************************************
                ! IDENTITY CASE
                ! ********************************************************
                case ( IDENTITY )
                    ! COMPUTING THE LOCAL PROJECTION
                    !> \todo problem avec le id=>id+1
                    call build_identity_field_Local(ao_ASS, ao_FEM, li_grids, li_id, li_elt &
                    , ai_field, lpr_fieldh)
                ! ********************************************************

                ! ********************************************************
                ! GRAD CASE
                ! ********************************************************
                case ( GRAD )
                    ! COMPUTING THE LOCAL PROJECTION
                    !> \todo problem avec le id=>id+1
                    call build_grad_field_Local(ao_ASS, ao_FEM, li_grids, li_id, li_elt &
                    , ai_field, lpr_fieldh)
                ! ********************************************************

                ! ********************************************************
                ! CURL CASE
                ! ********************************************************
                case ( CURL )
                    ! COMPUTING THE LOCAL PROJECTION
                    !> \todo problem avec le id=>id+1
                    call build_curl_field_Local(ao_ASS, ao_FEM, li_grids, li_id, li_elt &
                    , ai_field, lpr_fieldh)

            case Default

                call printlog("Type Projection-Operator Not Yet implemented"    &
                , ai_dtllevel = mi_dtllevel_base +  0)

            end select

            apr_fieldh(:,:,li_elt,:) = lpr_fieldh(1,:,:,:)

        END DO

        CALL printlog("assembly_field_on_grids  : End", ai_dtllevel = 1)

    end subroutine assembly_field_on_grids
    !----------------------------------------------------------------------------------------------
    subroutine assembly_field_on_grids_2D(ao_ASS, ao_FEM, ai_field, ai_patch, apr_fieldh &
        , ai_ndof, ai_nel, ai_maxdirnpts)
        implicit none
        type(ASSEMBLY) :: ao_ASS
        type(FEM) :: ao_FEM
        integer, intent(in) :: ai_field
        integer, intent(in) :: ai_patch
        integer, intent(in) :: ai_ndof
        integer, intent(in) :: ai_nel
        integer, intent(in) :: ai_maxdirnpts
        real*8, dimension(ai_ndof, ai_nel, ai_maxdirnpts, ai_maxdirnpts), intent(out) :: apr_fieldh
        ! LOCAL VARIABLES
        integer :: li_b
        integer :: li_A
        integer :: li_space
        integer :: li_elt
        integer :: li_nel
        integer :: li_patch
        integer :: li_id
        integer :: li_map
        integer :: li_grids
        integer :: li_npts
        integer :: li_pt
        integer :: li_realelt
        INTEGER, DIMENSION(:), pointer :: lpi_pti
        INTEGER, DIMENSION(:), pointer :: lpi_npts
        type(GRID_DATA), pointer :: lp_grid
        type(BBOX), pointer :: lp_bbox
        type(MAPPING), pointer :: lp_mapping
        type(SPACE), pointer :: lp_space
        TYPE(CONNECTIVITY), pointer :: lp_con

        CALL printlog("assembly_field_on_grids_2D : Start", ai_dtllevel = 1)

        li_space = ao_FEM % opi_InfoField(ai_field, INFOFIELD_SPACE)
        li_grids = ao_FEM % opi_InfoSpace(li_space, INFOSPACE_GRIDS)

        lp_bbox => ao_ASS % opo_bbox_sp(li_space)
        lp_space => ao_FEM % opo_spaces(li_space)
        lp_con => lp_space % oo_con

        ! LOOP THROUG PATCHS
        li_patch = ai_patch
        li_id = li_patch + 1

        li_nel = ao_FEM % opi_InfoPatch(li_grids, li_patch, INFOPATCH_NEL)

        lp_grid => ao_FEM % opo_grids(li_grids) % opo_grid(li_patch)

        ALLOCATE(lpi_pti(ao_FEM % opi_dim(li_grids)))
        ALLOCATE(lpi_npts(ao_FEM % opi_dim(li_grids)))

        apr_fieldh = 0.0_wp

        DO li_elt = 1, li_nel
            li_npts = lp_grid % opo_elts(li_elt) % oi_npts
            lpi_npts = lp_grid % opo_elts(li_elt) % opi_npts

            li_realelt = lp_con % opi_real_elts (li_patch,li_elt)

            !> \todo a changer: li_patch+1 en li_patch, ds la lib grids
            CALL assembly_logical_basis(lp_grid, lp_space % oo_mapping &
            , ao_ASS % oi_matrix_nderiv, li_id, li_elt, li_realelt &
            , ao_ASS % opo_bbox_sp(li_space) &
!            , ai_nen = lp_space % oi_maxnen &
            , lp_space % oi_maxnen &
!            , ai_ptw_evaluation = 0)
            , 0)

            !> \todo attention +1
            do li_b = 1, lp_con % opi_nen(li_id)

                li_A = lp_con % opi_LM(li_id, li_b, li_elt)

                if (li_A == 0) then
                    cycle
                end if

                !> \todo traiter le cas ou ndof > 1
                ! apr_Basis (1:li_ndof,1:li_npts)
                !> \todo on doit enlever 1 a cause des indices
                !> bbox % dBasis commence a 0

                DO li_pt = 0, li_npts - 1

                    CALL get_indices(lpi_npts, li_pt, lpi_pti)
                    lpi_pti = lpi_pti + 1
                    apr_fieldh(1, li_elt, lpi_pti(1), lpi_pti(2)) = apr_fieldh(1, li_elt, lpi_pti(1), lpi_pti(2)) &
                    +ao_FEM % opo_F(ai_field) % opr_c(li_A) &
                    *lp_bbox % opr_dBasis(0, li_b - 1, li_pt)

                END DO

            end do

        END DO

        DEALLOCATE(lpi_pti)
        DEALLOCATE(lpi_npts)

        CALL printlog("assembly_field_on_grids_2D : End", ai_dtllevel = 1)

    end subroutine assembly_field_on_grids_2D
    !----------------------------------------------------------------------------------------------
end module fields_tools_module

