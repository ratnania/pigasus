!     
! File:   identity_field.F90
! Author: root
!
! Created on March 8, 2012, 10:08 AM
!

module identity_field_module
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

    public :: build_identity_field_Local

#ifdef _DEBUG
    integer, parameter, private  :: mi_dtllevel_base = 0
#else
    integer, parameter, private  :: mi_dtllevel_base = 2
#endif

contains

    !----------------------------------------------------------------------------------------------
    subroutine build_identity_field_Local(ao_ASS, ao_FEM, ai_grids_id, ai_id, ai_elt, ai_field, apr_fieldh)
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
        integer :: li_A
        integer :: li_sp
        integer :: li_bloc
        integer :: li_ndof
        integer :: li_map
        INTEGER :: li_loc_id
        INTEGER :: li_nparam
        INTEGER :: li_conelt
        INTEGER :: li_field_id
        real(wp) :: lr_contribution
        INTEGER, PARAMETER :: MAXDIM = 3
        real(wp), dimension(ao_ASS % oi_maxndof,ao_ASS % oi_maxnpts) :: lpr_B
        real(wp), dimension(ao_ASS % oi_maxnpts) :: lpr_jacobian
        real(wp), DIMENSION(ao_ASS % oi_maxnpts) :: lpr_contribution
        real(wp), dimension(ao_ASS % oi_maxnpts) :: lpr_w
        TYPE(CONNECTIVITY), pointer :: lp_con
        TYPE(CONNECTIVITY), pointer :: lp_conprime
        TYPE(GRID_DATA), pointer :: lp_grid
        TYPE(PHYSICAL_BASIS), pointer :: lp_pBasis

#ifdef _DEBUG        
        call printlog("build_identity_field_Local : Begin", ai_dtllevel = mi_dtllevel_base + 1)
#endif
!print *, 'ai_elt=', ai_elt

        lpr_B = 0.0_wp

        !\todo prblm du id => id+1
        lp_grid => ao_FEM % opo_grids ( ai_grids_id ) % opo_grid ( ai_id-1 )
        li_npts = lp_grid % opo_elts ( ai_elt ) % oi_npts

        li_sp       = ao_FEM % opi_InfoField(ai_field, INFOFIELD_SPACE )
        li_ndof     = ao_FEM % opi_InfoField(ai_field, INFOFIELD_NDOF )
        li_loc_id   = ao_FEM % opi_InfoField(ai_field, INFOFIELD_LOCID )
!        li_nparam   = ao_FEM % opi_InfoField(ai_field, INFOFIELD_NPARAM )
        li_field_id = ao_FEM % opi_InfoField(ai_field, INFOFIELD_OPERANDE )

        lp_con          => ao_FEM % opo_spaces (li_sp) % oo_con
        lp_pBasis       => ao_ASS % opo_pBasis (li_sp)

        li_conelt = lp_con % opi_real_elts (ai_id-1,ai_elt)
        !> \todo attention +1
        do li_b = 1, lp_con % opi_nen(ai_id)

            li_A = lp_con % opi_LM(ai_id, li_b, li_conelt)
            
            if (li_A == 0) then
                cycle
            end if

            lpr_B (1,1:li_npts) = lp_pBasis % opr_B (li_b, 1:li_npts)

            !> \todo traiter le cas ou ndof > 1
            ! apr_Basis (1:li_ndof,1:li_npts)
            !> \todo on doit enlever 1 a cause des indices
            !> bbox % dBasis commence a 0

            apr_fieldh(ai_field,1,1, 1:li_npts) = apr_fieldh(ai_field,1,1, 1:li_npts) &
            + ao_FEM % opo_F(li_field_id) % opr_c(li_A) &
            * lpr_B (1,:)

        end do

#ifdef _DEBUG        
        call printlog("build_identity_field_Local : End", ai_dtllevel = mi_dtllevel_base + 1)
#endif
    end subroutine build_identity_field_Local

end module identity_field_module

