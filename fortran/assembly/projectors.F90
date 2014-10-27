!     
! File:   projectors.F90
! Author: root
!
! Created on January 27, 2012, 8:33 AM
!

module projectors_module
    use tracelog_module
    use fem_def
    use assembly_def
    
    use l2_projector_module
    use l2_projector_grad_module
    use l2_projector_curl_module
!    use l2_projector_vect_module
    implicit none

    interface build_field_projectors_Local
            module procedure build_field_projectors_Local
    end interface
    interface Assembly_Projections_from_elt
            module procedure Assembly_Projections_from_elt
    end interface

    interface build_l2_projector_Local
            module procedure build_l2_projector_Local
    end interface

    interface build_l2_projector_grad_Local
            module procedure build_l2_projector_grad_Local
    end interface

    interface build_l2_projector_curl_Local
            module procedure build_l2_projector_curl_Local
    end interface

!    interface build_l2_projector_vect_Local
!            module procedure build_l2_projector_vect_Local
!    end interface

contains
SUBROUTINE  build_field_projectors_Local(self, ao_FEM, ai_grids_id, ai_id, ai_elt &
                    , ai_field)
        IMPLICIT NONE
        TYPE(ASSEMBLY) :: self
        TYPE(FEM) :: ao_FEM
        INTEGER :: ai_grids_id
        INTEGER :: ai_id
        INTEGER :: ai_elt
        INTEGER :: ai_field
    
    select case ( ao_FEM % opi_InfoField(ai_field, INFOFIELD_OPERATOR))
        ! ********************************************************
        ! IDENTITY CASE
        ! ********************************************************
        case ( IDENTITY )
            ! COMPUTING THE LOCAL PROJECTION
            !> \todo problem avec le id=>id+1
            call build_l2_projector_Local(self, ao_FEM, ai_grids_id, ai_id+1, ai_elt &
            , ai_field, self % opr_Projection_elt)
        ! ********************************************************

        ! ********************************************************
        ! GRAD CASE
        ! ********************************************************
        case ( GRAD )
            ! COMPUTING THE LOCAL PROJECTION
            !> \todo problem avec le id=>id+1
            call build_l2_projector_grad_Local(self, ao_FEM, ai_grids_id, ai_id+1, ai_elt &
            , ai_field, self % opr_Projection_elt)
        ! ********************************************************

        ! ********************************************************
        ! CURL CASE
        ! ********************************************************
        case ( CURL )
            ! COMPUTING THE LOCAL PROJECTION
            !> \todo problem avec le id=>id+1
            call build_l2_projector_curl_Local(self, ao_FEM, ai_grids_id, ai_id+1, ai_elt &
            , ai_field, self % opr_Projection_elt)
        ! ********************************************************

!        ! ********************************************************
!        ! IDENTITY_VECT CASE
!        ! ********************************************************
!        case ( IDENTITY_VECT )
!            ! COMPUTING THE LOCAL PROJECTION
!            !> \todo problem avec le id=>id+1
!            call build_l2_projector_vect_Local(self, ao_FEM, ai_grids_id, ai_id+1, ai_elt &
!            , ai_field, self % opr_Projection_elt)
!        ! ********************************************************

        case Default

            call printlog("Type Projection-Operator Not Yet implemented"    &
            , ai_dtllevel = 0)

    end select


END SUBROUTINE  build_field_projectors_Local
!----------------------------------------------------------------------------------------------
    subroutine Assembly_Projections_from_elt(self, ao_FEM, ai_id, ai_elt, ai_field)
        implicit none
        type(ASSEMBLY) :: self
        type(FEM) :: ao_FEM
        integer :: ai_id
        integer :: ai_elt
        integer :: ai_field
        ! LOCAL VARIABLES
        integer :: li_b
        integer :: li_A
        integer :: li_sp
        integer :: li_conelt
        TYPE(CONNECTIVITY), pointer :: lp_con

        li_sp   = ao_FEM % opi_InfoField(ai_field, INFOFIELD_SPACE )
        lp_con  => ao_FEM % opo_spaces (li_sp) % oo_con

        li_conelt = lp_con % opi_real_elts (ai_id-1,ai_elt)

        do li_b = 1, lp_con % opi_nen (ai_id)

            li_A = lp_con % opi_LM ( ai_id, li_b, li_conelt )

            if (li_A == 0) then
                cycle
            end if

            ao_FEM % opo_F (ai_field) % opr_c (li_A) = ao_FEM % opo_F (ai_field) % opr_c (li_A) &
            + self % opr_Projection_elt (ai_field, li_b)

        end do

!#ifdef _DEBUG
!        print *, 'ao_FEM % opo_F (ai_field) % opr_c (:)=', ao_FEM % opo_F (ai_field) % opr_c (:)
!#endif

    end subroutine Assembly_Projections_from_elt
end module projectors_module

