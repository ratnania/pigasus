!     
! File:   fields_operators.F90
! Author: root
!
! Created on January 16, 2012, 9:27 AM
!

module fields_operators_module
    use tracelog_module
    use fem_def
    use assembly_def
    
    use identity_field_module
    use grad_field_module
    use curl_field_module
    use second_deriv_field_module
    use grad_s_field_module
    use second_deriv_s_field_module
    use hessian_field_module
    implicit none

    ! *******************************************
    !               LOCAL BUILDER
    ! *******************************************
    interface build_field_operators_local
            module procedure build_field_operators_local
    end interface
    interface Assembly_Field_Operator_from_elt
            module procedure Assembly_Field_Operator_from_elt
    end interface
    ! *******************************************
    
    ! *******************************************
    !               Identity_Field OPERATOR
    ! *******************************************
    interface build_Identity_Field_Local
            module procedure build_Identity_Field_Local
    end interface
    ! *******************************************

    ! *******************************************
    !               Grad_Field OPERATOR
    ! *******************************************
    interface build_Grad_Field_Local
            module procedure build_Grad_Field_Local
    end interface
    ! *******************************************

    ! *******************************************
    !               Curl_Field OPERATOR
    ! *******************************************
    interface build_Curl_Field_Local
            module procedure build_Curl_Field_Local
    end interface
    ! *******************************************

    ! *******************************************
    !               second_deriv_Field OPERATOR
    ! *******************************************
    interface build_second_deriv_field_Local
            module procedure build_second_deriv_field_Local
    end interface
    ! *******************************************

    ! *******************************************
    !               Grad_s_Field OPERATOR
    ! *******************************************
    interface build_Grad_s_Field_Local
            module procedure build_Grad_s_Field_Local
    end interface
    ! *******************************************

    ! *******************************************
    !               second_deriv_s_Field OPERATOR
    ! *******************************************
    interface build_second_deriv_s_field_Local
            module procedure build_second_deriv_s_field_Local
    end interface
    ! *******************************************

    ! *******************************************
    !               hessian_Field OPERATOR
    ! *******************************************
    interface hessian_field_Local
            module procedure build_hessian_field_Local
    end interface
    ! *******************************************
    
contains
    SUBROUTINE build_field_operators_local(self, ao_FEM, ai_grids_id, ai_id, ai_elt &
                    , ai_field)
        IMPLICIT NONE
        TYPE(ASSEMBLY) :: self
        TYPE(FEM) :: ao_FEM
        INTEGER :: ai_grids_id
        INTEGER :: ai_id
        INTEGER :: ai_elt
        INTEGER :: ai_field

    ! ********************************************************
    !   ASSEMBLING LOCAL FIELD OPERATOR
    ! ********************************************************
        select case ( ao_FEM % opi_InfoField(ai_field, INFOFIELD_OPERATOR))
            ! ********************************************************
            ! IDENTITY CASE
            ! ********************************************************
            case ( IDENTITY )
                ! COMPUTING THE LOCAL PROJECTION
                !> \todo problem avec le id=>id+1
                call build_identity_field_Local(self, ao_FEM, ai_grids_id, ai_id+1, ai_elt &
                , ai_field, self % opr_fieldh_elt)
            ! ********************************************************

            ! ********************************************************
            ! GRAD CASE
            ! ********************************************************
            case ( GRAD )
                ! COMPUTING THE LOCAL PROJECTION
                !> \todo problem avec le id=>id+1
                call build_grad_field_Local(self, ao_FEM, ai_grids_id, ai_id+1, ai_elt &
                , ai_field, self % opr_fieldh_elt)
            ! ********************************************************

            ! ********************************************************
            ! CURL CASE
            ! ********************************************************
            case ( CURL )
                ! COMPUTING THE LOCAL PROJECTION
                !> \todo problem avec le id=>id+1
                call build_curl_field_Local(self, ao_FEM, ai_grids_id, ai_id+1, ai_elt &
                , ai_field, self % opr_fieldh_elt)
            ! ********************************************************

            ! ********************************************************
            ! SECOND_DERIV CASE
            ! ********************************************************
            case ( SECOND_DERIV_FIELD )
                ! COMPUTING THE LOCAL PROJECTION
                !> \todo problem avec le id=>id+1
                call build_second_deriv_field_Local(self, ao_FEM, ai_grids_id, ai_id+1, ai_elt &
                , ai_field, self % opr_fieldh_elt)
            ! ********************************************************

            ! ********************************************************
            ! GRAD_S CASE
            ! ********************************************************
            case ( GRAD_S )
                ! COMPUTING THE LOCAL PROJECTION
                !> \todo problem avec le id=>id+1
                call build_grad_s_field_Local(self, ao_FEM, ai_grids_id, ai_id+1, ai_elt &
                , ai_field, self % opr_fieldh_elt)
            ! ********************************************************

            ! ********************************************************
            ! SECOND_DERIV_S CASE
            ! ********************************************************
            case ( SECOND_DERIV_S_FIELD )
                ! COMPUTING THE LOCAL PROJECTION
                !> \todo problem avec le id=>id+1
                call build_second_deriv_s_field_Local(self, ao_FEM, ai_grids_id, ai_id+1, ai_elt &
                , ai_field, self % opr_fieldh_elt)
            ! ********************************************************

            ! ********************************************************
            ! HESSIAN CASE
            ! ********************************************************
            case ( HESSIAN_FIELD )
                ! COMPUTING THE LOCAL PROJECTION
                !> \todo problem avec le id=>id+1
                call build_hessian_field_Local(self, ao_FEM, ai_grids_id, ai_id+1, ai_elt &
                , ai_field, self % opr_fieldh_elt)
            ! ********************************************************
            

        case Default
            call printlog("Type Field-Operator Not Yet implemented"    &
            , ai_dtllevel = 0)

        end select
        ! ********************************************************
    END SUBROUTINE build_field_operators_local
!----------------------------------------------------------------------------------------------
    subroutine Assembly_Field_Operator_from_elt(self, ao_FEM, ai_id, ai_elt, ai_field)
        implicit none
        type(ASSEMBLY) :: self
        type(FEM) :: ao_FEM
        integer :: ai_id
        integer :: ai_elt
        integer :: ai_field
        ! LOCAL VARIABLES
        integer :: li_ndof

        li_ndof = ao_FEM % opi_InfoField(ai_field, INFOFIELD_NDOF)

        ao_FEM % opo_F (ai_field) % opr_val(1:li_ndof,:,ai_elt,:) = self % opr_fieldh_elt(ai_field,1:li_ndof,:,:)

#ifdef _DEBUG
        print *, 'ao_FEM % opo_F (ai_field) % opr_val =', ao_FEM % opo_F (ai_field) % opr_val(1:li_ndof,:,ai_elt,:)
#endif

    end subroutine Assembly_Field_Operator_from_elt

end module fields_operators_module

