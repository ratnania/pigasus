!     
! File:   norms.F90
! Author: root
!
! Created on February 14, 2012, 11:53 AM
!

module norms_module
    use tracelog_module
    use fem_def
    use assembly_def
    
    use l2_norm_module
    use H1_norm_module
    implicit none

    interface build_norm_operators_local
            module procedure build_norm_operators_local
    end interface
    interface Assembly_Norms_from_elt
            module procedure Assembly_Norms_from_elt
    end interface

    interface build_l2_norm_Local
            module procedure build_l2_norm_Local
    end interface

    interface build_H1_norm_Local
            module procedure build_H1_norm_Local
    end interface
    
contains
    SUBROUTINE build_norm_operators_local(self, ao_FEM, ai_grids_id, ai_id, ai_elt &
                    , ai_norm)
        IMPLICIT NONE
        TYPE(ASSEMBLY) :: self
        TYPE(FEM) :: ao_FEM
        INTEGER :: ai_grids_id
        INTEGER :: ai_id
        INTEGER :: ai_elt
        INTEGER :: ai_norm
        
        select case ( ao_FEM % opi_InfoNorm(ai_norm, INFONORM_TYPE))

            ! ********************************************************
            ! L2-NORM CASE
            ! ********************************************************
            case ( NORM_L2 )
                ! COMPUTING THE LOCAL NORM
                !> \todo problem avec le id=>id+1
                call build_l2_norm_Local(self, ao_FEM, ai_grids_id, ai_id+1, ai_elt &
                , ai_norm, self % opr_Norm_elt)
            ! ********************************************************

            ! ********************************************************
            ! H1-NORM CASE
            ! ********************************************************
            case ( NORM_H1 )
                ! COMPUTING THE LOCAL NORM
                !> \todo problem avec le id=>id+1
                call build_H1_norm_Local(self, ao_FEM, ai_grids_id, ai_id+1, ai_elt &
                , ai_norm, self % opr_Norm_elt)
            ! ********************************************************

        case Default

            call printlog("Type Norm Not Yet implemented", ai_dtllevel = 0)

        end select

    END SUBROUTINE build_norm_operators_local
!----------------------------------------------------------------------------------------------
    subroutine Assembly_Norms_from_elt(self, ao_FEM, ai_id, ai_elt, ai_norm)
        implicit none
        type(ASSEMBLY) :: self
        type(FEM) :: ao_FEM
        integer :: ai_id
        integer :: ai_elt
        integer :: ai_norm
        ! LOCAL VARIABLES
        integer :: li_b
        integer :: li_A
        integer :: li_sp
        TYPE(CONNECTIVITY), pointer :: lp_con

        ao_FEM % opo_N (ai_norm) % opr_values (ai_id,ai_elt) = self % opr_Norm_elt(ai_norm)
!        print *, ai_elt  , self % opr_Norm_elt(ai_norm)

    end subroutine Assembly_Norms_from_elt
end module norms_module

