!     
! File:   operators.F90
! Author: root
!
! Created on January 16, 2012, 9:27 AM
!

module operators_module
    use tracelog_module
    use fem_def
    use assembly_def
    
    use mass_module
    use stiffness_module
    use advection_module
    use second_deriv_module

!    use mass_t2_module
!    use mass_vect_module
!    use rotational_module
!    use rotational_scalar_module
    implicit none

    ! *******************************************
    !               MATRIX BUILDER
    ! *******************************************
    interface build_matrix_operators_Local
            module procedure build_matrix_operators_Local
    end interface
    interface Assembly_Matrices_from_elt
            module procedure Assembly_Matrices_from_elt
    end interface
    ! *******************************************

contains

    SUBROUTINE build_matrix_operators_Local(self, ao_FEM, ai_grids_id, ai_id, ai_elt &
                    , ai_operator)
        IMPLICIT NONE
        TYPE(ASSEMBLY) :: self
        TYPE(FEM) :: ao_FEM
        INTEGER :: ai_grids_id
        INTEGER :: ai_id
        INTEGER :: ai_elt
        INTEGER :: ai_operator
        ! LOCAL
        INTEGER :: li_tensorlevel

        li_tensorlevel  = ao_fem % opi_InfoPatch ( ai_grids_id, ai_id, INFOPATCH_TENSOR )

        select case ( ao_FEM % opi_InfoOperator(ai_operator, INFOOPERATOR_TYPE))

            ! ********************************************************
            ! MASS MATRIX CASE
            ! ********************************************************
        case ( MASS)
            ! COMPUTING THE LOCAL MATRIX
            !> \todo problem avec le id=>id+1
            call build_Mass_Local(self, ao_FEM, ai_grids_id, ai_id+1, ai_elt &
            , ai_operator, self % opr_Matrix_elt)
!            IF (li_tensorlevel==1) THEN
!            ELSEIF (li_tensorlevel==2) THEN
!                call build_Mass_Local_t2(self, ao_FEM, ai_grids_id, ai_id+1, ai_elt, ai_operator)
!            END IF
            ! ********************************************************

            ! ********************************************************
            ! STIFFNESS MATRIX CASE
            ! ********************************************************
        case ( STIFFNESS)
            ! COMPUTING THE LOCAL MATRIX
            !> \todo problem avec le id=>id+1
            call build_Stiffness_Local(self, ao_FEM, ai_grids_id, ai_id+1, ai_elt &
            , ai_operator, self % opr_Matrix_elt)
            ! ********************************************************

            ! ********************************************************
            ! ADVECTION MATRIX CASE
            ! ********************************************************
        case ( ADVECTION)
            ! COMPUTING THE LOCAL MATRIX
            !> \todo problem avec le id=>id+1
            call build_Advection_Local(self, ao_FEM, ai_grids_id, ai_id+1, ai_elt &
            , ai_operator, self % opr_Matrix_elt)
            ! ********************************************************

            ! ********************************************************
            ! SECOND_DERIV MATRIX CASE
            ! ********************************************************
        case ( SECOND_DERIV)
            ! COMPUTING THE LOCAL MATRIX
            !> \todo problem avec le id=>id+1
            call build_second_deriv_Local(self, ao_FEM, ai_grids_id, ai_id+1, ai_elt &
            , ai_operator, self % opr_Matrix_elt)
            ! ********************************************************

!            ! ********************************************************
!            ! MASS_VECT MATRIX CASE
!            ! ********************************************************
!        case ( MASS_VECT)
!            ! COMPUTING THE LOCAL MATRIX
!            !> \todo problem avec le id=>id+1
!            call build_Mass_Vect_Local(self, ao_FEM, ai_grids_id, ai_id+1, ai_elt &
!            , ai_operator, self % opr_Matrix_elt)
!            ! ********************************************************
!
!            ! ********************************************************
!            ! ROTATIONAL MATRIX CASE
!            ! ********************************************************
!        case ( ROTATIONAL)
!            ! COMPUTING THE LOCAL MATRIX
!            !> \todo problem avec le id=>id+1
!            call build_Rotational_Local(self, ao_FEM, ai_grids_id, ai_id+1, ai_elt &
!            , ai_operator, self % opr_Matrix_elt)
!            ! ********************************************************
!
!            ! ********************************************************
!            ! ROTATIONAL_SCALAR MATRIX CASE
!            ! ********************************************************
!        case ( ROTATIONAL_SCALAR)
!            ! COMPUTING THE LOCAL MATRIX
!            !> \todo problem avec le id=>id+1
!            call build_rotational_scalar_Local(self, ao_FEM, ai_grids_id, ai_id+1, ai_elt &
!            , ai_operator, self % opr_Matrix_elt)
!            ! ********************************************************

        case Default

            call printlog("Type Matrix Not Yet implemented", ai_dtllevel = 0)
            call concatmsg("The type of the matrix with id ", ai_dtllevel = 0)
            call concatmsg(ai_operator, ai_dtllevel = 0)
            call printmsg(ai_dtllevel = 0)

        end select

    END SUBROUTINE build_matrix_operators_Local
    !----------------------------------------------------------------------------------------------
    subroutine Assembly_Matrices_from_elt(self, ao_FEM, ai_grids_id, ai_id, ai_elt, ai_operator)
        implicit none
        type(ASSEMBLY) :: self
        type(FEM) :: ao_FEM
        integer :: ai_grids_id
        integer :: ai_id
        integer :: ai_elt
        integer :: ai_operator
        ! LOCAL
        INTEGER :: li_tensorlevel

        li_tensorlevel = ao_fem % opi_InfoPatch ( ai_grids_id, ai_id, INFOPATCH_TENSOR )

        IF (li_tensorlevel==1) THEN

            CALL Assembly_Matrices_from_elt_basic(self, ao_FEM, ai_id+1, ai_elt, ai_operator)

        ELSE

!            CALL Assembly_Matrices_from_elt_t2(self, ao_FEM, ai_id+1, ai_elt, ai_operator)
                print *, 'Not used anymore'

        END IF

    end subroutine Assembly_Matrices_from_elt
    !----------------------------------------------------------------------------------------------
    subroutine Assembly_Matrices_from_elt_basic(self, ao_FEM, ai_id, ai_elt, ai_operator)
        implicit none
        type(ASSEMBLY) :: self
        type(FEM) :: ao_FEM
        integer :: ai_id
        integer :: ai_elt
        integer :: ai_operator
        ! LOCAL VARIABLES
        integer :: li_addtoMatrix
        integer :: li_b
        integer :: li_bprime
        integer :: li_A
        integer :: li_Aprime
        integer :: li_sp
        integer :: li_spprime
        integer :: li_conelt
        integer :: li_conprimeelt
        integer :: li_matrix
        integer :: li_nmatrices
        integer :: li_i
        TYPE(CONNECTIVITY), pointer :: lp_con
        TYPE(CONNECTIVITY), pointer :: lp_conprime
        TYPE(OPER), pointer :: lp_operator
        INTEGER :: ierr
        REAL(KIND=SPM_COEF_KIND) :: lr_val
        REAL(KIND=SPM_COEF_KIND) :: lr_scale

#ifdef _MURGE
        REAL(KIND=MURGE_COEF_KIND) :: val
        INTEGER(KIND=MURGE_INTS_KIND) :: id
        INTEGER(KIND=MURGE_INTS_KIND) :: ierr
#endif      

!        li_addtoMatrix = ao_FEM % opi_InfoOperator (ai_operator, INFOOPERATOR_ADDTO)

        lp_operator     => ao_FEM % opo_Op (ai_operator)

        li_sp       = ao_FEM % opi_InfoOperator (ai_operator, INFOOPERATOR_SPACE_1)
        li_spprime  = ao_FEM % opi_InfoOperator (ai_operator, INFOOPERATOR_SPACE_2)

        lp_con          => ao_FEM % opo_spaces (li_sp) % oo_con
        lp_conprime     => ao_FEM % opo_spaces (li_spprime) % oo_con

        li_conelt = lp_con % opi_real_elts (ai_id-1,ai_elt)
        li_conprimeelt = lp_conprime % opi_real_elts (ai_id-1,ai_elt)

!print *, 'elt =', li_conelt
        ! if we do not transpose the matrix
        IF ( ao_FEM % opi_InfoOperator (ai_operator, INFOOPERATOR_TRANSPOSE) == 0 ) THEN
            do li_b = 1, lp_con % opi_nen (ai_id)

                li_A = lp_con % opi_LM ( ai_id, li_b, li_conelt )

                if (li_A == 0) then
                    cycle
                end if

                do li_bprime = 1, lp_conprime % opi_nen (ai_id)

                    li_Aprime = lp_conprime % opi_LM ( ai_id, li_bprime, li_conprimeelt )

                    if (li_Aprime == 0) then
                        cycle
                    end if

                    li_nmatrices = lp_operator % opi_matrices_toassembly (0)
                    DO li_i = 1, li_nmatrices
                    li_matrix   = lp_operator % opi_matrices_toassembly (li_i)
                    lr_scale    = lp_operator % opr_scale (li_i)
                    lr_val      = lr_scale * self % opr_Matrix_elt(ai_operator, li_b, li_bprime)

                    ! WE INSERT THE CORRESPONDANT VALUE OF STIFFNESS, MASS MATRICES INTO THE CCR MATRIX
                    CALL SPM_ASSEMBLYSETVALUE(li_matrix, li_A, li_Aprime, lr_val, ierr)
                    END DO

                end do

            end do
        ELSE
            do li_b = 1, lp_con % opi_nen (ai_id)

                li_A = lp_con % opi_LM ( ai_id, li_b, li_conelt )

                if (li_A == 0) then
                    cycle
                end if

                do li_bprime = 1, lp_conprime % opi_nen (ai_id)

                    li_Aprime = lp_conprime % opi_LM ( ai_id, li_bprime, li_conprimeelt )

                    if (li_Aprime == 0) then
                        cycle
                    end if

                    li_nmatrices = lp_operator % opi_matrices_toassembly (0)
                    DO li_i = 1, li_nmatrices
                    li_matrix   = lp_operator % opi_matrices_toassembly (li_i)
                    lr_scale    = lp_operator % opr_scale (li_i)
                    lr_val      = lr_scale * self % opr_Matrix_elt(ai_operator, li_b, li_bprime)

                    ! WE INSERT THE CORRESPONDANT VALUE OF STIFFNESS, MASS MATRICES INTO THE CCR MATRIX
                    CALL SPM_ASSEMBLYSETVALUE(li_matrix, li_Aprime, li_A, lr_val, ierr)
                    END DO

                end do

            end do
        END IF

    end subroutine Assembly_Matrices_from_elt_basic

end module operators_module

