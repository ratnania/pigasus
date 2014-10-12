!**************************************************
!
!                   OPERATOR MODULE
!
!**************************************************
module operator_module
    use fem_def
    implicit none

    integer, parameter, private :: mi_dtllevel_base = 1

contains

    !---------------------------------------------------------------
    subroutine SET_INFOOPERATOR(self, ai_id, ai_param, ai_val)
        implicit none
        TYPE(FEM) :: self
        INTEGER :: ai_id
        INTEGER :: ai_param
        INTEGER :: ai_val

        self % opi_InfoOperator(ai_id, ai_param) = ai_val
    end subroutine SET_INFOOPERATOR    
    !---------------------------------------------------------------
    subroutine save_operator_toassembly(self, ai_operator)
        IMPLICIT NONE
        TYPE(FEM) :: self 
        INTEGER, intent(in)  :: ai_operator        
        ! LOCAL
        TYPE(OPER), POINTER :: lp_op

        CALL printlog("save_operator_toassembly: Start", ai_dtllevel = mi_dtllevel_base  + 2)

        lp_op => self % opo_Op (ai_operator)

        self % opo_Op (ai_operator) % opi_matrices_toassembly_tmp = self % opo_Op (ai_operator) % opi_matrices_toassembly

        CALL printlog("save_operator_toassembly : End", ai_dtllevel = mi_dtllevel_base  + 2)

    end subroutine save_operator_toassembly
    !---------------------------------------------------------------
    subroutine load_operator_toassembly(self, ai_operator)
        IMPLICIT NONE
        TYPE(FEM) :: self 
        INTEGER, intent(in)  :: ai_operator        
        ! LOCAL

        CALL printlog("load_operator_toassembly: Start", ai_dtllevel = mi_dtllevel_base  + 2)

        self % opo_Op (ai_operator) % opi_matrices_toassembly  = self % opo_Op (ai_operator) % opi_matrices_toassembly_tmp

        CALL printlog("load_operator_toassembly: End", ai_dtllevel = mi_dtllevel_base  + 2)

    end subroutine load_operator_toassembly
    !--------------------------------------------------------------
    !> \todo not finished yet
    subroutine allocate_operator(self, ai_id)
        implicit none
        TYPE(FEM) :: self
        !> OPERATOR REFERENCE IN THE DICTIONNARY
        integer :: ai_id
        ! LOCAL
        integer :: li_sp
        integer :: li_spprime
        integer :: li_size
        integer :: li_sizeprime
        integer :: li_grids
        integer :: li_gridsprime
        integer :: li_npatchs
        integer :: li_patch
        integer :: li_maxaddto
        INTEGER :: ierr

#ifdef _TRACE
        CALL printlog("allocate_operator : Start", ai_dtllevel = mi_dtllevel_base + 1)
#endif

        li_maxaddto     = self % oi_maxaddto

        li_sp           = self % opi_InfoOperator (ai_id, INFOOPERATOR_SPACE_1)
        li_spprime      = self % opi_InfoOperator (ai_id, INFOOPERATOR_SPACE_2)

        ALLOCATE(self % opo_Op (ai_id) % opr_scale(li_maxaddto))
        ALLOCATE(self % opo_Op (ai_id) % opi_matrices_toassembly(0:li_maxaddto))
        ALLOCATE(self % opo_Op (ai_id) % opi_matrices_toassembly_tmp(0:li_maxaddto))

        self % opo_Op (ai_id) % opr_scale = 1.0

#ifdef _TRACE
        CALL printlog("allocate_operator : End", ai_dtllevel = mi_dtllevel_base + 1)
#endif

    end subroutine allocate_operator
    !---------------------------------------------------------------
    subroutine allocate_operators(self)
        implicit none
        TYPE(FEM) :: self
        ! LOCAL
        INTEGER :: li_id
        INTEGER :: ierr 

#ifdef _TRACE
        CALL printlog("allocate_operators : Start", ai_dtllevel = mi_dtllevel_base + 1)
#endif

        DO li_id = 0, self % oi_nOperators - 1
        CALL allocate_operator(self, li_id)
        END DO

#ifdef _TRACE
        CALL printlog("allocate_operators : End", ai_dtllevel = mi_dtllevel_base + 1)
#endif

    end subroutine allocate_operators
    !---------------------------------------------------------------
    subroutine deallocate_operators(self)
        implicit none
        TYPE(FEM) :: self
        ! LOCAL
        INTEGER :: li_id
        INTEGER :: ierr

#ifdef _TRACE
        CALL printlog("deallocate_operators : Start", ai_dtllevel = mi_dtllevel_base + 1)
#endif

        DO li_id = 0, self % oi_nOperators - 1
        DEALLOCATE(self % opo_Op (li_id) % opr_scale)
        DEALLOCATE(self % opo_Op (li_id) % opi_matrices_toassembly)
        DEALLOCATE(self % opo_Op (li_id) % opi_matrices_toassembly_tmp)        
        END DO

#ifdef _TRACE
        CALL printlog("deallocate_operators : End", ai_dtllevel = mi_dtllevel_base + 1)
#endif

    end subroutine deallocate_operators


end module operator_module
!**************************************************
