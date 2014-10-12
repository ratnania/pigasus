!**************************************************
!
!                   SOLVER MODULE
!
!**************************************************
module solver_module
    use fem_def
    implicit none

    integer, parameter, private :: mi_dtllevel_base = 1
    integer, private :: mi_matrix
    integer, private :: mi_N
    INTEGER, PARAMETER :: MAXITER = 10000

contains

    !---------------------------------------------------------------
    subroutine SET_INFOSOLVER(self, ai_id, ai_param, ai_val)
        implicit none
        TYPE(FEM) :: self
        INTEGER :: ai_id
        INTEGER :: ai_param
        INTEGER :: ai_val

        self % opi_InfoSolver(ai_id, ai_param) = ai_val

    end subroutine SET_INFOSOLVER
    !---------------------------------------------------------------
    subroutine SET_SOLVER_MAXITER_fem(self, ai_id, ai_val)
        implicit none
        TYPE(FEM) :: self
        INTEGER :: ai_id
        INTEGER :: ai_val

        self % opo_solver (ai_id) % oi_maxiter = ai_val

    end subroutine SET_SOLVER_MAXITER_fem
    !---------------------------------------------------------------
    subroutine SET_SOLVER_RTOL_fem(self, ai_id, ar_val)
        implicit none
        TYPE(FEM) :: self
        INTEGER :: ai_id
        REAL(8) :: ar_val

        self % opo_solver (ai_id) % or_rtol = ar_val

    end subroutine SET_SOLVER_RTOL_fem
    !---------------------------------------------------------------
    subroutine SET_SOLVER_ATOL_fem(self, ai_id, ar_val)
        implicit none
        TYPE(FEM) :: self
        INTEGER :: ai_id
        REAL(8) :: ar_val

        self % opo_solver (ai_id) % or_atol = ar_val

    end subroutine SET_SOLVER_ATOL_fem
    !---------------------------------------------------------------
    subroutine SET_SOLVER_EPS_fem(self, ai_id, ar_val)
        implicit none
        TYPE(FEM) :: self
        INTEGER :: ai_id
        REAL(8) :: ar_val

        self % opo_solver (ai_id) % or_eps = ar_val

    end subroutine SET_SOLVER_EPS_fem
    !---------------------------------------------------------------
    subroutine GET_SOLVER_NITER_fem(self, ai_id, ai_val)
        implicit none
        TYPE(FEM) :: self
        INTEGER :: ai_id
        INTEGER :: ai_val

        ai_val = self % opo_solver (ai_id) % oi_niter

    end subroutine GET_SOLVER_NITER_fem
    !---------------------------------------------------------------
    subroutine GET_SOLVER_ERR_fem(self, ai_id, apr_val, ai_size)
        implicit none
        TYPE(FEM) :: self
        INTEGER :: ai_id
        INTEGER :: ai_size
        REAL(8), DIMENSION(ai_size) :: apr_val

        apr_val(1:ai_size) = self % opo_solver (ai_id) % opr_err(1:ai_size) 

    end subroutine GET_SOLVER_ERR_fem
    !---------------------------------------------------------------
    subroutine GET_SOLVER_RESNORM_fem(self, ai_id, ar_val)
        implicit none
        TYPE(FEM) :: self
        INTEGER :: ai_id
        REAL(8) :: ar_val

        ar_val = self % opo_solver (ai_id) % or_resnorm

    end subroutine GET_SOLVER_RESNORM_fem
    !---------------------------------------------------------------
    subroutine allocate_solvers(self)
        implicit none
        TYPE(FEM) :: self
        ! LOCAL
        INTEGER :: li_id
        INTEGER :: ierr 
        INTEGER :: li_solver
        INTEGER :: li_matrix
        INTEGER :: li_type
        INTEGER :: li_residual
        INTEGER :: li_N
        TYPE(SOLVER), POINTER :: lp_slv
        

#ifdef _TRACE
        CALL printlog("allocate_solvers : Start", ai_dtllevel = mi_dtllevel_base + 1)
#endif

        ALLOCATE(self % opo_solver(0:self%oi_nsolvers-1))

        DO li_id = 0, self % oi_nsolvers - 1
                lp_slv => self % opo_solver(li_id)
                lp_slv % ol_allocated = .FALSE.

                li_matrix = self % opi_InfoSolver(li_id, INFOSOLVER_MATRIX)
                li_type = self % opi_InfoMatrix(li_matrix, INFOMATRIX_TYPE)
                li_residual = self % opi_InfoSolver(li_id, INFOSOLVER_RESIDUAL)
                li_solver = self % opi_InfoSolver(li_id, INFOSOLVER_SOLVER)

                IF (li_type == GENERIC_MATRIX) THEN
                        CALL SPM_GetnR (li_matrix, li_N , ierr)
                        ALLOCATE(lp_slv % opr_U(li_N)) 
                        ALLOCATE(lp_slv % opr_B(li_N))
!                        print *, 'MATRIX = ', li_matrix
!                        print *, 'N = ', li_N 
                        lp_slv % ol_allocated = .TRUE.
                END IF

                ALLOCATE(lp_slv % opr_err(MAXITER))
                lp_slv % ol_allocatedErr = .TRUE.                

                lp_slv % or_eps = 1.0E-7
        END DO

#ifdef _TRACE
        CALL printlog("allocate_solvers : End", ai_dtllevel = mi_dtllevel_base + 1)
#endif

    end subroutine allocate_solvers
    !---------------------------------------------------------------
    subroutine deallocate_solvers(self)
        implicit none
        TYPE(FEM) :: self
        ! LOCAL
        INTEGER :: li_id
        INTEGER :: ierr 
        INTEGER :: li_solver
        INTEGER :: li_matrix
        INTEGER :: li_type        
        INTEGER :: li_residual
        LOGICAL :: ll_allocatedErr
        TYPE(SOLVER), POINTER :: lp_slv

#ifdef _TRACE
        CALL printlog("deallocate_solvers : Start", ai_dtllevel = mi_dtllevel_base + 1)
#endif

        DO li_id = 0, self % oi_nsolvers - 1
                lp_slv => self % opo_solver(li_id)

                li_matrix = self % opi_InfoSolver(li_id, INFOSOLVER_MATRIX)
                li_type = self % opi_InfoMatrix(li_matrix, INFOMATRIX_TYPE)
                li_residual = self % opi_InfoSolver(li_id, INFOSOLVER_RESIDUAL)
                li_solver = self % opi_InfoSolver(li_id, INFOSOLVER_SOLVER)

                ll_allocatedErr = lp_slv % ol_allocatedErr

                IF (lp_slv % ol_allocated) THEN
                        DEALLOCATE(lp_slv % opr_U) 
                        DEALLOCATE(lp_slv % opr_B)
                        lp_slv % ol_allocated = .FALSE.
                END IF

                IF (ll_allocatedErr) THEN
                        DEALLOCATE(lp_slv % opr_err)
                END IF

        END DO

        DEALLOCATE(self % opo_solver)

#ifdef _TRACE
        CALL printlog("deallocate_solvers : End", ai_dtllevel = mi_dtllevel_base + 1)
#endif

    end subroutine deallocate_solvers    
    !----------------------------------------------------------------------------------------------
    subroutine solver_getsolution_fem(self, ai_id, apr_val, ai_size)
        implicit none
        TYPE(FEM) :: self
        integer, intent(in) :: ai_id
        integer, intent ( in) :: ai_size
        real(wp), dimension ( ai_size), intent ( out) :: apr_val
        ! LOCAL
        INTEGER(KIND=SPM_INTS_KIND) :: root
        INTEGER(KIND=SPM_INTS_KIND) :: ierr

!        root = -1
!        CALL SPM_GETGLOBALSOLUTION(ai_matrix_id, apr_val, root, ierr)

        apr_val = self % opo_solver(ai_id) % opr_U

    end subroutine solver_getsolution_fem
    !----------------------------------------------------------------------------------------------
    subroutine solver_getfieldsolution_fem(self, ai_id, ai_field_id)
        implicit none
        TYPE(FEM) :: self
        integer, intent(in) :: ai_id
        integer, intent(in) :: ai_field_id
        ! LOCAL
        INTEGER(KIND=SPM_INTS_KIND) :: root
        INTEGER(KIND=SPM_INTS_KIND) :: ierr

!        root = -1
!        CALL SPM_GETGLOBALSOLUTION(ai_matrix_id, self % opo_F(ai_field_id) % opr_c(:), root, ierr)

        self % opo_F(ai_field_id) % opr_c = self % opo_solver(ai_id) % opr_U

    end subroutine solver_getfieldsolution_fem
    !----------------------------------------------------------------------------------------------
    subroutine solver_setvaluesrhs_fem(self, ai_id, apr_val, ai_size)
        implicit none
        TYPE(FEM) :: self
        integer, intent(in) :: ai_id
        integer, intent (in) :: ai_size
        real(wp), dimension (ai_size), intent ( in) :: apr_val
        ! LOCAL
        INTEGER(KIND=SPM_INTS_KIND) :: root
        INTEGER(KIND=SPM_INTS_KIND) :: ierr
        INTEGER :: li_matrix        
        INTEGER :: li_N
        TYPE(SOLVER), POINTER :: lp_slv

!        root = -1
!        CALL SPM_SETGLOBALRHS(ai_matrix_id, apr_val, root, SPM_ASSEMBLY_OVW, ierr)

        lp_slv => self % opo_solver(ai_id)
        IF (.not. lp_slv % ol_allocated ) THEN
                li_matrix = self % opi_InfoSolver(ai_id, INFOSOLVER_MATRIX)
                CALL SPM_GetnR (li_matrix, li_N , ierr)
                ALLOCATE(lp_slv % opr_U(li_N)) 
                ALLOCATE(lp_slv % opr_B(li_N))
                lp_slv % ol_allocated = .TRUE.
        END IF

!        print *, 'matrix: ', ai_id
!        print *, SIZE(apr_val,1) 
!        print *, SIZE(lp_slv % opr_B,1) 
        lp_slv % opr_B = apr_val        

    end subroutine solver_setvaluesrhs_fem
    !----------------------------------------------------------------------------------------------
    subroutine solver_setfieldrhs_fem(self, ai_id, ai_field_id)
        implicit none
        TYPE(FEM) :: self
        integer, intent(in) :: ai_id
        integer, intent(in) :: ai_field_id
        ! LOCAL
        INTEGER(KIND=SPM_INTS_KIND) :: root
        INTEGER(KIND=SPM_INTS_KIND) :: ierr

!        root = -1
!        CALL SPM_SETGLOBALRHS(ai_matrix_id, self % opo_F(ai_field_id) % opr_c(:) &
!        , root, SPM_ASSEMBLY_OVW, ierr)

        self % opo_solver(ai_id) % opr_B = self % opo_F(ai_field_id) % opr_c

    end subroutine solver_setfieldrhs_fem
    !----------------------------------------------------------------------------------------------
    subroutine solver_initialize_fem(self, ai_id)
        implicit none
        TYPE(FEM) :: self
        integer, intent(in) :: ai_id
        ! LOCAL
        INTEGER(KIND=SPM_INTS_KIND) :: ierr
        INTEGER :: li_matrix
        INTEGER :: li_residual
        INTEGER :: li_maxiter
        INTEGER :: li_solver
        REAL(8)    :: lr_eps
        TYPE(SOLVER), POINTER :: lp_slv

        lp_slv => self % opo_solver(ai_id)

        li_matrix       = self % opi_InfoSolver(ai_id, INFOSOLVER_MATRIX)
        li_residual     = self % opi_InfoSolver(ai_id, INFOSOLVER_RESIDUAL)
        li_solver       = self % opi_InfoSolver(ai_id, INFOSOLVER_SOLVER)

        li_maxiter      = lp_slv % oi_maxiter 
        lr_eps          = lp_slv % or_eps

        CALL SPM_SetOptionINT(li_matrix, SPM_IPARAM_SOLVER, li_solver, ierr)
        CALL SPM_SetOptionINT(li_matrix, SPM_IPARAM_SOLVER_MAXITER, li_maxiter, ierr)
        CALL SPM_SetOptionREAL(li_matrix, SPM_RPARAM_SOLVER_RTOL, lr_eps, ierr)

        IF (li_maxiter > MAXITER) THEN
                DEALLOCATE(lp_slv % opr_err)
                ALLOCATE(lp_slv % opr_err(li_maxiter))
        END IF

        CALL SPM_INITIALIZE_SOLVE(li_matrix, ierr)

    end subroutine solver_initialize_fem    
    !----------------------------------------------------------------------------------------------
    subroutine solver_run_fem(self, ai_id)
        implicit none
        TYPE(FEM) :: self
        integer, intent(in) :: ai_id
        ! LOCAL
        INTEGER(KIND=SPM_INTS_KIND) :: root
        INTEGER(KIND=SPM_INTS_KIND) :: ierr
        INTEGER :: li_matrix
        INTEGER :: li_residual
        INTEGER :: li_N
        TYPE(SOLVER), POINTER :: lp_slv

        lp_slv => self % opo_solver(ai_id)

        root = -1

        li_matrix       = self % opi_InfoSolver(ai_id, INFOSOLVER_MATRIX)
        li_residual     = self % opi_InfoSolver(ai_id, INFOSOLVER_RESIDUAL)

        CALL SPM_SOLVE(li_matrix, lp_slv % opr_U, lp_slv % opr_B, ierr)

        CALL SPM_GetOptionINT(li_matrix, SPM_IPARAM_SOLVER_NITER, lp_slv % oi_niter, ierr)
        CALL SPM_GetOptionREAL(li_matrix, SPM_RPARAM_SOLVER_RESNORM, lp_slv % or_resnorm, ierr)        
        ! TODO a remplir lp_slv % opr_err

    end subroutine solver_run_fem
    !----------------------------------------------------------------------------------------------
    subroutine solver_finalize_fem(self, ai_id)
        implicit none
        TYPE(FEM) :: self
        integer, intent(in) :: ai_id
        ! LOCAL
        INTEGER(KIND=SPM_INTS_KIND) :: ierr
        INTEGER :: li_matrix
        INTEGER :: li_residual

        li_matrix       = self % opi_InfoSolver(ai_id, INFOSOLVER_MATRIX)
        li_residual     = self % opi_InfoSolver(ai_id, INFOSOLVER_RESIDUAL)

        CALL SPM_FREE_SOLVE(li_matrix, ierr)

    end subroutine solver_finalize_fem

end module solver_module
!**************************************************

