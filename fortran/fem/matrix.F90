!**************************************************
!
!                   MATRIX MODULE
!
!**************************************************
module matrix_module
    use fem_def
    implicit none

    integer, parameter, private :: mi_dtllevel_base = 1

contains

    !---------------------------------------------------------------
    subroutine SET_INFOMATRIX(self, ai_id, ai_param, ai_val)
        implicit none
        TYPE(FEM) :: self
        INTEGER :: ai_id
        INTEGER :: ai_param
        INTEGER :: ai_val

        self % opi_InfoMatrix(ai_id, ai_param) = ai_val
    end subroutine SET_INFOMATRIX
    !----------------------------------------------------------------------------------------------
    !> \todo not finished yet
    subroutine allocate_matrix(self, ai_id)
        implicit none
        TYPE(FEM) :: self
        !> MATRIX REFERENCE IN THE DICTIONNARY
        integer :: ai_id
        ! LOCAL
        integer :: li_type
        integer :: li_graph
        INTEGER :: li_ijvassembly
        INTEGER, PARAMETER :: li_root = -1
        INTEGER :: ierr

#ifdef _TRACE
        CALL printlog("allocate_matrix : Start", ai_dtllevel = mi_dtllevel_base + 1)
#endif

        li_type  = self % opi_InfoMatrix (ai_id, INFOMATRIX_TYPE)
        IF (li_type==GENERIC_MATRIX) THEN
        li_graph = self % opi_InfoMatrix (ai_id, INFOMATRIX_GRAPH)
        li_ijvassembly = self % opi_InfoMatrix (ai_id, INFOMATRIX_IJVASSEMBLY)

        IF (li_ijvassembly==1) THEN
                CALL SPM_GraphGlobalIJV(ai_id &
                , self % opo_graph(li_graph) % oi_nR &
                , self % opo_graph(li_graph) % oi_nC &
                , li_root, ierr)
        ELSE
                CALL SPM_GraphGlobalCSR(ai_id &
                , self % opo_graph(li_graph) % oi_nC &
                , self % opo_graph(li_graph) % opi_ja &
                , self % opo_graph(li_graph) % opi_ia &
                , li_root, ierr)
        END IF
        END IF

#ifdef _TRACE
        CALL printlog("allocate_matrix : End", ai_dtllevel = mi_dtllevel_base + 1)
#endif

    end subroutine allocate_matrix
    !---------------------------------------------------------------
    subroutine allocate_matrices(self)
        implicit none
        TYPE(FEM) :: self
        ! LOCAL
        INTEGER :: li_id
        INTEGER :: ierr 

#ifdef _TRACE
        CALL printlog("allocate_matrices : Start", ai_dtllevel = mi_dtllevel_base + 1)
#endif

        ! Starting SPM
        CALL SPM_INITIALIZE(self % oi_nMatrices, ierr)
        IF (ierr /= SPM_SUCCESS) CALL abort()

        DO li_id = 0, self % oi_nMatrices - 1
        CALL allocate_matrix(self, li_id)
        END DO


#ifdef _TRACE
        CALL printlog("allocate_matrices : End", ai_dtllevel = mi_dtllevel_base + 1)
#endif

    end subroutine allocate_matrices
    !---------------------------------------------------------------
    subroutine deallocate_matrices(self)
        implicit none
        TYPE(FEM) :: self
        ! LOCAL
        INTEGER :: li_id
        INTEGER :: ierr

#ifdef _TRACE
        CALL printlog("deallocate_matrices : Start", ai_dtllevel = mi_dtllevel_base + 1)
#endif

        DO li_id = 0, self % oi_nMatrices - 1
        CALL SPM_CLEAN(li_id, ierr)
        END DO

        CALL SPM_FINALIZE(ierr)         

#ifdef _TRACE
        CALL printlog("deallocate_matrices : End", ai_dtllevel = mi_dtllevel_base + 1)
#endif

    end subroutine deallocate_matrices
    !---------------------------------------------------------------
    subroutine multop_matrix_array (self, ai_matrix_id, ai_field_id, apr_val, ai_size)
        implicit none
        TYPE(FEM) :: self
        INTEGER :: ai_matrix_id
        INTEGER :: ai_field_id
        INTEGER :: ai_size
        real(8), dimension (ai_size), intent (out) :: apr_val
        ! LOCAL
        INTEGER :: ierr

        IF (self % opi_InfoMatrix (ai_matrix_id, INFOMATRIX_IJVASSEMBLY) /= 1) THEN 
                CALL SPM_GEMV(ai_matrix_id, 0, self % opo_F(ai_field_id) % opr_c, apr_val, ierr)
        ELSE
                PRINT *, "multop_matrix_array : must stop. You are using IJV"
!                assembly and matrix operations are not handeled here. Please
!                Check that pyMurge is used'
                STOP
        END IF

    end subroutine multop_matrix_array
    !---------------------------------------------------------------
    subroutine multop_matrix_array2 (self, ai_matrix_id, apr_val_in, ai_size_in, apr_val_out, ai_size_out)
        implicit none
        TYPE(FEM) :: self
        INTEGER :: ai_matrix_id
        INTEGER :: ai_field_id
        INTEGER :: ai_size_in
        real(8), dimension (ai_size_in), intent (in) :: apr_val_in
        INTEGER :: ai_size_out
        real(8), dimension (ai_size_out), intent (out) :: apr_val_out
        ! LOCAL
        INTEGER :: ierr

        IF (self % opi_InfoMatrix (ai_matrix_id, INFOMATRIX_IJVASSEMBLY) /= 1) THEN 
                CALL SPM_GEMV(ai_matrix_id, 0, apr_val_in, apr_val_out, ierr)
        ELSE
                PRINT *, "multop_matrix_array : must stop. You are using IJV"
!                'assembly and matrix operations are not handeled here. Please' &
!                'Check that pyMurge is used'
                STOP
        END IF
    end subroutine multop_matrix_array2
    !---------------------------------------------------------------
    subroutine multop_matrix_field(self, ai_id, ai_j, ai_k)
        !> this routine implements the multiplication opperation :
        !> F(k) := M(id) *  F(j)
        !> where M(id) is a matrix and
        !> F(j), F(k) are 2 fields
        implicit none
        TYPE(FEM) :: self
        INTEGER :: ai_id
        INTEGER :: ai_j
        INTEGER :: ai_k

        PRINT *, 'multop_matrix_matrix: Not yet implemented'

!        CALL Mult_MV(self % opo_CCR_Matrix(ai_id) &
!        , self % opo_F(ai_j) % opr_c(:) &
!        , self % opo_F(ai_k) % opr_c(:))

    end subroutine multop_matrix_field
    !---------------------------------------------------------------
    subroutine addop_matrix_matrix(self, ai_i, ai_j, ai_k)
        !> this routine implements the addition opperation :
        !> M(k) := M(i) + M(j)
        !> where M(i), M(j), M(k) are a matrices
        implicit none
        TYPE(FEM) :: self
        INTEGER :: ai_i
        INTEGER :: ai_j
        INTEGER :: ai_k

        PRINT *, 'multop_matrix_matrix: Not yet implemented'
!        CALL add_csr_Matrix(self % opo_CCR_Matrix(ai_i) &
!        , self % opo_CCR_Matrix(ai_j) &
!        , self % opo_CCR_Matrix(ai_k))

    end subroutine addop_matrix_matrix
    !---------------------------------------------------------------
    subroutine multop_matrix_matrix(self, ai_i, ai_j, ai_k)
        !> this routine implements the addition opperation :
        !> M(k) := M(i) + M(j)
        !> where M(i), M(j), M(k) are a matrices
        implicit none
        TYPE(FEM) :: self
        INTEGER :: ai_i
        INTEGER :: ai_j
        INTEGER :: ai_k

        PRINT *, 'multop_matrix_matrix: Not yet implemented'
!        CALL mult_csr_Matrix(self % opo_CCR_Matrix(ai_i) &
!        , self % opo_CCR_Matrix(ai_j) &
!        , self % opo_CCR_Matrix(ai_k))

    end subroutine multop_matrix_matrix
    !---------------------------------------------------------------
    subroutine multop_matrix_scal (self, ai_i, ar_val, ai_k)
        implicit none
        TYPE(FEM) :: self
        INTEGER :: ai_i
        REAL(8) :: ar_val
        INTEGER :: ai_k
        PRINT *, 'multop_matrix_scal: Not yet implemented'

!        CALL Mult_Mscal(self % opo_CCR_Matrix(ai_i) &
!        , ar_val &
!        , self % opo_CCR_Matrix(ai_k))

    end subroutine multop_matrix_scal
    !---------------------------------------------------------------
    subroutine addop_matrix_scal (self, ai_i, ar_val, ai_k)
        implicit none
        TYPE(FEM) :: self
        INTEGER :: ai_i
        REAL(8) :: ar_val
        INTEGER :: ai_k
        PRINT *, 'addop_matrix_scal: Not yet implemented'

!        CALL Add_Mscal(self % opo_CCR_Matrix(ai_i) &
!        , ar_val &
!        , self % opo_CCR_Matrix(ai_k))

    end subroutine addop_matrix_scal

end module matrix_module
!**************************************************
