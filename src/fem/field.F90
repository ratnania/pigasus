!**************************************************
!
!                   FIELD MODULE
!
!**************************************************
module field_module
    use fem_def
    use spaces_module
    IMPLICIT NONE

    INTEGER, parameter, private :: mi_dtllevel_base = 1

contains

    !---------------------------------------------------------------
    SUBROUTINE SET_INFOFIELD(self, ai_id, ai_param, ai_val)
        IMPLICIT NONE
        TYPE(FEM) :: self
        INTEGER :: ai_id
        INTEGER :: ai_param
        INTEGER :: ai_val

        self % opi_InfoField(ai_id, ai_param) = ai_val

    END SUBROUTINE SET_INFOFIELD
    !---------------------------------------------------------------
    SUBROUTINE allocate_fields(self)
        IMPLICIT NONE
        TYPE(FEM) :: self
        ! LOCAL
        INTEGER :: li_Field
        INTEGER :: li_size
        INTEGER :: li_ndof
        INTEGER :: li_space
        INTEGER :: li_grids
        INTEGER :: li_dim
        INTEGER :: li_npatchs
        INTEGER :: li_maxnderiv
        INTEGER :: li_maxnel
        INTEGER :: li_maxnpts
#ifdef _TRACE
        CALL printlog("allocate_fields : Start", ai_dtllevel = mi_dtllevel_base + 1)
#endif

        ALLOCATE(self % opo_F(0:self % oi_nFields - 1))

        DO li_field = 0, self % oi_nFields - 1
            li_size     = self % opi_InfoField(li_Field, INFOFIELD_SIZE)
            li_ndof     = self % opi_InfoField(li_Field, INFOFIELD_NDOF)

            li_space    = self % opi_InfoField(li_Field, INFOFIELD_SPACE)
            li_grids    = self % opi_InfoSpace(li_space, INFOSPACE_GRIDS)
            li_dim      = self % opi_InfoGrids(li_grids, INFOGRIDS_DIM)
            li_npatchs  = self % opi_InfoGrids(li_grids, INFOGRIDS_NPATCHS)

            CALL get_maxnel_space(self, li_space, li_maxnel)

            li_maxnpts = MAXVAL (self % opi_InfoPatch(li_grids, 0:li_npatchs-1, INFOPATCH_MAXNPTS))
!            li_maxnpts = MIN(li_maxnpts, 1000)
!            print *, '----- grids ',li_grids,' ------'
!            print *, 'maxnpts : ', self % opi_InfoPatch(li_grids, 0:li_npatchs-1, INFOPATCH_MAXNPTS)

!            li_maxnderiv = li_dim ! this works for CURL, GRAD: ONLY 1st ORDER OPERATORS
            li_maxnderiv = li_dim**2 ! this works for CURL, GRAD and SECOND_ORDER_FIELD

            self % opo_F(li_field) % oi_size = li_size

            ALLOCATE(self % opo_F(li_field) % opr_c(li_size))
            ALLOCATE(self % opo_F(li_field) % opi_counter(li_size))

            CALL reset_field(self, li_field)

!            print *,"ndof      =", li_ndof
!            print *,"maxnderiv =", li_maxnderiv
!            print *,"maxnel    =", li_maxnel
!            print *,"maxnpts   =", li_maxnpts

            ALLOCATE(self % opo_F(li_field) % opr_val(li_ndof, li_maxnderiv, li_maxnel, li_maxnpts))
            self % opo_F(li_field) % opr_val (:,:,:,:) = 0.0_wp

        END DO
#ifdef _TRACE
        CALL printlog("allocate_fields : End", ai_dtllevel = mi_dtllevel_base + 1)
#endif

    END SUBROUTINE allocate_fields
    !---------------------------------------------------------------
    SUBROUTINE deallocate_fields(self)
        IMPLICIT NONE
        TYPE(FEM) :: self
        ! LOCAL
        INTEGER :: li_Field
#ifdef _TRACE
        CALL printlog("deallocate_fields : Start", ai_dtllevel = mi_dtllevel_base + 1)
#endif

        DO li_field = 0, self % oi_nFields - 1
            DEALLOCATE(self % opo_F(li_field) % opr_c)
            DEALLOCATE(self % opo_F(li_field) % opr_val)
        END DO

        DEALLOCATE(self % opo_F)
#ifdef _TRACE
        CALL printlog("deallocate_fields : End", ai_dtllevel = mi_dtllevel_base + 1)
#endif

    END SUBROUTINE deallocate_fields
    !---------------------------------------------------------------
    SUBROUTINE field_setval_coefs_fem(self, ai_id, ar_val)
        !> this routine implements the addition opperation :
        !> F(k) := F(i) + F(j)
        IMPLICIT NONE
        TYPE(FEM) :: self
        INTEGER :: ai_id
        REAL(wp) :: ar_val

        self % opo_F(ai_id) % opr_c(:) = ar_val

    END SUBROUTINE field_setval_coefs_fem
    !----------------------------------------------------------------------------------------------
    SUBROUTINE field_setarray_coefs_fem(self, ai_id, apr_val, ai_size)
        IMPLICIT NONE
        TYPE(FEM) :: self
        INTEGER, INTENT(IN) :: ai_id
        INTEGER, INTENT (IN) :: ai_size
        REAL(wp), DIMENSION (ai_size), INTENT ( in) :: apr_val

        self % opo_F(ai_id) % opr_c(1:self % opo_F(ai_id) % oi_size) = &
        apr_val(1:self % opo_F(ai_id) % oi_size)

    END SUBROUTINE field_setarray_coefs_fem
    !---------------------------------------------------------------
    SUBROUTINE field_setval_values_fem(self, ai_id, ar_val)
        !> this routine implements the addition opperation :
        !> F(k) := F(i) + F(j)
        IMPLICIT NONE
        TYPE(FEM) :: self
        INTEGER :: ai_id
        REAL(wp) :: ar_val

        self % opo_F(ai_id) % opr_val(:,:,:,:) = ar_val

    END SUBROUTINE field_setval_values_fem
    !----------------------------------------------------------------------------------------------
    SUBROUTINE field_setarray_values_fem(self, ai_id, apr_val, ai_ndof, ai_maxnderiv, ai_maxnel, ai_maxnpts)
        IMPLICIT NONE
        TYPE(FEM) :: self
        INTEGER, INTENT(IN) :: ai_id
        INTEGER, INTENT (IN) :: ai_ndof
        INTEGER, INTENT (IN) :: ai_maxnderiv
        INTEGER, INTENT (IN) :: ai_maxnel
        INTEGER, INTENT (IN) :: ai_maxnpts
        REAL(wp), DIMENSION (ai_ndof, ai_maxnderiv, ai_maxnel, ai_maxnpts), INTENT ( in) :: apr_val
        ! LOCAL VARIABLES
        INTEGER :: li_ndof
        INTEGER :: li_maxnderiv
        INTEGER :: li_maxnel
        INTEGER :: li_maxnpts
        INTEGER :: li_space
        INTEGER :: li_dim
        INTEGER :: li_grids

        li_ndof = self % opi_InfoField(ai_id, INFOFIELD_NDOF)
        li_space = self % opi_InfoField(ai_id, INFOFIELD_SPACE)
        li_grids = self % opi_InfoSpace(li_space, INFOSPACE_GRIDS)
        li_dim   = self % opi_InfoGrids(li_grids, INFOGRIDS_DIM)

        CALL get_maxnel_space(self, li_space, li_maxnel)

        li_maxnpts = MAXVAL (self % opi_InfoPatch(li_grids, :, INFOPATCH_MAXNPTS))

        li_maxnderiv = li_dim ! this works for CURL, GRAD: ONLY 1st ORDER OPERATORS

        self % opo_F(ai_id) % opr_val(1:li_ndof, 1:li_maxnderiv, 1:li_maxnel, 1:li_maxnpts) = &
        apr_val(1:li_ndof, 1:li_maxnderiv, 1:li_maxnel, 1:li_maxnpts)

    END SUBROUTINE field_setarray_values_fem
    !----------------------------------------------------------------------------------------------
    SUBROUTINE getfield_fem(self, ai_field, apr_result, ai_size)
        IMPLICIT NONE
        TYPE(FEM) :: self
        INTEGER, INTENT(IN) :: ai_field
        INTEGER, INTENT ( in) :: ai_size
        REAL(wp), DIMENSION ( ai_size), INTENT ( out) :: apr_result

!        print *, 'SIZE Field : ', SIZE(self % opo_F(ai_field) % opr_c , 1)
!        print *, 'ai_size : ', ai_size
        apr_result(1:ai_size) = self % opo_F(ai_field) % opr_c(1:ai_size)

    END SUBROUTINE getfield_fem
!----------------------------------------------------------------------------------------------
    SUBROUTINE get_fieldh ( self, ai_field, ai_patch, apr_fieldh, ai_k, ai_ndof, ai_maxder, ai_nel, ai_maxnpts )
        !> ai_k must be equal to ai_ndof * ai_maxder
        IMPLICIT NONE
        TYPE(FEM) :: self
        INTEGER, INTENT(IN)  :: ai_field
        INTEGER, INTENT(IN)  :: ai_patch
        INTEGER, INTENT(IN)  :: ai_k
        INTEGER, INTENT(IN)  :: ai_ndof
        INTEGER, INTENT(IN)  :: ai_maxder
        INTEGER, INTENT(IN)  :: ai_nel
        INTEGER, INTENT(IN)  :: ai_maxnpts
        REAL*8 , DIMENSION( ai_k, ai_nel, ai_maxnpts), INTENT(out)  :: apr_fieldh
        ! LOCAL VARIABLES
        INTEGER :: li_i
        INTEGER :: li_j
        INTEGER :: li_index
#ifdef _TRACE
        CALL printlog("get_fieldh : Start", ai_dtllevel = 1)
#endif

        DO li_i = 1, ai_ndof
            DO li_j = 1, ai_maxder
                li_index = ( li_i - 1 ) * ai_maxder + li_j
                apr_fieldh (li_index, 1:ai_nel, 1:ai_maxnpts) =   &
                self % opo_F (ai_field) % opr_val(li_i, li_j, 1:ai_nel, 1:ai_maxnpts)
!                print *, self % opo_F (ai_field) % opr_val(li_i, li_j, 1:ai_nel, 1:ai_maxnpts)
            END DO
        END DO
#ifdef _TRACE
        CALL printlog("get_fieldh : End", ai_dtllevel = 1)
#endif

    END SUBROUTINE get_fieldh
!----------------------------------------------------------------------------------------------
    SUBROUTINE field_to_matrix_1D_fem ( self, ai_field, ai_patch, apr_matrix, ai_n_1 )
        IMPLICIT NONE
        TYPE(FEM) :: self
        INTEGER, INTENT(IN)  :: ai_field
        INTEGER, INTENT(IN)  :: ai_patch
        INTEGER, INTENT(IN)  :: ai_n_1
        REAL*8 , DIMENSION(ai_n_1), INTENT(out)  :: apr_matrix
        ! LOCAL VARIABLES
        INTEGER  :: li_elt
        INTEGER  :: li_conelt
        INTEGER  :: li_b
        INTEGER  :: li_i
        INTEGER  :: li_j
        INTEGER  :: li_A
        INTEGER  :: li_P
        INTEGER  :: li_nel
        INTEGER  :: li_id
        INTEGER  :: li_map
        INTEGER  :: li_grids
        INTEGER  :: li_space
        INTEGER  :: li_npts
        INTEGER  :: li_pt
        type(GRID_DATA), POINTER :: lp_grid
        type(SPACE), POINTER :: lp_space
        TYPE(CONNECTIVITY), POINTER :: lp_con
#ifdef _TRACE
        CALL printlog("field_to_matrix_1D_fem : Start", ai_dtllevel = 1)
#endif

        ! LOOP THROUG PATCHS
        li_id = ai_patch + 1

        li_space = self % opi_InfoField ( ai_field, INFOFIELD_SPACE )
        li_grids = self % opi_InfoSpace ( li_space, INFOSPACE_GRIDS )

        lp_space => self % opo_spaces ( li_space )
        lp_con   => lp_space % oo_con

        li_nel  = self % opi_InfoPatch ( li_grids, li_id-1, INFOPATCH_NEL)
        DO li_elt = 1, li_nel

            li_conelt = lp_con % opi_REAL_elts (li_id-1,li_elt)

            DO li_b = 1, lp_con % opi_nen (li_id)

                li_P = lp_con % opi_LM ( li_id, li_b, li_conelt )
            ! ...
            ! computing the i,j indices.
            ! for multi-patchs, we need to compute the right A, using the same technic as for i and j
            ! ...
                li_A = lp_con % opi_IEN ( li_id, li_b, li_conelt )

                li_i = li_A + 1
            ! ...
                IF (li_P == 0) THEN
                    apr_matrix ( li_i ) = 0.0
                ELSE
                    apr_matrix ( li_i ) = self % opo_F (ai_field) % opr_c (li_P)
                END IF

            END DO

        END DO

#ifdef _TRACE
        CALL printlog("field_to_matrix_1D_fem : End", ai_dtllevel = 1)
#endif

    END SUBROUTINE field_to_matrix_1D_fem      
!----------------------------------------------------------------------------------------------
    SUBROUTINE field_from_matrix_1D_fem ( self, ai_field, ai_patch, ar_alpha, apr_matrix, ai_n_1 )
    ! does the following operation F <- alpha F + M
        IMPLICIT NONE
        TYPE(FEM) :: self
        INTEGER, INTENT(IN)  :: ai_field
        INTEGER, INTENT(IN)  :: ai_patch
        REAL*8 , INTENT(IN)  :: ar_alpha
        INTEGER, INTENT(IN)  :: ai_n_1
        REAL*8 , DIMENSION(ai_n_1), INTENT(IN)  :: apr_matrix
        ! LOCAL VARIABLES
        INTEGER  :: li_i
        INTEGER  :: li_P
        INTEGER  :: li_A
        INTEGER  :: li_space
        type(SPACE), POINTER :: lp_space
        TYPE(CONNECTIVITY), POINTER :: lp_con
        REAL*8, DIMENSION(:), POINTER  :: lpr_c
        INTEGER, DIMENSION(:), POINTER  :: lpi_counter

#ifdef _TRACE
        CALL printlog("field_from_matrix_1D_fem: Start", ai_dtllevel = 1)
#endif

        li_space    = self % opi_InfoField ( ai_field, INFOFIELD_SPACE )
        lp_space    => self % opo_spaces ( li_space )
        lp_con      => lp_space % oo_con        
        lpr_c       => self % opo_F (ai_field) % opr_c
        lpi_counter => self % opo_F (ai_field) % opi_counter

        DO li_i = 1, ai_n_1

                li_P = li_i
                li_A = lp_con % opi_info (ai_patch) % opi_ID_loc ( li_P )

                IF (li_A > 0) THEN
                    lpr_c (li_A) = ar_alpha * lpr_c(li_A) + apr_matrix ( li_i )
                    lpi_counter(li_A) = lpi_counter(li_A) + 1
                END IF

        END DO

#ifdef _TRACE
        CALL printlog("field_from_matrix_1D_fem: End", ai_dtllevel = 1)
#endif

    END SUBROUTINE field_from_matrix_1D_fem
!----------------------------------------------------------------------------------------------
    SUBROUTINE field_to_matrix_2D_fem ( self, ai_field, ai_patch, apr_matrix, ai_n_1, ai_n_2 )
        IMPLICIT NONE
        TYPE(FEM) :: self
        INTEGER, INTENT(IN)  :: ai_field
        INTEGER, INTENT(IN)  :: ai_patch
        INTEGER, INTENT(IN)  :: ai_n_1
        INTEGER, INTENT(IN)  :: ai_n_2
        REAL*8 , DIMENSION(ai_n_1, ai_n_2), INTENT(out)  :: apr_matrix
        ! LOCAL VARIABLES
        INTEGER  :: li_elt
        INTEGER  :: li_conelt
        INTEGER  :: li_b
        INTEGER  :: li_i
        INTEGER  :: li_j
        INTEGER  :: li_A
        INTEGER  :: li_P
        INTEGER  :: li_nel
        INTEGER  :: li_id
        INTEGER  :: li_map
        INTEGER  :: li_grids
        INTEGER  :: li_space
        INTEGER  :: li_npts
        INTEGER  :: li_pt
        type(GRID_DATA), POINTER :: lp_grid
        type(SPACE), POINTER :: lp_space
        TYPE(CONNECTIVITY), POINTER :: lp_con
#ifdef _TRACE
        CALL printlog("field_to_matrix_2D_fem : Start", ai_dtllevel = 1)
#endif

        ! LOOP THROUG PATCHS
        li_id = ai_patch + 1

        li_space = self % opi_InfoField ( ai_field, INFOFIELD_SPACE )
        li_grids = self % opi_InfoSpace ( li_space, INFOSPACE_GRIDS )

        lp_space => self % opo_spaces ( li_space )
        lp_con   => lp_space % oo_con

        li_nel  = self % opi_InfoPatch ( li_grids, li_id-1, INFOPATCH_NEL)
        DO li_elt = 1, li_nel

            li_conelt = lp_con % opi_REAL_elts (li_id-1,li_elt)

            DO li_b = 1, lp_con % opi_nen (li_id)

                li_P = lp_con % opi_LM ( li_id, li_b, li_conelt )
            ! ...
            ! computing the i,j indices.
            ! for multi-patchs, we need to compute the right A, using the same technic as for i and j
            ! ...
                li_A = lp_con % opi_IEN ( li_id, li_b, li_conelt )

                li_j = int (li_A / ai_n_1) + 1
                if (li_j==1) then
                    li_i = li_A + 1
                else
                    li_i = li_A - ( li_j - 1 ) * ai_n_1 + 1
                END if
            ! ...
                IF (li_P == 0) THEN
                    apr_matrix ( li_i, li_j ) = 0.0
                ELSE
                    apr_matrix ( li_i, li_j ) = self % opo_F (ai_field) % opr_c (li_P)
                END IF

            END DO

        END DO

#ifdef _TRACE
        CALL printlog("field_to_matrix_2D_fem : End", ai_dtllevel = 1)
#endif

    END SUBROUTINE field_to_matrix_2D_fem
!----------------------------------------------------------------------------------------------
    SUBROUTINE field_from_matrix_2D_fem( self, ai_field, ai_patch, ar_alpha, apr_matrix, ai_n_1, ai_n_2 )
    ! does the following operation F <- alpha F + M
        IMPLICIT NONE
        TYPE(FEM) :: self
        INTEGER, INTENT(IN)  :: ai_field
        INTEGER, INTENT(IN)  :: ai_patch
        REAL*8 , INTENT(IN)  :: ar_alpha
        INTEGER, INTENT(IN)  :: ai_n_1
        INTEGER, INTENT(IN)  :: ai_n_2
        REAL*8 , DIMENSION(ai_n_1, ai_n_2), INTENT(IN)  :: apr_matrix
        ! LOCAL VARIABLES
        INTEGER  :: li_i
        INTEGER  :: li_j
        INTEGER  :: li_A
        INTEGER  :: li_P
        INTEGER  :: li_space
        type(SPACE), POINTER :: lp_space
        TYPE(CONNECTIVITY), POINTER :: lp_con
        INTEGER, DIMENSION(:), POINTER :: lpi_ID_loc
        REAL*8, DIMENSION(:), POINTER  :: lpr_c
        INTEGER, DIMENSION(:), POINTER  :: lpi_counter

#ifdef _TRACE
        CALL printlog("field_from_matrix_2D_fem: Start", ai_dtllevel = 1)
#endif

        li_space    = self % opi_InfoField ( ai_field, INFOFIELD_SPACE )
        lp_space    => self % opo_spaces ( li_space )
        lp_con      => lp_space % oo_con  
        lpi_ID_loc  => lp_con % opi_info (ai_patch) % opi_ID_loc 
        lpr_c       => self % opo_F (ai_field) % opr_c
        lpi_counter => self % opo_F (ai_field) % opi_counter

        DO li_j = 1, ai_n_2
        DO li_i = 1, ai_n_1

                li_P = ( li_j - 1 ) * ai_n_1 + li_i 

                li_A = lpi_ID_loc ( li_P )
                IF (li_A > 0) THEN
                    lpr_c (li_A) = ar_alpha * lpr_c (li_A) + apr_matrix (li_i, li_j)
                    lpi_counter(li_A) = lpi_counter(li_A) + 1
                END IF

        END DO
        END DO

#ifdef _TRACE
        CALL printlog("field_from_matrix_2D_fem: End", ai_dtllevel = 1)
#endif

    end SUBROUTINE field_from_matrix_2D_fem    
!----------------------------------------------------------------------------------------------
    SUBROUTINE field_to_matrix_3D_fem ( self, ai_field, ai_patch, apr_matrix, ai_n_1, ai_n_2 , ai_n_3 )
        IMPLICIT NONE
        TYPE(FEM) :: self
        INTEGER, INTENT(IN)  :: ai_field
        INTEGER, INTENT(IN)  :: ai_patch
        INTEGER, INTENT(IN)  :: ai_n_1
        INTEGER, INTENT(IN)  :: ai_n_2
        INTEGER, INTENT(IN)  :: ai_n_3
        REAL*8 , DIMENSION(ai_n_1, ai_n_2, ai_n_3), INTENT(out)  :: apr_matrix
        ! LOCAL VARIABLES
        INTEGER  :: li_elt
        INTEGER  :: li_conelt
        INTEGER  :: li_b
        INTEGER  :: li_i
        INTEGER  :: li_j
        INTEGER  :: li_k
        INTEGER  :: li_r
        INTEGER  :: li_A
        INTEGER  :: li_P
        INTEGER  :: li_nel
        INTEGER  :: li_id
        INTEGER  :: li_map
        INTEGER  :: li_grids
        INTEGER  :: li_space
        INTEGER  :: li_npts
        INTEGER  :: li_pt
        type(GRID_DATA), POINTER :: lp_grid
        type(SPACE), POINTER :: lp_space
        TYPE(CONNECTIVITY), POINTER :: lp_con
#ifdef _TRACE
        CALL printlog("field_to_matrix_3D_fem : Start", ai_dtllevel = 1)
#endif

        ! LOOP THROUG PATCHS
        li_id = ai_patch + 1

        li_space = self % opi_InfoField ( ai_field, INFOFIELD_SPACE )
        li_grids = self % opi_InfoSpace ( li_space, INFOSPACE_GRIDS )

        lp_space => self % opo_spaces ( li_space )
        lp_con   => lp_space % oo_con

        li_nel  = self % opi_InfoPatch ( li_grids, li_id-1, INFOPATCH_NEL)

        DO li_elt = 1, li_nel

            li_conelt = lp_con % opi_REAL_elts (li_id-1,li_elt)

            DO li_b = 1, lp_con % opi_nen (li_id)

                li_P = lp_con % opi_LM ( li_id, li_b, li_conelt )

            ! ...
            ! computing the i,j,k indices.
            ! for multi-patchs, we need to compute the right A, using the same technic as for i and j,k
            ! ...
                li_A = lp_con % opi_IEN ( li_id, li_b, li_conelt )

                li_k = int (li_A / ai_n_1 * ai_n_2 ) + 1
                li_r = mod (li_A, li_k-1) + 1
                li_j = int (li_r / ai_n_1) + 1
                li_i = mod (li_r, li_j-1) + 1
            ! ...

                IF (li_P == 0) THEN
                    apr_matrix ( li_i, li_j, li_k ) = 0.0
                ELSE
                    apr_matrix ( li_i, li_j, li_k ) = self % opo_F (ai_field) % opr_c (li_P)
                END IF

            END DO

        END DO

#ifdef _TRACE
        CALL printlog("field_to_matrix_3D_fem : End", ai_dtllevel = 1)
#endif

    end SUBROUTINE field_to_matrix_3D_fem    
!----------------------------------------------------------------------------------------------
    SUBROUTINE field_from_matrix_3D_fem( self, ai_field, ai_patch, ar_alpha, apr_matrix, ai_n_1, ai_n_2 , ai_n_3 )
    ! does the following operation F <- alpha F + M
        IMPLICIT NONE
        TYPE(FEM) :: self
        INTEGER, INTENT(IN)  :: ai_field
        INTEGER, INTENT(IN)  :: ai_patch
        REAL*8 , INTENT(IN)  :: ar_alpha
        INTEGER, INTENT(IN)  :: ai_n_1
        INTEGER, INTENT(IN)  :: ai_n_2
        INTEGER, INTENT(IN)  :: ai_n_3
        REAL*8 , DIMENSION(ai_n_1, ai_n_2, ai_n_3), INTENT(IN)  :: apr_matrix
        ! LOCAL VARIABLES
        INTEGER  :: li_i
        INTEGER  :: li_j
        INTEGER  :: li_k
        INTEGER  :: li_A
        INTEGER  :: li_P
        INTEGER  :: li_space
        type(SPACE), POINTER :: lp_space
        TYPE(CONNECTIVITY), POINTER :: lp_con
        INTEGER, DIMENSION(:), POINTER :: lpi_ID_loc
        REAL*8, DIMENSION(:), POINTER  :: lpr_c
        INTEGER, DIMENSION(:), POINTER  :: lpi_counter

#ifdef _TRACE
        CALL printlog("field_from_matrix_3D_fem: Start", ai_dtllevel = 1)
#endif

        li_space    = self % opi_InfoField ( ai_field, INFOFIELD_SPACE )
        lp_space    => self % opo_spaces ( li_space )
        lp_con      => lp_space % oo_con   
        lpi_ID_loc  => lp_con % opi_info (ai_patch) % opi_ID_loc 
        lpr_c       => self % opo_F (ai_field) % opr_c
        lpi_counter => self % opo_F (ai_field) % opi_counter

        DO li_k = 1, ai_n_3
        DO li_j = 1, ai_n_2
        DO li_i = 1, ai_n_1

                li_P = ( li_k - 1 ) * ai_n_1 * ai_n_2 + ( li_j - 1 ) * ai_n_1 + li_i 

                li_A = lpi_ID_loc ( li_P )
                IF (li_A > 0) THEN
                    lpr_c (li_A) = ar_alpha * lpr_c (li_A) + apr_matrix ( li_i, li_j, li_k )
                    lpi_counter(li_A) = lpi_counter(li_A) + 1
                END IF

        END DO
        END DO
        END DO

#ifdef _TRACE
        CALL printlog("field_from_matrix_3D_fem: End", ai_dtllevel = 1)
#endif

    end SUBROUTINE field_from_matrix_3D_fem        
    !----------------------------------------------------------------------------
    SUBROUTINE reset_field(self, ai_field)
        IMPLICIT NONE
        TYPE(FEM) :: self
        INTEGER, INTENT(IN) :: ai_field

        self % opo_F (ai_field) % opr_c         = 0.0_wp
        self % opo_F (ai_field) % opi_counter   = 0

    end SUBROUTINE reset_field
    !----------------------------------------------------------------------------
    SUBROUTINE set_initialize_field(self, ai_field)
        IMPLICIT NONE
        TYPE(FEM) :: self
        INTEGER, INTENT(IN) :: ai_field

        CALL reset_field(self, ai_field)

    end SUBROUTINE set_initialize_field
    !----------------------------------------------------------------------------
    SUBROUTINE set_finalize_field(self, ai_field)
        IMPLICIT NONE
        TYPE(FEM) :: self
        INTEGER, INTENT(IN) :: ai_field
        ! LOCAL
        INTEGER :: li_size
        INTEGER :: li_A
        INTEGER, DIMENSION(:), POINTER :: lpi_counter
        REAL*8 , DIMENSION(:), POINTER :: lpr_c

        li_size     = self % opo_F(ai_field) % oi_size
        lpr_c       => self % opo_F (ai_field) % opr_c
        lpi_counter => self % opo_F (ai_field) % opi_counter

!        print *, lpi_counter

        DO li_A = 1, li_size
        IF (lpi_counter(li_A) > 0) THEN
                lpr_c(li_A) = lpr_c(li_A) / lpi_counter(li_A)
        END IF
        END DO 

    end SUBROUTINE set_finalize_field
    !----------------------------------------------------------------------------
    SUBROUTINE addop_field2(self, ai_i, ai_j, apr_result, ai_size)
        IMPLICIT NONE
        TYPE(FEM) :: self
        INTEGER, INTENT(IN) :: ai_i
        INTEGER, INTENT(IN) :: ai_j
        INTEGER, INTENT ( in) :: ai_size
        REAL(8), DIMENSION ( ai_size), INTENT ( out) :: apr_result

        apr_result(:) = self % opo_F(ai_i) % opr_c(:) &
        +self % opo_F(ai_j) % opr_c(:)

    end SUBROUTINE addop_field2
    !----------------------------------------------------------------------------
    SUBROUTINE multop_field2(self, ai_i, ai_j, apr_result, ai_size)
        IMPLICIT NONE
        TYPE(FEM) :: self
        INTEGER, INTENT(IN) :: ai_i
        INTEGER, INTENT(IN) :: ai_j
        INTEGER, INTENT ( in) :: ai_size
        REAL(8), DIMENSION ( ai_size), INTENT ( out) :: apr_result

        apr_result(:) = self % opo_F(ai_i) % opr_c(:) &
        * self % opo_F(ai_j) % opr_c(:)

    end SUBROUTINE multop_field2
    !----------------------------------------------------------------------------
    SUBROUTINE addscalop_field2(self, ai_i, ar_val, apr_result, ai_size)
        IMPLICIT NONE
        TYPE(FEM) :: self
        INTEGER, INTENT(IN) :: ai_i
        REAL(8), INTENT(IN) :: ar_val
        INTEGER, INTENT ( in) :: ai_size
        REAL(8), DIMENSION ( ai_size), INTENT ( out) :: apr_result

        apr_result(:) = self % opo_F(ai_i) % opr_c(:) &
        + ar_val

    end SUBROUTINE addscalop_field2
    !----------------------------------------------------------------------------------------------
    SUBROUTINE multscalop_field2(self, ai_i, ar_val, apr_result, ai_size)
        IMPLICIT NONE
        TYPE(FEM) :: self
        INTEGER, INTENT(IN) :: ai_i
        REAL(8), INTENT(IN) :: ar_val
        INTEGER, INTENT ( in) :: ai_size
        REAL(8), DIMENSION ( ai_size), INTENT ( out) :: apr_result

        apr_result(:) = self % opo_F(ai_i) % opr_c(:) &
        * ar_val

    end SUBROUTINE multscalop_field2
    !---------------------------------------------------------------
    SUBROUTINE addop_field(self, ai_i, ai_j, ai_k)
        !> this routine implements the addition opperation :
        !> F(k) := F(i) + F(j)
        IMPLICIT NONE
        TYPE(FEM) :: self
        INTEGER :: ai_i
        INTEGER :: ai_j
        INTEGER :: ai_k

        self % opo_F(ai_k) % opr_c(:) = self % opo_F(ai_i) % opr_c(:) &
        +self % opo_F(ai_j) % opr_c(:)

    end SUBROUTINE addop_field
    !---------------------------------------------------------------
    SUBROUTINE multop_field(self, ai_i, ai_j, ai_k)
        !> this routine implements the addition opperation :
        !> F(k) := F(i) + F(j)
        IMPLICIT NONE
        TYPE(FEM) :: self
        INTEGER :: ai_i
        INTEGER :: ai_j
        INTEGER :: ai_k

        self % opo_F(ai_k) % opr_c(:) = self % opo_F(ai_i) % opr_c(:) &
        *self % opo_F(ai_j) % opr_c(:)

    end SUBROUTINE multop_field
    !---------------------------------------------------------------
    SUBROUTINE addop_field_scal(self, ai_i, ar_val, ai_k)
        !> this routine implements the addition opperation :
        !> F(k) := F(i) + F(j)
        IMPLICIT NONE
        TYPE(FEM) :: self
        INTEGER :: ai_i
        REAL(wp) :: ar_val
        INTEGER :: ai_k

        self % opo_F(ai_k) % opr_c(:) = self % opo_F(ai_i) % opr_c(:) + ar_val

    end SUBROUTINE addop_field_scal
    !---------------------------------------------------------------
    SUBROUTINE multop_field_scal(self, ai_i, ar_val, ai_k)
        !> this routine implements the addition opperation :
        !> F(k) := F(i) + F(j)
        IMPLICIT NONE
        TYPE(FEM) :: self
        INTEGER :: ai_i
        REAL(wp) :: ar_val
        INTEGER :: ai_k

        self % opo_F(ai_k) % opr_c(:) = ar_val * self % opo_F(ai_i) % opr_c(:)

    end SUBROUTINE multop_field_scal
    !---------------------------------------------------------------
    SUBROUTINE multop_field_array (self, ai_field_id, apr_val_in, ai_size_in)
        IMPLICIT NONE
        TYPE(FEM) :: self
        INTEGER :: ai_field_id
        INTEGER :: ai_size_in
        REAL(8), DIMENSION (ai_size_in), INTENT (IN) :: apr_val_in

        self % opo_F(ai_field_id) % opr_c(:) = apr_val_in * self % opo_F(ai_field_id) % opr_c(:)

    end SUBROUTINE multop_field_array

end module field_module
!**************************************************
