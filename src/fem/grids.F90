!**************************************************
!
!                   GRIDS MODULE
!
!**************************************************
module grids_module
!    use grids_def
!    USE grids, ONLY : create_elements, free_elements
    use fem_def
    implicit none

    integer, parameter, private :: mi_dtllevel_base = 1

contains

    !---------------------------------------------------------------
    subroutine SET_INFOGRIDS(self, ai_id, ai_param, ai_val)
        implicit none
        TYPE(FEM) :: self
        INTEGER :: ai_id
        INTEGER :: ai_param
        INTEGER :: ai_val

        self % opi_InfoGrids(ai_id, ai_param) = ai_val
    end subroutine SET_INFOGRIDS
    !----------------------------------------------------------------------------------------------
    !> \todo not finished yet
    subroutine allocate_grids(self, ai_id)
        implicit none
        TYPE(FEM) :: self
        integer :: ai_id
        ! LOCAL
        integer :: li_npatchs
        integer :: li_dim
        integer :: li_patch
        integer :: li_nel
        integer :: li_npts
        integer :: li_dirnpts
        integer :: li_dof
        integer :: li_nfields
        integer :: li_nmatrices
        integer :: li_nnorms
        integer :: li_tensorlevel
        integer, dimension(:), pointer :: lpi_npts
#ifdef _TRACE
        CALL printlog("allocate_grids : Start", ai_dtllevel = mi_dtllevel_base + 1)
#endif

        li_npatchs = self % opi_InfoGrids(ai_id, INFOGRIDS_NPATCHS)
        li_dof = self % opi_InfoGrids(ai_id, INFOGRIDS_DOF)
        li_nfields = self % opi_InfoGrids(ai_id, INFOGRIDS_NFIELDS)
        li_nmatrices = self % opi_InfoGrids(ai_id, INFOGRIDS_NMATRICES)
        li_nnorms = self % opi_InfoGrids(ai_id, INFOGRIDS_NNORMS)

#ifdef _DEBUG
call concatmsg("ai_id = ", ai_dtllevel = mi_dtllevel_base + 3)
call concatmsg(ai_id, ai_dtllevel = mi_dtllevel_base + 3)

call concatmsg("li_npatchs = ", ai_dtllevel = mi_dtllevel_base + 3)
call concatmsg(li_npatchs, ai_dtllevel = mi_dtllevel_base + 3)

call concatmsg("li_dof = ", ai_dtllevel = mi_dtllevel_base + 3)
call concatmsg(li_dof, ai_dtllevel = mi_dtllevel_base + 3)

call printmsg(ai_dtllevel = mi_dtllevel_base + 3)
#endif

        self % opo_grids(ai_id) % oi_nPatchs = li_npatchs
        ALLOCATE ( self % opo_grids(ai_id) % opo_grid(0:li_npatchs - 1))
!        print *, "%"
        ! ...
        ! creation of the grids
        ! ...
        li_dim = self % opi_dim(ai_id)
        ALLOCATE(lpi_npts(li_dim))
        DO li_patch = 0, li_npatchs - 1
!            print *, "%%"
            li_nel = self % opi_InfoPatch(ai_id, li_patch, INFOPATCH_NEL)
!            print *, 'ai_id, li_patch, li_nel = ', ai_id, li_patch, li_nel
            li_npts = self % opi_InfoPatch(ai_id, li_patch, INFOPATCH_MAXNPTS)
            li_tensorlevel = self % opi_InfoPatch(ai_id, li_patch, INFOPATCH_TENSOR)
            li_dirnpts = self % opi_InfoPatch(ai_id, li_patch, INFOPATCH_DIRMAXNPTS)
            lpi_npts    = li_npts

#ifdef _DEBUG
call concatmsg("li_dim = ", ai_dtllevel = mi_dtllevel_base + 3)
call concatmsg(li_dim, ai_dtllevel = mi_dtllevel_base + 3)

call concatmsg("li_dof = ", ai_dtllevel = mi_dtllevel_base + 3)
call concatmsg(li_dof, ai_dtllevel = mi_dtllevel_base + 3)

call concatmsg("li_nel = ", ai_dtllevel = mi_dtllevel_base + 3)
call concatmsg(li_nel, ai_dtllevel = mi_dtllevel_base + 3)

call concatmsg("li_npts = ", ai_dtllevel = mi_dtllevel_base + 3)
call concatmsg(li_npts, ai_dtllevel = mi_dtllevel_base + 3)

call concatmsg("li_dirnpts = ", ai_dtllevel = mi_dtllevel_base + 3)
call concatmsg(li_dirnpts, ai_dtllevel = mi_dtllevel_base + 3)

call concatmsg("li_tensorlevel = ", ai_dtllevel = mi_dtllevel_base + 3)
call concatmsg(li_tensorlevel, ai_dtllevel = mi_dtllevel_base + 3)

call concatmsg("li_nfields = ", ai_dtllevel = mi_dtllevel_base + 3)
call concatmsg(li_nfields, ai_dtllevel = mi_dtllevel_base + 3)

call concatmsg("li_nmatrices = ", ai_dtllevel = mi_dtllevel_base + 3)
call concatmsg(li_nmatrices, ai_dtllevel = mi_dtllevel_base + 3)

call concatmsg("li_nnorms = ", ai_dtllevel = mi_dtllevel_base + 3)
call concatmsg(li_nnorms, ai_dtllevel = mi_dtllevel_base + 3)

call concatmsg("self % oi_maxnparams_matrices = ", ai_dtllevel = mi_dtllevel_base + 3)
call concatmsg(self % oi_maxnparams_matrices, ai_dtllevel = mi_dtllevel_base + 3)

call concatmsg("self % oi_maxnparams_fields = ", ai_dtllevel = mi_dtllevel_base + 3)
call concatmsg(self % oi_maxnparams_fields, ai_dtllevel = mi_dtllevel_base + 3)

call concatmsg("self % oi_maxnparams_norms = ", ai_dtllevel = mi_dtllevel_base + 3)
call concatmsg(self % oi_maxnparams_norms, ai_dtllevel = mi_dtllevel_base + 3)

call printmsg(ai_dtllevel = mi_dtllevel_base + 3)
#endif

        ! *************************************************
        !                 TENSOR-LEVEL   1
        ! *************************************************
            IF (li_tensorlevel==1 ) THEN
                CALL create_elements(self % opo_grids(ai_id) % opo_grid(li_patch) &
                , li_dim, li_nel    &
    !            , ai_dof = li_dof    &
                , li_dof    &
    !            , ai_nfields   = li_nfields    &
                , li_nfields    &
    !            , ai_nmatrices = li_nmatrices &
                , li_nmatrices &
    !            , ai_nnorms = li_nnorms &
                , li_nnorms &
    !            , ai_npts = li_npts &
                , li_npts &
    !            , ai_dirnpts = li_dirnpts   &
                , li_dirnpts   &
    !            , al_tensor = ll_tensor &
                , li_tensorlevel &
    !            , ai_maxnparams_matrices = self % oi_maxnparams_matrices &
                , self % oi_maxnparams_operators &
        !            , ai_maxnparams_fields = self % oi_maxnparams_fields & ! TODO il faut recalculer le nparam et reallouer si necessaire
    !            , ai_maxnparams_fields = MAX(4, self % oi_maxnparams_fields) &
                , MAX(4, self % oi_maxnparams_fields) &
    !            , ai_maxnparams_norms = self % oi_maxnparams_norms )
                , self % oi_maxnparams_norms )
            END IF
        ! *************************************************

        ! *************************************************
        !                 TENSOR-LEVEL   2
        ! *************************************************
!            IF (li_tensorlevel==2 ) THEN
!                CALL create_elements(self % opo_grids(ai_id) % opo_grid(li_patch) &
!                , li_dim, li_nel    &
!    !            , ai_dof = li_dof    &
!                , li_dof    &
!    !            , ai_nfields   = li_nfields    &
!                , li_nfields    &
!    !            , ai_nmatrices = li_nmatrices &
!                , li_nmatrices &
!    !            , ai_nnorms = li_nnorms &
!                , li_nnorms &
!    !            , api_npts = lpi_npts &
!!                , lpi_npts &
!    !            , ai_npts = li_npts &
!                , li_npts &
!    !            , ai_dirnpts = li_dirnpts   &
!                , li_dirnpts   &
!    !            , al_tensor = ll_tensor &
!                , li_tensorlevel &
!    !            , ai_maxnparams_matrices = self % oi_maxnparams_matrices &
!                , self % oi_maxnparams_operators &
!        !            , ai_maxnparams_fields = self % oi_maxnparams_fields & ! TODO il faut recalculer le nparam et reallouer si necessaire
!    !            , ai_maxnparams_fields = MAX(4, self % oi_maxnparams_fields) &
!                , MAX(4, self % oi_maxnparams_fields) &
!    !            , ai_maxnparams_norms = self % oi_maxnparams_norms )
!                , self % oi_maxnparams_norms )
!            END IF
        END DO

        DEALLOCATE(lpi_npts)

        ! ...
#ifdef _TRACE
        CALL printlog("allocate_grids : End", ai_dtllevel = mi_dtllevel_base + 1)
#endif

    end subroutine allocate_grids
    !---------------------------------------------------------------
    subroutine deallocate_grids(self, ai_id)
        implicit none
        TYPE(FEM) :: self
        INTEGER :: ai_id
        ! LOCAL
        integer :: li_id
#ifdef _TRACE
        CALL printlog("deallocate_grids : Start", ai_dtllevel = mi_dtllevel_base + 1)
#endif

        ! ...
        DO li_id = 0, self % opo_grids(ai_id) % oi_nPatchs - 1
            CALL free_elements(self % opo_grids(ai_id) % opo_grid(li_id))
        END DO
        ! ...

        DEALLOCATE ( self % opo_grids(ai_id) % opo_grid)
#ifdef _TRACE
        CALL printlog("deallocate_grids : End", ai_dtllevel = mi_dtllevel_base + 1)
#endif

    end subroutine deallocate_grids
    !---------------------------------------------------------------
    subroutine allocate_all_grids(self)
        implicit none
        TYPE(FEM) :: self
        ! LOCAL
        INTEGER :: li_id
#ifdef _TRACE
        CALL printlog("allocate_all_grids : Start", ai_dtllevel = mi_dtllevel_base + 1)
#endif

        ALLOCATE(self % opo_grids(0:self % oi_nGrids - 1))

        DO li_id = 0, self % oi_nGrids - 1
            CALL allocate_grids(self, li_id)
        END DO
#ifdef _TRACE
        CALL printlog("allocate_all_grids : End", ai_dtllevel = mi_dtllevel_base + 1)
#endif

    end subroutine allocate_all_grids
    !---------------------------------------------------------------
    subroutine deallocate_all_grids(self)
        implicit none
        TYPE(FEM) :: self
        ! LOCAL
        INTEGER :: li_id
#ifdef _TRACE
        CALL printlog("deallocate_all_grids : Start", ai_dtllevel = mi_dtllevel_base + 1)
#endif

        DO li_id = 0, self % oi_nGrids - 1
            CALL deallocate_grids(self, li_id)
        END DO
        DEALLOCATE ( self % opo_grids)
#ifdef _TRACE
        CALL printlog("deallocate_all_grids : End", ai_dtllevel = mi_dtllevel_base + 1)
#endif

    end subroutine deallocate_all_grids

end module grids_module
!**************************************************
