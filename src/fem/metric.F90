!**************************************************
!
!                   METRIC MODULE
!
!**************************************************
module metric_module
    use fem_def
    implicit none

    private :: find_grids_metric
    integer, parameter, private :: mi_dtllevel_base = 1

contains
    !----------------------------------------------------------------------------------------------
    INTEGER FUNCTION find_grids_metric (self, ai_id)
    IMPLICIT NONE
        TYPE(FEM) :: self
        integer :: ai_id
          ! LOCAL
        integer :: li_usemetric
        integer :: li_metric_id
        integer :: li_grids
#ifdef _TRACE
        CALL printlog("find_grids_metric : Start", ai_dtllevel = mi_dtllevel_base + 1)
#endif
        
    ! ...
    ! looking for the correct space_id
    ! ...
        li_grids = 0
        li_usemetric = self % opi_InfoGrids(li_grids, INFOGRIDS_USEMETRIC )
        IF (li_usemetric == 1) THEN
            li_metric_id = self % opi_InfoGrids(li_grids, INFOGRIDS_METRIC_ID )
        ELSE
            li_metric_id = -1
        END IF

        DO WHILE ( ( li_grids < self % oi_ngrids ) .AND. ( li_metric_id .NE. ai_id ) )
            li_grids = li_grids + 1
            li_usemetric = self % opi_InfoGrids(li_grids, INFOGRIDS_USEMETRIC )
            IF (li_usemetric == 1) THEN
                li_metric_id = self % opi_InfoGrids(li_grids, INFOGRIDS_METRIC_ID )
            ELSE
                li_metric_id = -1
            END IF
        END DO
    ! ...

    ! ...
    ! Check if there is no metric used
    ! ...
        IF (li_metric_id== -1) THEN
            CALL printlog("Error find_grids_metric : No Metric Is used", ai_dtllevel = mi_dtllevel_base + 1)
        END IF
    ! ...
#ifdef _TRACE
        CALL printlog("find_grids_metric : End", ai_dtllevel = mi_dtllevel_base + 1)
#endif

        find_grids_metric = li_grids
        
    END FUNCTION find_grids_metric
    !----------------------------------------------------------------------------------------------
    subroutine allocate_metric(self, ai_id)
        implicit none
        TYPE(FEM) :: self
        integer :: ai_id
        ! LOCAL
        integer :: li_nderiv
        integer :: li_space
        integer :: li_grids
        integer :: li_npts
        integer :: li_nel
        integer :: li_npatchs
#ifdef _TRACE
        CALL printlog("allocate_metric : Start", ai_dtllevel = mi_dtllevel_base + 1)
#endif

        li_grids  = find_grids_metric (self, ai_id)

        li_npatchs = self % opi_InfoGrids(li_grids, INFOGRIDS_NPATCHS)
        li_nel  = MAXVAL ( self % opi_InfoPatch(li_grids, :, INFOPATCH_NEL) )
        li_npts = MAXVAL ( self % opi_InfoPatch(li_grids, :, INFOPATCH_MAXNPTS) )
        
        li_nderiv = self % opi_Rd(li_grids)

        ALLOCATE(self % opo_metrics(ai_id) % opr_points (0:li_npatchs-1, li_nel, li_npts, 0:li_nderiv, self % opi_Rd(li_grids)) )
#ifdef _TRACE
        CALL printlog("allocate_metric : End", ai_dtllevel = mi_dtllevel_base + 1)
#endif

    end subroutine allocate_metric
    !---------------------------------------------------------------
    subroutine deallocate_metric(self, ai_id)
        implicit none
        TYPE(FEM) :: self
        INTEGER :: ai_id
        ! LOCAL
#ifdef _TRACE
        CALL printlog("deallocate_metric : Start", ai_dtllevel = mi_dtllevel_base + 1)
#endif

        DEALLOCATE(self % opo_metrics(ai_id) % opr_points)
#ifdef _TRACE
        CALL printlog("deallocate_metric : End", ai_dtllevel = mi_dtllevel_base + 1)
#endif

    end subroutine deallocate_metric
    !----------------------------------------------------------------------------------------------
    subroutine allocate_all_metrics(self)
        implicit none
        TYPE(FEM) :: self
        ! LOCAL
        integer :: li_id
#ifdef _TRACE
        CALL printlog("allocate_all_metrics : Start", ai_dtllevel = mi_dtllevel_base + 1)
#endif

        DO li_id = 0, self % oi_nmetrics - 1
            CALL allocate_metric(self, li_id)
        END DO
#ifdef _TRACE
        CALL printlog("allocate_all_metrics : End", ai_dtllevel = mi_dtllevel_base + 1)
#endif

    end subroutine allocate_all_metrics
    !----------------------------------------------------------------------------------------------
    subroutine deallocate_all_metrics(self)
        implicit none
        TYPE(FEM) :: self
        ! LOCAL
        integer :: li_id
#ifdef _TRACE
        CALL printlog("deallocate_all_metrics : Start", ai_dtllevel = mi_dtllevel_base + 1)
#endif

        DO li_id = 0, self % oi_nmetrics - 1
            CALL deallocate_metric(self, li_id)
        END DO
#ifdef _TRACE
        CALL printlog("deallocate_all_metrics : End", ai_dtllevel = mi_dtllevel_base + 1)
#endif

    end subroutine deallocate_all_metrics
end module metric_module
!**************************************************
