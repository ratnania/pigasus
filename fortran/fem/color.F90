!**************************************************
!
!                   COLOR MODULE
!
!**************************************************
module color_module
    use fem_def
    implicit none

    integer, parameter, private :: mi_dtllevel_base = 1

contains

    !---------------------------------------------------------------
    subroutine SET_INFOCOLOR(self, ai_id, ai_param, ai_val)
        implicit none
        TYPE(FEM) :: self
        INTEGER :: ai_id
        INTEGER :: ai_param
        INTEGER :: ai_val

        self % opi_InfoColor(ai_id, ai_param) = ai_val
    end subroutine SET_INFOCOLOR    
    !---------------------------------------------------------------
    subroutine save_color_toassembly(self, ai_color)
        IMPLICIT NONE
        TYPE(FEM) :: self 
        INTEGER, intent(in)  :: ai_color        
        ! LOCAL

        CALL printlog("save_color_toassembly: Start", ai_dtllevel = mi_dtllevel_base  + 2)

        self % opo_colors (ai_color) % opi_objects_toassembly_tmp = self % opo_colors (ai_color) % opi_objects_toassembly

        CALL printlog("save_color_toassembly : End", ai_dtllevel = mi_dtllevel_base  + 2)

    end subroutine save_color_toassembly
    !---------------------------------------------------------------
    subroutine load_color_toassembly(self, ai_color)
        IMPLICIT NONE
        TYPE(FEM) :: self 
        INTEGER, intent(in)  :: ai_color        
        ! LOCAL

        CALL printlog("load_color_toassembly: Start", ai_dtllevel = mi_dtllevel_base  + 2)

        self % opo_colors (ai_color) % opi_objects_toassembly  = self % opo_colors (ai_color) % opi_objects_toassembly_tmp

        CALL printlog("load_color_toassembly: End", ai_dtllevel = mi_dtllevel_base  + 2)

    end subroutine load_color_toassembly
    !--------------------------------------------------------------
    !> \todo not finished yet
    subroutine allocate_color(self, ai_id)
        implicit none
        TYPE(FEM) :: self
        !> COLOR REFERENCE IN THE DICTIONNARY
        integer :: ai_id
        ! LOCAL
        integer :: li_maxaddto

#ifdef _TRACE
        CALL printlog("allocate_color : Start", ai_dtllevel = mi_dtllevel_base + 1)
#endif

        li_maxaddto     = self % oi_maxcoloraddto

        ALLOCATE(self % opo_colors (ai_id) % opi_objects_toassembly(0:li_maxaddto))
        ALLOCATE(self % opo_colors (ai_id) % opi_objects_toassembly_tmp(0:li_maxaddto))

        self % opo_colors (ai_id) % opi_objects_toassembly      = 0
        self % opo_colors (ai_id) % opi_objects_toassembly_tmp  = 0

#ifdef _TRACE
        CALL printlog("allocate_color : End", ai_dtllevel = mi_dtllevel_base + 1)
#endif

    end subroutine allocate_color
    !--------------------------------------------------------------
    !> \todo not finished yet
    subroutine deallocate_color(self, ai_id)
        implicit none
        TYPE(FEM) :: self
        !> COLOR REFERENCE IN THE DICTIONNARY
        integer :: ai_id
        ! LOCAL

#ifdef _TRACE
        CALL printlog("deallocate_color : Start", ai_dtllevel = mi_dtllevel_base + 1)
#endif

        DEALLOCATE(self % opo_colors (ai_id) % opi_objects_toassembly)
        DEALLOCATE(self % opo_colors (ai_id) % opi_objects_toassembly_tmp)

#ifdef _TRACE
        CALL printlog("deallocate_color : End", ai_dtllevel = mi_dtllevel_base + 1)
#endif

    end subroutine deallocate_color
    !---------------------------------------------------------------
    subroutine allocate_colors(self)
        implicit none
        TYPE(FEM) :: self
        ! LOCAL
        INTEGER :: li_id
        INTEGER :: ierr 

#ifdef _TRACE
        CALL printlog("allocate_colors : Start", ai_dtllevel = mi_dtllevel_base + 1)
#endif

        PRINT *, 'self % oi_maxnColors = ', self % oi_maxnColors

        DO li_id = 0, self % oi_maxnColors - 1
        CALL allocate_color(self, li_id)
        END DO

#ifdef _TRACE
        CALL printlog("allocate_colors : End", ai_dtllevel = mi_dtllevel_base + 1)
#endif

    end subroutine allocate_colors
    !---------------------------------------------------------------
    subroutine deallocate_colors(self)
        implicit none
        TYPE(FEM) :: self
        ! LOCAL
        INTEGER :: li_id
        INTEGER :: ierr

#ifdef _TRACE
        CALL printlog("deallocate_colors : Start", ai_dtllevel = mi_dtllevel_base + 1)
#endif

        DO li_id = 0, self % oi_maxnColors - 1
        CALL deallocate_color(self, li_id)
        END DO

#ifdef _TRACE
        CALL printlog("deallocate_colors : End", ai_dtllevel = mi_dtllevel_base + 1)
#endif

    end subroutine deallocate_colors


end module color_module
!**************************************************
