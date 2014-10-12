!**************************************************
!
!                   NORM MODULE
!
!**************************************************
module norm_module
    use fem_def
    use spaces_module
    implicit none

    integer, parameter, private :: mi_dtllevel_base = 1

contains

    !---------------------------------------------------------------
    subroutine SET_INFONORM(self, ai_id, ai_param, ai_val)
        implicit none
        TYPE(FEM) :: self
        INTEGER :: ai_id
        INTEGER :: ai_param
        INTEGER :: ai_val

        self % opi_InfoNorm(ai_id, ai_param) = ai_val

    end subroutine SET_INFONORM
    !---------------------------------------------------------------
    subroutine allocate_norms(self)
        implicit none
        TYPE(FEM) :: self
        ! LOCAL
        INTEGER :: li_Norm
        INTEGER :: li_field
        INTEGER :: li_space
        INTEGER :: li_grids
        INTEGER :: li_npatchs
        INTEGER :: li_maxnel
#ifdef _TRACE
        CALL printlog("allocate_norms : Start", ai_dtllevel = mi_dtllevel_base + 1)
#endif

        ALLOCATE(self % opo_N(0:self % oi_nNorms - 1))

        DO li_Norm = 0, self % oi_nNorms - 1
            li_field = self % opi_InfoNorm(li_Norm, INFONORM_FIELD)
            li_space = self % opi_InfoField(li_field, INFOFIELD_SPACE)
            li_grids = self % opi_InfoSpace(li_space, INFOSPACE_GRIDS)
            li_npatchs = self % opi_InfoGrids(li_grids, INFOGRIDS_NPATCHS)
            
            CALL get_maxnel_space(self, li_space, li_maxnel)

            ALLOCATE(self % opo_N(li_Norm) % opr_values(0:li_npatchs-1, li_maxnel))
            self % opo_N(li_Norm) % opr_values (:,:) = 0.0_wp

        END DO
#ifdef _TRACE
        CALL printlog("allocate_norms : End", ai_dtllevel = mi_dtllevel_base + 1)
#endif

    end subroutine allocate_norms
    !---------------------------------------------------------------
    subroutine deallocate_norms(self)
        implicit none
        TYPE(FEM) :: self
        ! LOCAL
        INTEGER :: li_Norm
#ifdef _TRACE
        CALL printlog("deallocate_norms : Start", ai_dtllevel = mi_dtllevel_base + 1)
#endif

        DO li_Norm = 0, self % oi_nNorms - 1
            DEALLOCATE(self % opo_N(li_Norm) % opr_values)
        END DO

        DEALLOCATE(self % opo_N)
#ifdef _TRACE
        CALL printlog("deallocate_norms : End", ai_dtllevel = mi_dtllevel_base + 1)
#endif

    end subroutine deallocate_norms
    !----------------------------------------------------------------------------------------------
    subroutine getGlobalNorm_fem(self, ai_norm, ar_result)
        implicit none
        TYPE(FEM) :: self
        integer, intent(in) :: ai_norm
        real(wp), intent ( out) :: ar_result
        ! LOCAL
        integer :: li_patch
        integer :: li_elt
        real(wp) :: lr_norm
        INTEGER :: li_field
        INTEGER :: li_space
        INTEGER :: li_grids
        INTEGER :: li_npatchs
        INTEGER :: li_nel
        
        li_field = self % opi_InfoNorm(ai_norm, INFONORM_FIELD)
        li_space = self % opi_InfoField(li_field, INFOFIELD_SPACE)
        li_grids = self % opi_InfoSpace(li_space, INFOSPACE_GRIDS)
        li_npatchs = self % opi_InfoGrids(li_grids, INFOGRIDS_NPATCHS)

        lr_norm = 0.0_wp
        DO li_patch = 0, li_npatchs - 1
            li_nel = self % opi_InfoPatch(li_grids, li_patch, INFOPATCH_NEL)
            DO li_elt = 1,li_nel
                lr_norm = lr_norm + self % opo_N (ai_norm) % opr_values (li_patch, li_elt)
            END DO
        END DO

        lr_norm = SQRT (lr_norm)

        ar_result = lr_norm

    end subroutine getGlobalNorm_fem
    !----------------------------------------------------------------------------------------------
    subroutine getPatchNorm_fem(self, ai_norm, apr_result, ai_npatchs)
        implicit none
        TYPE(FEM) :: self
        integer, intent(in) :: ai_norm
        integer, intent ( in) :: ai_npatchs
        real(wp), dimension ( ai_npatchs), intent ( out) :: apr_result
        ! LOCAL
        integer :: li_elt
        INTEGER :: li_patch
        real(wp) :: lr_norm
        INTEGER :: li_field
        INTEGER :: li_space
        INTEGER :: li_grids
        INTEGER :: li_nel
        
        li_field = self % opi_InfoNorm(ai_norm, INFONORM_FIELD)
        li_space = self % opi_InfoField(li_field, INFOFIELD_SPACE)
        li_grids = self % opi_InfoSpace(li_space, INFOSPACE_GRIDS)

        DO li_patch = 0, ai_npatchs - 1
            li_nel = self % opi_InfoPatch(li_grids, li_patch, INFOPATCH_NEL)
            
            lr_norm = 0.0_wp
            DO li_elt = 1,li_nel
                lr_norm = lr_norm + self % opo_N (ai_norm) % opr_values (li_patch, li_elt)
            END DO
            lr_norm = SQRT (lr_norm)

            apr_result (li_patch+1) = lr_norm
        END DO

    end subroutine getPatchNorm_fem
    !----------------------------------------------------------------------------------------------
    subroutine getElementNorm_fem(self, ai_norm, ai_patch, apr_result, ai_nel)
        implicit none
        TYPE(FEM) :: self
        integer, intent(in) :: ai_norm
        integer, intent(in) :: ai_patch
        integer, intent(in) :: ai_nel
        real(wp), dimension (ai_nel), intent ( out) :: apr_result
        ! LOCAL

        apr_result(1:ai_nel) = SQRT(self % opo_N (ai_norm) % opr_values (ai_patch, 1:ai_nel))

    end subroutine getElementNorm_fem

end module norm_module
!**************************************************
