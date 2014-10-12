!     
! File:   connectivity.F90
! Author: root
!
! Created on December 20, 2011, 3:36 PM
!

module connectivity_module
    use tracelog_module
    use connectivities_def
    implicit none

#ifdef _DEBUG
    integer, parameter, private  :: mi_dtllevel_base = 0
#else
    integer, parameter, private  :: mi_dtllevel_base = 2
#endif

contains

    !---------------------------------------------------------------
    subroutine create_connectivity(self, ai_npatch, ai_nen, ai_nelt)
        implicit none
        type(CONNECTIVITY), intent(inout) :: self
        integer :: ai_npatch
        integer :: ai_nen
        integer :: ai_nelt
        ! LOCAL

        CALL printlog("create_connectivity : Start", ai_dtllevel = 1)

        ALLOCATE(self % opi_info(0:ai_npatch-1))
        ALLOCATE(self % opi_nen(ai_npatch))

        self % opi_nen = ai_nen

        self % oi_npatch = ai_npatch

        CALL printlog("create_connectivity : End", ai_dtllevel = 1)
    end subroutine create_connectivity
    !---------------------------------------------------------------
    subroutine set_connectivity_id(self, api_ID, ai_N)
        implicit none
        type(CONNECTIVITY), intent(inout) :: self
        integer :: ai_N
        integer, dimension(ai_N) :: api_ID
        ! LOCAL

        CALL printlog("set_connectivity_id : Start", ai_dtllevel = 1)

        ALLOCATE(self % opi_ID(ai_N))
        self % opi_ID = api_ID

        CALL printlog("set_connectivity_id : End", ai_dtllevel = 1)

    end subroutine set_connectivity_id    
    !---------------------------------------------------------------
    subroutine set_connectivity_lm(self, api_LM, ai_npatch, ai_maxnen, ai_nelt)
        implicit none
        type(CONNECTIVITY), intent(inout) :: self
        integer :: ai_npatch
        integer :: ai_maxnen
        integer :: ai_nelt
        integer, dimension(ai_npatch, ai_maxnen, ai_nelt) :: api_LM
        ! LOCAL

        CALL printlog("set_connectivity_lm : Start", ai_dtllevel = 1)

        ALLOCATE(self % opi_LM(ai_npatch, ai_maxnen, ai_nelt))
        self % opi_LM = api_LM

        CALL printlog("set_connectivity_lm : End", ai_dtllevel = 1)

    end subroutine set_connectivity_lm    
    !---------------------------------------------------------------
    subroutine set_connectivity_ien(self, api_IEN, ai_npatch, ai_maxnen, ai_nelt)
        implicit none
        type(CONNECTIVITY), intent(inout) :: self
        integer :: ai_npatch
        integer :: ai_maxnen
        integer :: ai_nelt
        integer, dimension(ai_npatch, ai_maxnen, ai_nelt) :: api_IEN
        ! LOCAL

        CALL printlog("set_connectivity_ien : Start", ai_dtllevel = 1)

        ALLOCATE(self % opi_IEN(ai_npatch, ai_maxnen, ai_nelt))
        self % opi_IEN = api_IEN

        CALL printlog("set_connectivity_ien : End", ai_dtllevel = 1)

    end subroutine set_connectivity_ien      
    !---------------------------------------------------------------
    subroutine set_real_elts(self, api_real_elts, ai_npatch, ai_nelt)
        implicit none
        type(CONNECTIVITY), intent(inout) :: self
        integer :: ai_npatch
        integer :: ai_nelt
        integer, dimension(:,:) :: api_real_elts
        ! LOCAL

        CALL printlog("set_real_elts: Start", ai_dtllevel = 1)

        ALLOCATE(self % opi_real_elts(0:ai_npatch-1, ai_nelt))

        self % opi_real_elts = api_real_elts

        CALL printlog("set_real_elts: End", ai_dtllevel = 1)
    end subroutine set_real_elts
    !---------------------------------------------------------------
    subroutine create_ID_loc(self, ai_patch, ai_n, api_id_loc)
        implicit none
        type(CONNECTIVITY), intent(inout) :: self
        integer :: ai_patch
        integer :: ai_n
        integer, dimension(ai_n) :: api_id_loc
        ! LOCAL

        CALL printlog("create_ID_loc: Start", ai_dtllevel = 1)

        ALLOCATE(self % opi_info (ai_patch) % opi_ID_loc(ai_n))

        self % opi_info (ai_patch) % opi_ID_loc = api_id_loc

        CALL printlog("create_ID_loc: End", ai_dtllevel = 1)
    end subroutine create_ID_loc
    !---------------------------------------------------------------
    subroutine free_connectivity(self)
        implicit none
        type(CONNECTIVITY), intent(inout) :: self
        ! LOCAL
        INTEGER :: li_i

        CALL printlog("free_connectivity : Start", ai_dtllevel = 1)

        DEALLOCATE(self % opi_IEN)
        DEALLOCATE(self % opi_ID)
        DEALLOCATE(self % opi_LM)
        DEALLOCATE(self % opi_real_elts)
        DEALLOCATE(self % opi_nen)
        DO li_i = 0, self % oi_npatch - 1
                DEALLOCATE(self % opi_info (li_i) % opi_ID_loc)
        END DO
        DEALLOCATE(self % opi_info)
        
        CALL printlog("free_connectivity : End", ai_dtllevel = 1)
    end subroutine free_connectivity


end module connectivity_module
