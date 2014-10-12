!
! File:   bbox.F90
! Author: ratnani
!
! Created on December 1, 2011, 10:26 AM
!

module bbox_module
    use tracelog_module
    use bbox_def
!    use assl_module
    implicit none

    interface create_bbox
            module procedure create_bbox_t1 !, create_bbox_t2
    end interface

#ifdef _DEBUG
    integer, parameter, private  :: mi_dtllevel_base = 0
#else
    integer, parameter, private  :: mi_dtllevel_base = 2
#endif
    
contains
    !---------------------------------------------------------------
    subroutine create_bbox_t1(self, ai_dim, ai_nderiv, ai_nen, ai_npts, ai_maxnpts, ai_maxp)
        implicit none
        type(BBOX), intent(inout) :: self
        integer :: ai_dim
        integer :: ai_nderiv
        integer :: ai_nen
        !> the number of points in the element
        !> we must give the max-value among all elements
        integer :: ai_npts
        !> the maximum number of points among all directions and all elements
        integer, optional :: ai_maxnpts
        integer, optional :: ai_maxp
        ! LOCAL
        
#ifdef _TRACE
        CALL printlog("create_bbox_t1 : Start", ai_dtllevel = mi_dtllevel_base + 1)
#endif

        self % oi_nderiv = ai_nderiv
        self % oi_nen = ai_nen
        self % oi_npts = ai_npts
        self % oi_dim = ai_dim
        CALL get_deriv_code(ai_dim, ai_nderiv, self % oi_nderiv_code)

        self % oi_tensorlevel = 1

        if ( present(ai_maxnpts)) then
            self % oi_maxnpts = ai_maxnpts
        end if

        if ( present(ai_maxp)) then
            self % oi_maxp = ai_maxp
        end if

        ALLOCATE(self % opi_leftmk(ai_dim))
        ALLOCATE(self % opr_dB(ai_dim, 0:ai_nderiv, 0:ai_maxp, 1:ai_maxnpts))
        ALLOCATE(self % opr_dBasis(0:self % oi_nderiv_code, 0:ai_nen-1, 1:ai_npts))
        
#ifdef _TRACE
        CALL printlog("create_bbox_t1 : End", ai_dtllevel = mi_dtllevel_base + 1)
#endif
        
    end subroutine create_bbox_t1
    !---------------------------------------------------------------
!    subroutine create_bbox_t2(self, ai_dim, ai_nderiv, ai_n, ai_nel, api_k, api_g)
!        implicit none
!        type(BBOX), intent(inout) :: self
!        integer :: ai_dim
!        integer :: ai_nderiv
!        integer :: ai_n
!        integer :: ai_nel
!        integer, dimension(:), pointer :: api_k
!        integer, dimension(:), pointer :: api_g
!        ! LOCAL
!        integer :: li_d, li_x
!        integer, dimension(2) :: lpi_k
!        integer, dimension(2) :: lpi_g
!
!#ifdef _TRACE
!        CALL printlog("create_bbox_t2 : Start", ai_dtllevel = mi_dtllevel_base + 1)
!#endif
!
!        self % oi_tensorlevel = 2
!        self % oi_dim = ai_dim
!        self % oi_nderiv = ai_nderiv
!
!        self % oi_nen = 1
!        self % oi_npts = 1
!
!        DO li_d = 1, ai_dim
!            self % oi_nen = self % oi_nen * api_k(li_d)
!            self % oi_npts = self % oi_npts * api_g(li_d)
!        END DO
!
!        CALL get_deriv_code(ai_dim, ai_nderiv, self % oi_nderiv_code)
!
!        ! ...
!        ! create assl for geometry computation
!        ! ...
!        lpi_k(1) = api_k(1)
!        lpi_g(1) = api_g(1)
!        lpi_k (2:ai_dim) = 1
!        lpi_g (2:ai_dim) = 1
!        DO li_d = 2, ai_dim
!            lpi_k(2) = lpi_k(2) * api_k(li_d)
!            lpi_g(2) = lpi_g(2) * api_g(li_d)
!        END DO
!        call assl_create(self % oo_asslx, ai_n &
!                            , lpi_g(1), lpi_k(1) &
!                            , lpi_g(2), lpi_k(2) )
!!                            , lpi_k(1), lpi_g(1) &
!!                            , lpi_k(2), lpi_g(2) )
!        ! ...
!
!        ! ...
!        ! allocate local basis
!        ! ...
!        ! il faudrait utiliser le nderiv_code pour l'hyperplan
!        ALLOCATE(self % opr_lBasis  (0:self % oi_nderiv_code, lpi_k(2) * lpi_g(2)              , ai_nel))
!        ALLOCATE(self % opr_lBasis_t(0:self % oi_nderiv_code, lpi_k(2) * lpi_k(2) * lpi_g(2), ai_nel))
!        self % opr_lBasis = 0.0
!        self % opr_lBasis_t = 0.0
!        ! ...
!
!        ! ...
!        ! create assl for the assembling process
!        ! ...
!        lpi_k = lpi_k * lpi_k
!        call assl_create(self % oo_assl, ai_n &
!!                            , lpi_g(1), lpi_k(1) &
!!                            , lpi_g(2), lpi_k(2) )
!                            , lpi_k(1), lpi_g(1) &
!                            , lpi_k(2), lpi_g(2) )
!
!        ! ...
!
!!        print *, 'asslx % annz =', self % oo_asslx % annz
!!        print *, 'asslx % rnnz =', self % oo_asslx % rnnz
!!        print *, 'asslx % bnnz =', self % oo_asslx % bnnz
!!        print *, 'asslx % arbnnz =', self % oo_asslx % arbnnz
!!
!!        print *, 'assl % annz =', self % oo_assl % annz
!!        print *, 'assl % rnnz =', self % oo_assl % rnnz
!!        print *, 'assl % bnnz =', self % oo_assl % bnnz
!!        print *, 'assl % arbnnz =', self % oo_assl % arbnnz
!        ! ...
!        ! allocate global basis
!        ! ...
!#ifdef _TRACE
!        CALL printlog("TODO : a revoir la taille de l allocation self % oo_asslx % annz", ai_dtllevel = mi_dtllevel_base + 1)
!#endif
!
!        ALLOCATE(self % opr_gBasis  (0:self % oi_nderiv_code, assl_get_annz(self % oo_asslx)))
!        self % opr_gBasis = 0.0
!        ALLOCATE(self % opr_gBasis_t(0:self % oi_nderiv_code, assl_get_annz(self % oo_assl)))
!        self % opr_gBasis_t = 0.0
!        ! ...
!        
!#ifdef _TRACE
!        CALL printlog("create_bbox_t2 : End", ai_dtllevel = mi_dtllevel_base + 1)
!#endif
!
!    end subroutine create_bbox_t2
    !---------------------------------------------------------------
    subroutine free_bbox(self)
        implicit none
        type(BBOX), intent(inout) :: self
        ! local
        integer :: li_d, li_x        
        
#ifdef _TRACE
        CALL printlog("free_bbox : Start", ai_dtllevel = mi_dtllevel_base + 1)
#endif
        IF (self % oi_tensorlevel == 1) THEN
            DEALLOCATE(self % opi_leftmk)
            DEALLOCATE(self % opr_dB)
            DEALLOCATE(self % opr_dBasis)
        END IF

!        IF (self % oi_tensorlevel == 2) THEN
!            DEALLOCATE(self % opr_lBasis)
!            DEALLOCATE(self % opr_lBasis_t)
!            DEALLOCATE(self % opr_gBasis)
!            DEALLOCATE(self % opr_gBasis_t)
!            call assl_free(self % oo_asslx)
!            call assl_free(self % oo_assl)
!!            print *, 'free_bbox : probleme pour desallouer assl'
!        END IF

#ifdef _TRACE
        CALL printlog("free_bbox : End", ai_dtllevel = mi_dtllevel_base + 1)
#endif
    end subroutine free_bbox
    !---------------------------------------------------------------
    subroutine get_deriv_code(ai_dim, ai_nderiv, ai_nderiv_code)
        implicit none
        integer, intent(in) :: ai_dim
        integer, intent(in) :: ai_nderiv
        integer, intent(inout) :: ai_nderiv_code

        if (ai_dim == 1) then
            ai_nderiv_code = ai_nderiv
            return
        end if

        if (ai_dim == 2) then
            select case ( ai_nderiv)
                case ( 0 )
                    ai_nderiv_code = 0; return
                case ( 1 )
                    ai_nderiv_code = 2; return
                case ( 2 )
                    ai_nderiv_code = 5; return
                case Default
                    print*,"get_deriv_code: Type code Not Yet implemented"
                    return
            end select
        end if

        if (ai_dim == 3) then
            select case ( ai_nderiv)
                case ( 0 )
                    ai_nderiv_code = 0; return
                case ( 1 )
                    ai_nderiv_code = 3; return
                case ( 2 )
                    ai_nderiv_code = 9; return
                case Default
                    print*,"get_deriv_code: Type code Not Yet implemented"
                    return
            end select
        end if

    end subroutine get_deriv_code
end module bbox_module
