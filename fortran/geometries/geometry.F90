!
! File:   geometry.F90
! Author: ratnani
!
! Created on December 1, 2011, 10:26 AM
!

module geometry_module
    use used_precision
    use tracelog_module
    use geometries_def
    use geometry_tools
    implicit none

    !    private
    !
    !    public :: create_geometry, free_geometry, get_controlpoint, get_weight, hrefine_geometry_array &
    !    , hrefine_geometry_array_1d, hrefine_geometry_array_2d

#ifdef _DEBUG
    integer, parameter, private :: mi_dtllevel_base = 0
#else
    integer, parameter, private :: mi_dtllevel_base = 2
#endif

contains
    !---------------------------------------------------------------
    subroutine create_geometry(self, ai_dim, ai_Rd, ai_rational, api_N, api_P, apr_u, apr_P, apr_W)
        !> self routine creates and allocates memory for the object geometry
        implicit none
        type(GEOMETRY), intent(inout) :: self
        integer :: ai_dim
        integer :: ai_Rd
        integer :: ai_rational
        integer, dimension(:) :: api_N
        integer, dimension(:) :: api_P
        real(wp), dimension(:,:), optional :: apr_u
        real(wp), dimension(:,:), optional :: apr_P
        real(wp), dimension(:), optional :: apr_W
        ! LOCAL
        integer :: li_npts
        integer :: li_maxnk
        integer :: li_maxn
        integer :: li_j
        integer :: li_n1
        integer :: li_prodk
        integer :: li_prodn
        integer :: i, j, k, r 
        integer, dimension(3) :: lpi_N
        
#ifdef _TRACE
        CALL printlog("create_geometry : Start", ai_dtllevel = mi_dtllevel_base + 1)
#endif
        
        self % oi_rational = ai_rational

        lpi_N = 1 
        lpi_N(1:ai_dim) = api_N(1:ai_dim)

        self % ol_elt_allocated = .FALSE.

        !        print *,"%"
        li_npts = 1
        do li_j = 1, ai_dim
            li_npts = li_npts * api_N(li_j)
        end do
        !        print *,"%%"

!        li_maxnk = MAXVAL(api_N(:) + api_P(:) + 1)
        li_maxnk = MAXVAL(api_N(:)) + MAXVAL(api_P(:) + 1)
        li_maxn = MAXVAL(api_N(:))

!        print *, "li_maxnk=", li_maxnk
!        print *, "ai_dim=", ai_dim
!        print *,"%%%"

        ALLOCATE (self % opi_N(ai_dim))
        ALLOCATE (self % opi_P(ai_dim))
        ALLOCATE (self % opr_u(ai_dim, li_maxnk))
        ALLOCATE (self % opr_P(ai_Rd, 0:li_npts - 1))
        ALLOCATE (self % opr_W(0:li_npts - 1))
        ALLOCATE (self % opr_L(ai_dim, li_maxn))
        ALLOCATE (self % opi_ndiffu(ai_dim))
        ALLOCATE (self % opi_Real_Index(ai_dim, li_maxnk))

        ALLOCATE (self % opi_Nw(3))
        ALLOCATE (self % opr_Ww(0:lpi_N(3)-1, 0:lpi_N(2)-1, 0:lpi_N(1)-1))

        self % opi_Nw = lpi_N

        !        print *,"%%%%"
        self % oi_dim = ai_dim
        self % oi_Rd = ai_Rd
        self % oi_npts = li_npts
        !        print *,"%%%%%"
        self % opi_N(1:ai_dim) = api_N(1:ai_dim)
        !        print *,"%%%%%%"
        self % opi_P(1:ai_dim) = api_P(1:ai_dim)

        self % opr_L (:,:) = 1.0_wp
        !        print *,"%%%%%%%"

        !        print *, "SIZE(self % opr_u,1)=", SIZE(self % opr_u,1)
        !        print *, "SIZE(self % opr_u,2)=", SIZE(self % opr_u,2)
        !        print *, "SIZE(apr_u,1)=", SIZE(apr_u,1)
        !        print *, "SIZE(apr_u,2)=", SIZE(apr_u,2)
        if (present(apr_u)) then
            self % opr_u = apr_u
            CALL compute_diff_knots(self)
            CALL set_real_indices(self)
            CALL create_elt_index(self)

            li_n1 = self % opi_ndiffu(1) - 1
            li_prodk = 1
            do li_j = 1, ai_dim
                li_prodk = li_prodk * ( api_P(li_j) + 1 )
            end do
            li_prodn = 1
            do li_j = 2, ai_dim
                li_prodn = li_prodn * ( self % opi_ndiffu(li_j) - 1 )
            end do
            ALLOCATE (self % opr_X(ai_Rd, li_prodn, li_prodk * li_n1))
!            print *, 'opr_X allocated with ', li_prodk * li_n1
        end if

        if (present(apr_P)) then
            self % opr_P = apr_P
        end if

        if (present(apr_W)) then
            self % opr_W(0:li_npts - 1) = apr_W(1:li_npts)

            r = 0
            do k=0, lpi_N(3)-1 
            do j=0, lpi_N(2)-1
            do i=0, lpi_N(1)-1
               self % opr_Ww(k,j,i) = self % opr_W (r)
               r = r + 1
            end do
            end do            
            end do            
        end if

#ifdef _TRACE
        CALL printlog("create_geometry : End", ai_dtllevel = mi_dtllevel_base + 1)
#endif
    end subroutine create_geometry
    !---------------------------------------------------------------
    subroutine free_geometry(self)
        !> self routine frees and deallocates memory for the object geometry
        implicit none
        type(GEOMETRY), intent(inout) :: self
#ifdef _TRACE
        CALL printlog("free_geometry : Start", ai_dtllevel = mi_dtllevel_base + 1)
#endif

        DEALLOCATE (self % opi_N)
        DEALLOCATE (self % opi_P)
        DEALLOCATE (self % opr_u)
        DEALLOCATE (self % opr_P)
        DEALLOCATE (self % opr_W)
        DEALLOCATE (self % opr_L)
        DEALLOCATE (self % opi_ndiffu)
        DEALLOCATE (self % opi_Real_Index)
        DEALLOCATE (self % opr_X)
        DEALLOCATE (self % opi_Nw)
        DEALLOCATE (self % opr_Ww)
        IF (self % ol_elt_allocated) THEN
            DEALLOCATE (self % opi_ELT_INDEX)
            DEALLOCATE (self % opi_ELT)
        END IF
#ifdef _TRACE
        CALL printlog("free_geometry : End", ai_dtllevel = mi_dtllevel_base + 1)
#endif
    end subroutine free_geometry
    !---------------------------------------------------------------
    subroutine update_points_geometry(self, apr_P, apr_W)
        !> self routine creates and allocates memory for the object geometry
        implicit none
        type(GEOMETRY), intent(inout) :: self
        real(wp), dimension(:,:) :: apr_P
        real(wp), dimension(:) :: apr_W
        ! LOCAL
        integer :: li_npts
        integer :: li_maxnk
        integer :: li_maxn
        integer :: li_j
        integer :: li_n1
        integer :: i, j, k, r 
        
#ifdef _TRACE
        CALL printlog("update_points_geometry: Start", ai_dtllevel = mi_dtllevel_base + 1)
#endif
        li_npts = 1
        do li_j = 1, self % oi_dim 
            li_npts = li_npts * self % opi_N(li_j)
        end do

        self % opr_P = apr_P
        self % opr_W(0:li_npts - 1) = apr_W(1:li_npts)

!        r = 0
!        do k=0, self % opi_N(3)-1 
!        do j=0, self % opi_N(2)-1
!        do i=0, self % opi_N(1)-1
!           self % opr_Ww(k,j,i) = self % opr_W (r)
!           r = r + 1
!        end do
!        end do            
!        end do            

#ifdef _TRACE
        CALL printlog("update_points_geometry: End", ai_dtllevel = mi_dtllevel_base + 1)
#endif
    end subroutine update_points_geometry 
    !---------------------------------------------------------------
    subroutine copy_geometry(self, ao_geo)
        implicit none
        type(GEOMETRY), intent(in) :: self
        type(GEOMETRY), intent(inout) :: ao_geo
        ! LOCAL
        integer :: li_N
        integer :: li_P
        integer :: li_npts
        integer :: li_d
        integer :: li_Rd
        integer :: li_dim
#ifdef _TRACE
        CALL printlog("copy_geometry : Start", ai_dtllevel = mi_dtllevel_base + 5)
#endif

        !        ao_geo % oi_dim = self % oi_dim
        !        ao_geo % oi_Rd = self % oi_Rd
        !        ao_geo % oi_npts = self % oi_npts

        li_npts = self % oi_npts
        li_dim = self % oi_dim
        li_Rd = self % oi_Rd

        do li_d = 1, li_dim
            li_N = self % opi_N(li_d)
            li_P = self % opi_P(li_d)
            ao_geo % opr_u(li_d, 1: li_N + li_P + 1) = self % opr_u(li_d, 1: li_N + li_P + 1)
        end do

        ao_geo % opr_P(1:li_Rd, 0: li_npts - 1) = self % opr_P(1:li_Rd, 0: li_npts - 1)
        ao_geo % opr_W(0: li_npts - 1) = self % opr_W(0: li_npts - 1)
#ifdef _TRACE
        CALL printlog("copy_geometry : End", ai_dtllevel = mi_dtllevel_base + 5)
#endif
    end subroutine copy_geometry
    !---------------------------------------------------------------
    subroutine compute_diff_knots(self)
        !> self routine allows us to insert an array of new knots in the d-dimension of
        !> the vector knot of the geometry
        implicit none
        type(GEOMETRY), intent(inout) :: self
        ! LOCAL 
        integer :: li_d
        integer :: li_i
        integer :: li_j
        integer :: li_index
#ifdef _TRACE
        CALL printlog("compute_diff_knots : Start", ai_dtllevel = mi_dtllevel_base + 5)
#endif

        do li_d = 1, self % oi_dim
            li_index = 0
            li_i = 1
            do while (li_i <= self % opi_N(li_d) + self % opi_p(li_d))
                li_j = li_i + 1
                do while ((self % opr_u(li_d, li_j) == self % opr_u(li_d, li_i)) &
                    .AND. (li_j < li_i + self % opi_p(li_d) + 1) &
                    .AND. (li_j <= self % opi_N(li_d) + self % opi_p(li_d)))
                    li_j = li_j + 1
                end do
                li_index = li_index + 1
                li_i = li_j
            end do
            self % opi_ndiffu(li_d) = li_index
        end do
#ifdef _TRACE
        CALL printlog("compute_diff_knots : End", ai_dtllevel = mi_dtllevel_base + 5)
#endif

    end subroutine compute_diff_knots
    !---------------------------------------------------------------
    subroutine set_real_indices(self)
        !> self routine allows us to insert an array of new knots in the d-dimension of
        !> the vector knot of the geometry
        implicit none
        type(GEOMETRY), intent(inout) :: self
        ! LOCAL
        integer :: li_d
        integer :: li_i
#ifdef _TRACE
        CALL printlog("set_real_indices : Start", ai_dtllevel = mi_dtllevel_base + 5)
#endif

        do li_d = 1, self % oi_dim
            self % opi_Real_Index(li_d, 1) = 1
            do li_i = 2, self % opi_N(li_d) + self % opi_P(li_d) + 1
                if (self % opr_u(li_d, li_i) == self % opr_u(li_d, li_i - 1)) then
                    self % opi_Real_Index(li_d, li_i) = self % opi_Real_Index(li_d, li_i - 1)
                else
                    self % opi_Real_Index(li_d, li_i) = self % opi_Real_Index(li_d, li_i - 1) + 1
                end if
            end do
        end do
#ifdef _TRACE
        CALL printlog("set_real_indices : End", ai_dtllevel = mi_dtllevel_base + 5)
#endif

    end subroutine set_real_indices
    !---------------------------------------------------------------
    !> self ROUTINE CREATES AND INITIALIZES THE opi_ELT_INDEX ARRAY
    !> \todo this works only for open knot vector
    subroutine create_elt_index(self)
        implicit none
        type(GEOMETRY), intent(inout) :: self
        !LOCAL VARIABLES
#ifdef _TRACE
        CALL printlog("create_elt_index : Start", ai_dtllevel = mi_dtllevel_base + 1)
#endif
        
        select case ( self % oi_dim )
            case ( 1 )
                CALL create_elt_index_1D(self)
            case ( 2 )
                CALL create_elt_index_2D(self)
            case Default
                call printlog("Not Yet implemented", ai_dtllevel = mi_dtllevel_base +  0)
        end select            
#ifdef _TRACE
        CALL printlog("create_elt_index : End", ai_dtllevel = mi_dtllevel_base + 1)
#endif

    end subroutine create_elt_index
    !---------------------------------------------------------------
    !> self ROUTINE CREATES AND INITIALIZES THE opi_ELT_INDEX ARRAY
    !> \todo this works only for open knot vector
    subroutine create_elt_index_1D(self)
        implicit none
        type(GEOMETRY), intent(inout) :: self
        !LOCAL VARIABLES
        integer :: li_e_loc
        integer :: li_nel
        integer :: li_i
        integer :: li_e

        ! first loop to compute the number of elt and allocate memeory
        li_e_loc = 0
        li_e = 0
        do li_i = self % opi_P(1) + 1, self % opi_N(1)
            li_e = li_e + 1
            ! WE CHECK IF THE ELEMENT HAS ZERO MEASURE
            if (self % opr_u(1, li_i) /= self % opr_u(1, li_i + 1)) then

                li_e_loc = li_e_loc + 1

            end if

        end do

        li_nel = li_e_loc

        ALLOCATE(self % opi_ELT_INDEX(li_nel, self % oi_dim))
        ALLOCATE(self % opi_ELT(li_e))
        self % ol_elt_allocated = .TRUE.

        ! initialize data
        li_e_loc = 0
        li_e = 0
        do li_i = self % opi_P(1) + 1, self % opi_N(1)

            li_e = li_e + 1	
            ! WE CHECK IF THE ELEMENT HAS ZERO MEASURE
            if (self % opr_u(1, li_i) /= self % opr_u(1, li_i + 1)) then

                li_e_loc = li_e_loc + 1

                self % opi_ELT ( li_e ) = li_e_loc
                self % opi_ELT_INDEX(li_e_loc, 1) = li_i

            end if

        end do

    end subroutine create_elt_index_1D
    !---------------------------------------------------------------
    !> self ROUTINE CREATES AND INITIALIZES THE opi_ELT_INDEX ARRAY
    subroutine create_elt_index_2D(self)
        implicit none
        type(GEOMETRY), intent(inout) :: self
        !LOCAL VARIABLES
        integer :: li_e_loc
        integer :: li_e
        integer :: li_i, li_j
        integer :: li_nel

        ! first loop to compute the number of elt and allocate memeory
        li_e_loc = 0
        li_e = 0
        do li_j = self % opi_P(2) + 1, self % opi_N(2)

            do li_i = self % opi_P(1) + 1, self % opi_N(1)

                li_e = li_e + 1
                
                ! WE CHECK IF THE ELEMENT HAS ZERO MEASURE
                if ((self % opr_u(2, li_j) /= self % opr_u(2, li_j + 1))) then

                    if ((self % opr_u(1, li_i) /= self % opr_u(1, li_i + 1))) then

                        li_e_loc = li_e_loc + 1

                    end if

                end if

            end do

        end do

        li_nel = li_e_loc

        ALLOCATE(self % opi_ELT_INDEX(li_nel, self % oi_dim))
        ALLOCATE(self % opi_ELT(li_e))
        self % ol_elt_allocated = .TRUE.

        ! initialize data
        li_e_loc = 0
        li_e = 0
        do li_j = self % opi_P(2) + 1, self % opi_N(2)

            do li_i = self % opi_P(1) + 1, self % opi_N(1)

                li_e = li_e + 1

                ! WE CHECK IF THE ELEMENT HAS ZERO MEASURE
                if ((self % opr_u(2, li_j) /= self % opr_u(2, li_j + 1))) then

                    if ((self % opr_u(1, li_i) /= self % opr_u(1, li_i + 1))) then

                        li_e_loc = li_e_loc + 1

                        self % opi_ELT ( li_e ) = li_e_loc

                        self % opi_ELT_INDEX(li_e_loc, 1) = li_i
                        self % opi_ELT_INDEX(li_e_loc, 2) = li_j

                    end if

                end if

            end do

        end do

    end subroutine create_elt_index_2D
    !---------------------------------------------------------------
    subroutine print_geometry(self, ai_dtllevel)
        implicit none
        type(GEOMETRY) :: self
        integer, optional :: ai_dtllevel
        ! LOCAL VARAIABLES
        integer :: li_i
        integer :: li_d
        integer :: li_detail
#ifdef _TRACE
        CALL printlog("print_geometry : Start", ai_dtllevel = mi_dtllevel_base + 1)
#endif

        li_detail = 4
        if (present(ai_dtllevel)) then
            li_detail = ai_dtllevel
        end if

        if (li_detail == 0) then
            return
        end if

        if (li_detail >= 1) then
            print*, '======= params ======='
            print*, 'dim =', self % oi_dim
            print*, 'Rd =', self % oi_Rd
            print*, 'N =', self % opi_N(:)
            print*, 'P =', self % opi_P(:)
        end if

        if (li_detail >= 2) then
            print*, '======= knots ======='
            do li_d = 1, self % oi_dim
                print *, "-> dim = ", li_d
                print*, self % opr_u(li_d, 1:self % opi_N(li_d) + self % opi_P(li_d) + 1)
            end do
        end if

        if (li_detail >= 3) then
            print*, '======= Points ======='
            do li_d = 1, self % oi_Rd
                print *, "-> dim = ", li_d
                print*, self % opr_P(li_d,:)
            end do
        end if

        if (li_detail >= 4) then
            print*, '======= Weights ======='
            print*, self % opr_W(:)
        end if
#ifdef _TRACE
        CALL printlog("print_geometry : End", ai_dtllevel = mi_dtllevel_base + 1)
#endif

    end subroutine print_geometry

end module geometry_module
