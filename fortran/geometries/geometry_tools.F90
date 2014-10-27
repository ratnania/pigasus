!     
! File:   geometry_tools.F90
! Author: ratnani
!
! Created on December 5, 2011, 2:11 PM
!

module geometry_tools
    use bsp
    use used_precision
    use tracelog_module
    use geometries_def
    implicit none

#ifdef _DEBUG
    integer, parameter, private  :: mi_dtllevel_base = 0
#else
    integer, parameter, private  :: mi_dtllevel_base = 2
#endif

contains

subroutine interv ( xt, lxt, x, left, mflag )

!*****************************************************************************80
!
!! INTERV brackets a real value in an ascending vector of values.
!
!  Discussion:
!
!    The XT array is a set of increasing values.  The goal of the routine
!    is to determine the largest index I so that XT(I) <= X.
!
!    The routine is designed to be efficient in the common situation
!    that it is called repeatedly, with X taken from an increasing
!    or decreasing sequence.
!
!    This will happen when a piecewise polynomial is to be graphed.
!    The first guess for LEFT is therefore taken to be the value
!    returned at the previous call and stored in the local variable ILO.
!
!    A first check ascertains that ILO < LXT.  This is necessary
!    since the present call may have nothing to do with the previous
!    call.  Then, if
!
!      XT(ILO) <= X < XT(ILO+1),
!
!    we set LEFT = ILO and are done after just three comparisons.
!
!    Otherwise, we repeatedly double the difference ISTEP = IHI - ILO
!    while also moving ILO and IHI in the direction of X, until
!
!      XT(ILO) <= X < XT(IHI)
!
!    after which we use bisection to get, in addition, ILO + 1 = IHI.
!    The value LEFT = ILO is then returned.
!
!  Modified:
!
!    14 February 2007
!
!  Author:
!
!    Carl DeBoor
!
!  Reference:
!
!    Carl DeBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663.
!
!  Parameters:
!
!    Input, real(wp) XT(LXT), a nondecreasing sequence of values.
!
!    Input, integer LXT, the dimension of XT.
!
!    Input, real(wp) X, the point whose location with
!    respect to the sequence XT is to be determined.
!
!    Output, integer LEFT, the index of the bracketing value:
!      1     if             X  <  XT(1)
!      I     if   XT(I)  <= X  < XT(I+1)
!      LXT   if  XT(LXT) <= X
!
!    Output, integer MFLAG, indicates whether X lies within the
!    range of the data.
!    -1:            X  <  XT(1)
!     0: XT(I)   <= X  < XT(I+1)
!    +1: XT(LXT) <= X
!
  implicit none

  integer lxt

  integer left
  integer mflag
  integer ihi
  integer, save :: ilo = 1
  integer istep
  integer middle
  real(wp) x
  real(wp) xt(lxt)

  ihi = ilo + 1

  if ( lxt <= ihi ) then

    if ( xt(lxt) <= x ) then
      go to 110
    end if

    if ( lxt <= 1 ) then
      mflag = -1
      left = 1
      return
    end if

    ilo = lxt - 1
    ihi = lxt

  end if

  if ( xt(ihi) <= x ) then
    go to 20
  end if

  if ( xt(ilo) <= x ) then
    mflag = 0
    left = ilo
    return
  end if
!
!  Now X < XT(ILO).  Decrease ILO to capture X.
!
  istep = 1

10 continue

  ihi = ilo
  ilo = ihi - istep

  if ( 1 < ilo ) then
    if ( xt(ilo) <= x ) then
      go to 50
    end if
    istep = istep * 2
    go to 10
  end if

  ilo = 1

  if ( x < xt(1) ) then
    mflag = -1
    left = 1
    return
  end if

  go to 50
!
!  Now XT(IHI) <= X.  Increase IHI to capture X.
!
20 continue

  istep = 1

30 continue

  ilo = ihi
  ihi = ilo + istep

  if ( ihi < lxt ) then

    if ( x < xt(ihi) ) then
      go to 50
    end if

    istep = istep * 2
    go to 30

  end if

  if ( xt(lxt) <= x ) then
    go to 110
  end if
!
!  Now XT(ILO) < = X < XT(IHI).  Narrow the interval.
!
  ihi = lxt

50 continue

  do

    middle = ( ilo + ihi ) / 2

    if ( middle == ilo ) then
      mflag = 0
      left = ilo
      return
    end if
!
!  It is assumed that MIDDLE = ILO in case IHI = ILO+1.
!
    if ( xt(middle) <= x ) then
      ilo = middle
    else
      ihi = middle
    end if

  end do
!
!  Set output and return.
!
110 continue

  mflag = 1

  if ( x == xt(lxt) ) then
    mflag = 0
  end if

  do left = lxt, 1, -1
    if ( xt(left) < xt(lxt) ) then
      return
    end if
  end do

  return
end subroutine interv


    !---------------------------------------------------------------
    subroutine get_info(self, ao_geoinfo)
        implicit none
        type(GEOMETRY), intent(in) :: self
        type(GEOMETRY_INFO), intent(inout) :: ao_geoinfo
        ! LOCAL
        integer :: li_d
        integer :: li_nen

        li_nen = 1
        do li_d = 1, self % oi_dim
            li_nen = li_nen * (self % opi_P(li_d) + 1)
        end do
        ao_geoinfo % oi_nen = li_nen
    end subroutine get_info
    !---------------------------------------------------------------
    subroutine get_nen_geometry(self, ai_elt, ai_nen)
        !> \todo this must be generalized to geometries with diff nen on elts
        !> \todo we must add a specific array for this
        implicit none
        type(GEOMETRY), intent(in) :: self
        INTEGER :: ai_elt
        INTEGER, intent(inout) :: ai_nen
        ! LOCAL
        integer :: li_d
        integer :: li_nen

        li_nen = 1
        do li_d = 1, self % oi_dim
            li_nen = li_nen * (self % opi_P(li_d) + 1)
        end do

        ai_nen = li_nen        

    end subroutine get_nen_geometry
    !---------------------------------------------------------------
    subroutine assembly_deriv_basis_vect(self, ai_nderiv, api_npts, apr_works, apr_dBasis, api_leftmk)
        !> apr_dB contains values of all basis and their derivatives at apr_work
        !> for each direction li_d, all values ar in apr_B(li_d,1:nderiv,1:liP+1, 1:li_dimx)
        !> it is supposed that all works are in the same IR^d element
        implicit none
        type(GEOMETRY), intent(in) :: self
        integer, intent(in) :: ai_nderiv
        integer, dimension(:), intent(in) :: api_npts
        real(wp), dimension(:,:), intent(in) :: apr_works
        real(wp), dimension(0:,0:,:), intent(out) :: apr_dBasis
        integer, dimension(:), intent(out) :: api_leftmk
        ! LOCAL
        integer :: li_d
        integer :: li_index
        integer :: li_left
        integer :: li_leftmk
        integer :: li_maxP
        integer :: li_mflag
        integer :: li_dimx
        integer :: li_i
        integer :: li_p
        integer :: li_der
        real(wp) :: lr_w
        real(wp) :: lr_dw
        real(wp) :: lr_d2w
        real(wp), dimension(:,:), pointer :: lpr_a
        real(wp), dimension(:,:), pointer :: lpr_deltal
        real(wp), dimension(:,:), pointer :: lpr_deltar
        real(wp), dimension(:), pointer :: lpr_saved
        real(wp), dimension(:), pointer :: lpr_term
!        real(wp), dimension(:,:,:), pointer :: lpr_dB
        real(wp), dimension(:,:,:), pointer :: lpr_dB
        integer, dimension(self%oi_dim) :: lpi_left
        integer, dimension(self%oi_dim) :: lpi_i
        integer :: nx, ny, px, py, rx, ry
        integer :: N, d
        integer :: i, j, k 
        real   (kind=8), dimension(:), pointer :: Ux, Uy, X, Y 
        real   (kind=8), dimension(:,:), pointer :: Ww
        real   (kind=8), dimension(:,:,:), pointer :: Basis 

#ifdef _TRACE
        CALL printlog("assembly_deriv_basis_vect : Start", ai_dtllevel = mi_dtllevel_base + 1)
#endif

        do li_d = 1, self % oi_dim
            call interv(self % opr_u(li_d,:) &
            , self % opi_N(li_d) + 1 &
            , apr_works(li_d,1) &
            , lpi_left(li_d) &
            , li_mflag)

            api_leftmk(li_d) = lpi_left(li_d) - self % opi_P(li_d) - 1
        end do

        ! ****************************************************************
        !                           1D CASE
        ! ****************************************************************
        if (self % oi_dim==1) then
            N = ai_nderiv

            nx = self % opi_N(1) - 1
            px = self % opi_P(1)
            rx = api_npts(1) - 1

            call AssembleBasis1(ai_nderiv,N, self % oi_rational &
            , nx, px, self % opr_u(1,1:nx+px+2) &
            , self % opr_Ww(0, 0, 0:nx) &
            , rx, apr_works(1,1:rx+1) &
            , apr_dBasis (0:N,0:(px+1)-1,1:(rx+1)))
        end if
        ! ****************************************************************

        ! ****************************************************************
        !                           2D CASE
        ! ****************************************************************
        if (self % oi_dim==2) then
            if (ai_nderiv==1) then
                N = 3-1
            else if(ai_nderiv==2) then
                N = 6-1 
            end if
            nx = self % opi_N(1) - 1
            ny = self % opi_N(2) - 1
            px = self % opi_P(1)
            py = self % opi_P(2)
            rx = api_npts(1) - 1
            ry = api_npts(2) - 1

            call AssembleBasis2(ai_nderiv,N, self % oi_rational &
            , nx, px, self % opr_u(1,1:nx+px+2) &
            , ny, py, self % opr_u(2,1:ny+py+2) &
            , self % opr_Ww(0, 0:ny, 0:nx) &
            , rx, apr_works(1,1:rx+1) &
            , ry, apr_works(2,1:ry+1) &
            , apr_dBasis (0:N,0:(px+1)*(py+1)-1,1:(rx+1)*(ry+1)))
        end if
        ! ****************************************************************

#ifdef _TRACE
        CALL printlog("assembly_deriv_basis_vect : End", ai_dtllevel = mi_dtllevel_base + 1)
#endif
    end subroutine assembly_deriv_basis_vect
    !---------------------------------------------------------------
    subroutine assembly_deriv_point_vect_new(self, ai_elt, ai_npts, apr_dBasis, api_leftmk, apr_P)
        implicit none
        type(GEOMETRY), intent(in) :: self
        integer, intent(in) :: ai_elt
        integer, intent(in) :: ai_npts
        real(wp), dimension(0:,0:,:), intent(in) :: apr_dBasis
        integer, dimension(:), intent(in) :: api_leftmk
        real(wp), dimension(0:,:,:), intent(out) :: apr_P
        ! LOCAL
        integer :: li_index
        integer :: li_b
        integer :: li_d
        integer :: li_deriv
        integer :: li_dim
        integer :: li_nderiv_code
        integer, dimension(self % oi_dim) :: lpi_i
        integer, dimension(self % oi_dim) :: lpi_bi
        integer, dimension(self % oi_dim) :: lpi_Pp1
        type(GEOMETRY_INFO) :: lo_geoinfo
#ifdef _TRACE
        CALL printlog("assembly_deriv_point_vect_new : Start", ai_dtllevel = mi_dtllevel_base + 1)
#endif

        CALL get_info(self, lo_geoinfo)

        li_nderiv_code = SIZE(apr_P, 1) - 1
!        print *, 'li_nderiv_code = ', li_nderiv_code

        li_dim = self % oi_dim
        lpi_Pp1 = self % opi_P + 1
        apr_P = 0.0_wp
        do li_b = 0, lo_geoinfo % oi_nen - 1

            lpi_bi (1:li_dim) = self % opi_bi (ai_elt,li_b+1,1:li_dim)

            lpi_i = api_leftmk + lpi_bi
            li_index = get_index(self % opi_N, lpi_i)

            !**************************************************************************************
            !	COMPUTING THE POINT
            !**************************************************************************************
!            if (self % oi_dim==1) then
!                print *, '=================='
!            endif

            do li_deriv = 0, li_nderiv_code
!            if (self % oi_dim==1) then
!                print *, 'li_deriv ', li_deriv
!                    print *, apr_dBasis (li_deriv, li_b,1:ai_npts) 
!            endif

                ! we take the minimum in the case of the geometry is stored in 3D
                do li_d = 1, self % oi_Rd
                    apr_P(li_deriv,li_d,1:ai_npts) = apr_P(li_deriv,li_d,1:ai_npts)     &
                    + apr_dBasis (li_deriv, li_b,1:ai_npts) * self % opr_P(li_d,li_index)
                end do
            end do
            !**************************************************************************************

        end do

#ifdef _TRACE
        CALL printlog("assembly_deriv_point_vect_new : End", ai_dtllevel = mi_dtllevel_base + 1)
#endif
    end subroutine assembly_deriv_point_vect_new    
    !---------------------------------------------------------------
    subroutine compute_leftmk(self, apr_works, api_leftmk)
        implicit none
        type(GEOMETRY), intent(in) :: self
        real(wp), dimension(:), intent(in) :: apr_works
        integer, dimension(:), intent(out) :: api_leftmk
        ! LOCAL
        integer :: li_d
        integer :: li_mflag
        integer, dimension(self%oi_dim) :: lpi_left

#ifdef _TRACE
        CALL printlog("compute_leftmk : Start", ai_dtllevel = mi_dtllevel_base + 1)
#endif

        do li_d = 1, self % oi_dim
            call interv(self % opr_u(li_d,:) &
            , self % opi_N(li_d) + 1 &
            , apr_works(li_d) &
            , lpi_left(li_d) &
            , li_mflag)

            api_leftmk(li_d) = lpi_left(li_d) - self % opi_P(li_d) - 1
        end do

#ifdef _TRACE
        CALL printlog("compute_leftmk : End", ai_dtllevel = mi_dtllevel_base + 1)
#endif
    end subroutine compute_leftmk
    !---------------------------------------------------------------
    integer function get_index(api_N, api_i)
        !> this routine gives the (i,j) weight
        implicit none
        integer, dimension (:), intent(in) :: api_N
        integer, dimension (:), intent(in) :: api_i
        ! LOCAL
        integer :: li_index
        integer :: li_dim
        integer :: li_d
        integer :: li_j
        integer :: li_prod
#ifdef _TRACE
        CALL printlog("get_index : Start", ai_dtllevel = mi_dtllevel_base + 5)
#endif
        li_dim = SIZE(api_N)

        li_prod = 1
        do li_j = 1, li_dim
            li_prod = li_prod * api_N(li_j)
        end do

        li_index = api_i(1)
        do li_d = li_dim, 2, -1
            li_prod = li_prod / api_N(li_d)
            li_index = li_index + api_i(li_d) * li_prod
        end do

        get_index = li_index
#ifdef _TRACE
        CALL printlog("get_index : End", ai_dtllevel = mi_dtllevel_base + 5)
#endif
    end function get_index
    !---------------------------------------------------------------
    subroutine get_indices(api_N, ai_index, api_i)
        !> this routine gives the (i,j) weight
        !> api_i must be of dimension self%oi_dim
        implicit none
        integer, dimension (:), intent(in) :: api_N
        integer, intent(in) :: ai_index
        integer, dimension (:), intent(out) :: api_i
        ! LOCAL
        integer :: li_dim
        integer :: li_d
        integer :: li_j
        integer :: li_r
        integer :: li_current
        integer :: li_prod
#ifdef _TRACE
        CALL printlog("get_indices : Start", ai_dtllevel = mi_dtllevel_base + 5)
#endif
        li_dim = SIZE(api_N)

        li_prod = 1
        do li_j = 1, li_dim
            li_prod = li_prod * api_N(li_j)
        end do

        li_current = ai_index
        do li_d = li_dim, 1, -1

            li_prod = li_prod / api_N(li_d)

            api_i(li_d) = int(li_current/ li_prod)
            li_current = mod(li_current, li_prod)

        end do
#ifdef _TRACE
        CALL printlog("get_indices : End", ai_dtllevel = mi_dtllevel_base + 5)
#endif
    end subroutine get_indices    
    !---------------------------------------------------------------
    subroutine set_geometry_weights( self, ai_direction, apr_values, ai_dim)
        implicit none
        type(GEOMETRY) :: self
        INTEGER :: ai_direction
        INTEGER :: ai_dim
        real(wp), dimension(ai_dim) :: apr_values
        ! LOCAL

        CALL printlog("set_geometry_weights : Start", ai_dtllevel = mi_dtllevel_base + 1)

        self % opr_L (ai_direction, 1:ai_dim) = apr_values (1:ai_dim)

        CALL printlog("set_geometry_weights : End", ai_dtllevel = mi_dtllevel_base + 1)

    end subroutine set_geometry_weights
        !---------------------------------------------------------------
        subroutine geo_from_code_to_nderiv(ai_dim, ai_nderiv_code, api_deriv)
            !> compute the corresponding deriv indices
            !> to put in a new file grids_tools
            implicit none
            integer, intent(in) :: ai_dim
            integer, intent(in) :: ai_nderiv_code
            integer, dimension (:), intent(out) :: api_deriv
            ! LOCAL

            if (ai_dim == 1) then
                api_deriv ( 1 ) = ai_nderiv_code
                return
            end if

            if (ai_dim == 2) then
                select case ( ai_nderiv_code)
                    case ( 0 )
                        api_deriv ( 1 ) = 0
                        api_deriv ( 2 ) = 0
                        return
                    case ( 1 )
                        api_deriv ( 1 ) = 1
                        api_deriv ( 2 ) = 0
                        return
                    case ( 2 )
                        api_deriv ( 1 ) = 0
                        api_deriv ( 2 ) = 1
                        return
                    case ( 3 )
                        api_deriv ( 1 ) = 2
                        api_deriv ( 2 ) = 0
                        return
                    case ( 4 )
                        api_deriv ( 1 ) = 1
                        api_deriv ( 2 ) = 1
                        return
                    case ( 5 )
                        api_deriv ( 1 ) = 0
                        api_deriv ( 2 ) = 2
                        return

                    case Default
                        print*,"get_deriv_code: Type code Not Yet implemented"
                        return
                end select
            end if

            if (ai_dim == 3) then
                select case ( ai_nderiv_code)
                    case ( 0 )
                        api_deriv ( 1 ) = 0
                        api_deriv ( 2 ) = 0
                        api_deriv ( 3 ) = 0
                        return
                    case ( 1 )
                        api_deriv ( 1 ) = 1
                        api_deriv ( 2 ) = 0
                        api_deriv ( 3 ) = 0
                        return
                    case ( 2 )
                        api_deriv ( 1 ) = 0
                        api_deriv ( 2 ) = 1
                        api_deriv ( 3 ) = 0
                        return
                    case ( 3 )
                        api_deriv ( 1 ) = 0
                        api_deriv ( 2 ) = 0
                        api_deriv ( 3 ) = 1
                        return
                    case ( 4 )
                        api_deriv ( 1 ) = 2
                        api_deriv ( 2 ) = 0
                        api_deriv ( 3 ) = 0
                        return
                    case ( 5 )
                        api_deriv ( 1 ) = 1
                        api_deriv ( 2 ) = 1
                        api_deriv ( 3 ) = 0
                        return
                    case ( 6 )
                        api_deriv ( 1 ) = 1
                        api_deriv ( 2 ) = 0
                        api_deriv ( 3 ) = 1
                        return
                    case ( 7 )
                        api_deriv ( 1 ) = 0
                        api_deriv ( 2 ) = 2
                        api_deriv ( 3 ) = 0
                        return
                    case ( 8 )
                        api_deriv ( 1 ) = 0
                        api_deriv ( 2 ) = 1
                        api_deriv ( 3 ) = 1
                        return
                    case ( 9 )
                        api_deriv ( 1 ) = 0
                        api_deriv ( 2 ) = 0
                        api_deriv ( 3 ) = 2
                        return

                    case Default
                        print*,"get_deriv_code: Type code Not Yet implemented"
                        return
                end select
            end if

        end subroutine geo_from_code_to_nderiv
end module geometry_tools
