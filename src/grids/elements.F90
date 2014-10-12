!     
! File:   elements.F90
! Author: root
!
! Created on January 2, 2012, 5:10 PM
!

module elements_module
    use grids_def
    use tracelog_module
    use geometries_def
    use geometry_module
    use geometries_module
    implicit none

contains
    !----------------------------------------------------------------------------------------------
    subroutine create_elts_elements(self, api_npts)
        implicit none
        type(GRID_DATA) :: self
        !> the number of points per element/cell
        INTEGER, dimension(:) :: api_npts
        ! LOCAL VARIABLES
        integer :: li_elt
        integer :: li_npts

        CALL printlog("create_elts_elements : Start", ai_dtllevel = 1)

        ALLOCATE ( self % opo_elts ( self % oi_nel ) )

        print*, 'api_npts=', api_npts
        do li_elt = 1, self % oi_nel

            li_npts = api_npts(li_elt)
            
            ALLOCATE ( self % opo_elts (li_elt) % opr_pts(li_npts, self % oi_dim))
            ALLOCATE ( self % opo_elts (li_elt) % opr_w(li_npts))
            ! not done yet
            !ALLOCATE ( self % opo_elts ( li_elt )%opl_activ_pts ( li_npts) )

            ! for the moment, there is no local refinement, so all elements are activated
            self % opo_elts (li_elt) % ol_activ = .TRUE.

            self % opo_elts (li_elt) % oi_npts = li_npts

!            ALLOCATE ( self % opo_elts (li_elt) % opr_values(0:self % oi_nfields-1, self % oi_dof, li_npts))
        end do

        CALL printlog("create_elts_elements : End", ai_dtllevel = 1)

    end subroutine create_elts_elements
    !----------------------------------------------------------------------------------------------
    subroutine free_elements_elts(self)
        implicit none
        type(GRID_DATA) :: self
        ! LOCAL VARIABLES
        integer :: li_elt

        do li_elt = 1, self % oi_nel

            DEALLOCATE ( self % opo_elts (li_elt) % opr_pts)
            DEALLOCATE ( self % opo_elts (li_elt) % opr_w)
!            DEALLOCATE ( self % opo_elts (li_elt) % opr_values)

        end do

        DEALLOCATE ( self % opo_elts )

    end subroutine free_elements_elts
    !----------------------------------------------------------------------------------------------
    subroutine create_elts_from_tensor_product_1D(self, apr_localx, apr_localw)
        implicit none
        type(GRID_DATA) :: self
        !> apr_localx contains for each element a local grid
        !> access : apr_localx ( li_elt, li_pt )
        !> in apr_localx ( li_elt, 0 ) we must put the number of pts inside the current element
        real(wp), dimension(:,0:) :: apr_localx
        !> local weight
        real(wp), dimension(:,:) :: apr_localw
        ! LOCAL
        integer :: li_elt
        integer :: li_i
        integer :: li_npts
        integer :: li_index

        li_index = 0
        DO li_elt = 1, self % oi_nel

            ! getting the number of points
            li_npts = int(apr_localx(li_elt,0))

            DO li_i = 1, li_npts

                li_index = li_index + 1

                self % opo_elts (li_elt) % opr_pts(li_index, 1) = apr_localx(li_elt,li_i)
                
                self % opo_elts (li_elt) % opr_w(li_index) = apr_localw(li_elt,li_i)

            END DO

            if ( self % opo_elts (li_elt) % oi_npts /= li_index ) then
                print *, 'Error create_elts_from_tensor_product_1D: given a wrong number of points'
            end if
            li_index = 0
            
        END DO
        
    end subroutine create_elts_from_tensor_product_1D
    !----------------------------------------------------------------------------------------------
    subroutine create_elts_from_tensor_product_2D(self, apr_localx, apr_localy, apr_localwx, apr_localwy)
        implicit none
        type(GRID_DATA) :: self
        !> apr_localx contains for each element a local grid
        !> access : apr_localx ( li_elt, li_pt )
        !> in apr_localx ( li_elt, 0 ) we must put the number of pts inside the current element
        real(wp), dimension(:,0:) :: apr_localx
        real(wp), dimension(:,0:) :: apr_localy
        real(wp), dimension(:,:) :: apr_localwx
        real(wp), dimension(:,:) :: apr_localwy
        ! LOCAL
        integer :: li_elt
        integer :: li_ex
        integer :: li_ey
        integer :: li_nex
        integer :: li_ney
        integer :: li_i
        integer :: li_j
        integer :: li_nx
        integer :: li_ny
        integer :: li_index

        li_nex = SIZE(apr_localx,1)
        li_ney = SIZE(apr_localy,1)

        li_index = 0
        li_elt = 0
        DO li_ey = 1, li_ney

            DO li_ex = 1, li_nex

                li_elt = li_elt + 1

                ! getting the number of points
                li_nx = int(apr_localx(li_ex,0))
                li_ny = int(apr_localy(li_ey,0))
                
                DO li_j = 1, li_ny

                    DO li_i = 1, li_nx

                        li_index = li_index + 1

                        self % opo_elts (li_elt) % opr_pts(li_index, 1) = apr_localx(li_ex,li_i)
                        self % opo_elts (li_elt) % opr_pts(li_index, 2) = apr_localy(li_ey,li_j)

                        self % opo_elts (li_elt) % opr_w(li_index) = apr_localwx(li_ex,li_i) * apr_localwy(li_ey,li_j)

                    END DO

                END DO

                if ( self % opo_elts (li_elt) % oi_npts /= li_index ) then
                    print *, 'Error create_elts_from_tensor_product_2D: given a wrong number of points'
                end if
                li_index = 0

            END DO

        END DO
        
    end subroutine create_elts_from_tensor_product_2D
    !----------------------------------------------------------------------------------------------
    subroutine create_elts_from_tensor_product_3D(self, apr_localx, apr_localy, apr_localz, apr_localwx, apr_localwy, apr_localwz)
        implicit none
        type(GRID_DATA) :: self
        !> apr_localx contains for each element a local grid
        !> access : apr_localx ( li_elt, li_pt )
        !> in apr_localx ( li_elt, 0 ) we must put the number of pts inside the current element
        real(wp), dimension(:,0:) :: apr_localx
        real(wp), dimension(:,0:) :: apr_localy
        real(wp), dimension(:,0:) :: apr_localz
        real(wp), dimension(:,:) :: apr_localwx
        real(wp), dimension(:,:) :: apr_localwy
        real(wp), dimension(:,:) :: apr_localwz
        ! LOCAL
        integer :: li_elt
        integer :: li_ex
        integer :: li_ey
        integer :: li_ez
        integer :: li_nex
        integer :: li_ney
        integer :: li_nez
        integer :: li_i
        integer :: li_j
        integer :: li_k
        integer :: li_nx
        integer :: li_ny
        integer :: li_nz
        integer :: li_index

        li_nex = SIZE(apr_localx,1)
        li_ney = SIZE(apr_localy,1)
        li_nez = SIZE(apr_localz,1)

        li_index = 0
        li_elt = 0
        DO li_ez = 1, li_nez

            DO li_ey = 1, li_ney

                DO li_ex = 1, li_nex

                    li_elt = li_elt + 1

                    ! getting the number of points
                    li_nx = int(apr_localx(li_ex,0))
                    li_ny = int(apr_localy(li_ey,0))
                    li_nz = int(apr_localz(li_ez,0))

                    DO li_k = 1, li_nz

                        DO li_j = 1, li_ny

                            DO li_i = 1, li_nx

                                li_index = li_index + 1

                                self % opo_elts (li_elt) % opr_pts(li_index, 1) = apr_localx(li_ex,li_i)
                                self % opo_elts (li_elt) % opr_pts(li_index, 2) = apr_localy(li_ey,li_j)
                                self % opo_elts (li_elt) % opr_pts(li_index, 3) = apr_localz(li_ez,li_k)

                                self % opo_elts (li_elt) % opr_w(li_index) = apr_localwx(li_ex,li_i) &
                                * apr_localwy(li_ey,li_j) &
                                * apr_localwz(li_ez,li_k)

                            END DO

                        END DO

                    END DO

                    if ( self % opo_elts (li_elt) % oi_npts /= li_index ) then
                        print *, 'Error create_elts_from_tensor_product_3D: given a wrong number of points'
                    end if
                    li_index = 0
                    
                END DO

            END DO

        END DO
        
    end subroutine create_elts_from_tensor_product_3D
    !----------------------------------------------------------------------------------------------
    !> this routine gives the id of the element where belongs the points pt
    !> \todo to implement
    subroutine find_element_elts (self,apr_pt)
        implicit none
        type(GRID_DATA) :: self
        real(wp), dimension(:) :: apr_pt
        print *, 'find_element_elts : NOT YET IMPLEMENTED'
        return
    end subroutine find_element_elts
    !----------------------------------------------------------------------------------------------
    !> this routine gives the id of the element where belongs the points pt
    !> \todo to implement
    subroutine insert_point_elts (self,apr_pt)
        implicit none
        type(GRID_DATA) :: self
        real(wp), dimension(:) :: apr_pt
        print *, 'insert_point_elts : NOT YET IMPLEMENTED'
        return
    end subroutine insert_point_elts
    !----------------------------------------------------------------------------------------------
    subroutine print_grid_elts(self)
        implicit none
        type(GRID_DATA) :: self
        ! LOCAL VARIABLES
        integer :: li_elt
        integer :: li_npts
        integer :: li_d

        CALL printlog("print_grid_elts : Start", ai_dtllevel = 1)

        print *, '----------------'
        print *, 'nel = ', self % oi_nel
        print *, 'dim = ', self % oi_dim
        print *, 'ngrid = ', self % oi_ngrid
        print *, '----------------'

        do li_elt = 1, self % oi_nel

            li_npts = self % opo_elts (li_elt) % oi_npts
            print *, '**********************'
            print *, 'element = ', li_elt
            print *, '**********************'
            print *, 'npts = ', li_npts
            print *, 'points = '
            do li_d = 1, self % oi_dim
                print *, 'd = ', li_d
                print *, self % opo_elts (li_elt) % opr_pts(1:li_npts, li_d)
            end do
            print *, 'weights = ', self % opo_elts (li_elt) % opr_w(1:li_npts)

        end do

        CALL printlog("print_grid_elts : End", ai_dtllevel = 1)

    end subroutine print_grid_elts
    !---------------------------------------------------------------
    subroutine assembly_elts(self, ao_geos, ai_nderiv, ai_id, ai_elt, apr_points)
        implicit none
        type(GRID_DATA) :: self
        type(GEOMETRIES) :: ao_geos
        integer :: ai_nderiv
        integer :: ai_id !id of the geometry (in multi-patchs)
        integer :: ai_elt ! the id of the element
        real(wp), dimension(:,:,:), intent(out) :: apr_points !(ai_nderiv,ai_Rd,ai_npts)
        ! LOCAL
        real(wp), dimension(:,:), pointer :: lpr_sites
        integer :: li_Rd
        integer :: li_npts

        CALL printlog("assembly_elts : Start", ai_dtllevel = 1)

        li_npts = self % opo_elts (ai_elt) % oi_npts
        li_Rd = ao_geos % opo_geo(ai_id) % oi_Rd

        ALLOCATE(lpr_sites(self % oi_dim,li_npts))

        apr_points = 0.0

        print *, 'NOT YET IMPLEMENTED'

        DEALLOCATE(lpr_sites)

        CALL printlog("assembly_elts : End", ai_dtllevel = 1)

    end subroutine assembly_elts
end module elements_module


