!     
! File:   elements_tensor.F90
! Author: root
!
! Created on January 4, 2012, 10:56 PM
!

module elements_tensor_module
    use grids_def
    use tracelog_module
    use geometries_def
    use geometry_module
    use geometries_module
    use bbox_def
    implicit none

#ifdef _DEBUG
    INTEGER, PARAMETER, PRIVATE  :: mi_dtllevel_base = 0
#else
    INTEGER, PARAMETER, PRIVATE  :: mi_dtllevel_base = 2
#endif

contains
    !----------------------------------------------------------------------------------------------
    subroutine create_elts_t_elements(self, api_npts, ai_totalnpts &
        , ai_maxnparams_matrices, ai_maxnparams_fields, ai_maxnparams_norms)
        implicit none
        TYPE(GRID_DATA) :: self
        !> the number of points per element/cell
        INTEGER, DIMENSION(:) :: api_npts
        INTEGER :: ai_totalnpts
        INTEGER, optional :: ai_maxnparams_matrices
        INTEGER, optional :: ai_maxnparams_fields
        INTEGER, optional :: ai_maxnparams_norms
        ! LOCAL VARIABLES
        INTEGER :: li_dim
        INTEGER :: li_elt
        INTEGER :: li_npts
#ifdef _TRACE
        CALL printlog("create_elts_t_elements : Start", ai_dtllevel = mi_dtllevel_base + 1)
#endif

        li_dim = self % oi_dim
        ALLOCATE ( self % opo_elts (self % oi_nel) )

        do li_elt = 1, self % oi_nel

            self % opo_elts (li_elt) % ol_stored = .FALSE.
            
            li_npts = api_npts(li_elt)
!print*,'xxxx'
            ALLOCATE ( self % opo_elts (li_elt) % opr_pts(li_npts, li_dim))
!print*,'xxxx'
            ALLOCATE ( self % opo_elts (li_elt) % opr_w(ai_totalnpts))
!print*,'xxxx'
            ALLOCATE ( self % opo_elts (li_elt) % opi_npts(li_dim))
!print*,'xxxx'
            if ( self % ol_use_values_f ) then
                ALLOCATE ( self % opo_elts (li_elt) % opr_values_f(0:self % oi_nfields-1, ai_maxnparams_fields, ai_totalnpts))
            end if
!print*,'xxxx'
            if ( self % ol_use_values_m ) then
                ALLOCATE ( self % opo_elts (li_elt) % opr_values_m(0:self % oi_nmatrices-1, ai_maxnparams_matrices, ai_totalnpts))
            end if

!print*,'xxxx'
            if ( self % ol_use_values_n ) then
                ALLOCATE ( self % opo_elts (li_elt) % opr_values_n(0:self % oi_nnorms-1, ai_maxnparams_norms, ai_totalnpts))
            end if
!print*,'xxxx'

            self % opo_elts (li_elt) % opi_npts = li_npts
            
            ! not done yet
            !ALLOCATE ( self ( li_elt )%opl_activ_pts ( li_npts) )

            ! for the moment, there is no local refinement, so all elements are activated
            self % opo_elts (li_elt) % ol_activ = .TRUE.
!print*,'xxxx'
            self % opo_elts (li_elt) % oi_npts = ai_totalnpts
!            print *, 'li_npts, ai_totalnpts=', li_npts, ai_totalnpts
        end do
        !> used only for tensor
        self % oi_maxdir_npts = MAXVAL(api_npts)        

#ifdef _TRACE
        CALL printlog("create_elts_t_elements : End", ai_dtllevel = mi_dtllevel_base + 1)
#endif

    end subroutine create_elts_t_elements
    !----------------------------------------------------------------------------------------------
    subroutine free_elements_elts_t(self)
        implicit none
        TYPE(GRID_DATA) :: self
        ! LOCAL VARIABLES
        INTEGER :: li_elt
#ifdef _TRACE
        CALL printlog("free_elements_elts_t : Start", ai_dtllevel = mi_dtllevel_base + 1)
#endif

        do li_elt = 1, self % oi_nel

            DEALLOCATE ( self % opo_elts (li_elt) % opr_pts)
            DEALLOCATE ( self % opo_elts (li_elt) % opr_w)
            DEALLOCATE ( self % opo_elts (li_elt) % opi_npts)
            if ( self % ol_use_values_f ) then
                DEALLOCATE ( self % opo_elts (li_elt) % opr_values_f)
            end if
            if ( self % ol_use_values_m ) then
                DEALLOCATE ( self % opo_elts (li_elt) % opr_values_m)
            end if
            if ( self % ol_use_values_n ) then
                DEALLOCATE ( self % opo_elts (li_elt) % opr_values_n)
            end if

        end do

        DEALLOCATE ( self % opo_elts )

        IF ( ALLOCATED(self % opr_dBasisatx)) THEN
            DEALLOCATE( self % opr_dBasisatx )
        END IF
#ifdef _TRACE
        CALL printlog("free_elements_elts_t : End", ai_dtllevel = mi_dtllevel_base + 1)
#endif

    end subroutine free_elements_elts_t
    !----------------------------------------------------------------------------------------------
    subroutine create_elts_t_from_tensor_product_1D(self, apr_localx, apr_localw)
        implicit none
        TYPE(GRID_DATA) :: self
        !> apr_localx contains for each element a local grid
        !> access : apr_localx ( li_elt, li_pt )
        !> in apr_localx ( li_elt, 0 ) we must put the number of pts inside the current element
        REAL(wp), DIMENSION(:,0:) :: apr_localx
        !> local weight
        REAL(wp), DIMENSION(:,:) :: apr_localw
        ! LOCAL
        INTEGER :: li_elt
        INTEGER :: li_i
        INTEGER :: li_npts
        INTEGER :: li_index
#ifdef _TRACE
        CALL printlog("create_elts_t_from_tensor_product_1D : Start", ai_dtllevel = mi_dtllevel_base + 1)
#endif

        li_index = 0
        DO li_elt = 1, self % oi_nel

            ! getting the number of points
            li_npts = int(apr_localx(li_elt,0))
!            print*,'li_elt,li_npts=',li_elt,li_npts
            self % opo_elts (li_elt) % opi_npts (1) = li_npts

            DO li_i = 1, li_npts

                li_index = li_index + 1

                self % opo_elts (li_elt) % opr_pts(li_index, 1) = apr_localx(li_elt,li_i)

                self % opo_elts (li_elt) % opr_w(li_index) = apr_localw(li_elt,li_i)

            END DO
!        print*,'opr_pts =', self % opo_elts (li_elt) % opr_pts(:,:)

            if ( self % opo_elts (li_elt) % oi_npts /= li_index ) then
                print *, 'Error create_elts_t_from_tensor_product_1D: given a wrong number of points'
            end if
            li_index = 0

        END DO
#ifdef _TRACE
        CALL printlog("create_elts_t_from_tensor_product_1D : End", ai_dtllevel = mi_dtllevel_base + 1)
#endif

    end subroutine create_elts_t_from_tensor_product_1D
    !----------------------------------------------------------------------------------------------
    subroutine create_elts_t_from_tensor_product_2D(self, apr_localx, apr_localy, apr_localwx, apr_localwy)
        implicit none
        TYPE(GRID_DATA) :: self
        !> apr_localx contains for each element a local grid
        !> access : apr_localx ( li_elt, li_pt )
        !> in apr_localx ( li_elt, 0 ) we must put the number of pts inside the current element
        REAL(WP), DIMENSION(:,0:) :: apr_localx
        REAL(WP), DIMENSION(:,0:) :: apr_localy
        REAL(WP), DIMENSION(:,:) :: apr_localwx
        REAL(WP), DIMENSION(:,:) :: apr_localwy
        ! LOCAL
        INTEGER :: li_elt
        INTEGER :: li_ex
        INTEGER :: li_ey
        INTEGER :: li_nex
        INTEGER :: li_ney
        INTEGER :: li_i
        INTEGER :: li_j
        INTEGER :: li_nx
        INTEGER :: li_ny
        INTEGER :: li_index
#ifdef _TRACE
        CALL printlog("create_elts_t_from_tensor_product_2D : Start", ai_dtllevel = mi_dtllevel_base + 1)
#endif

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

                self % opo_elts (li_elt) % opi_npts (1) = li_nx
                self % opo_elts (li_elt) % opi_npts (2) = li_ny

                self % opo_elts (li_elt) % oi_npts = li_nx * li_ny

                self % opo_elts (li_elt) % opr_pts(1:li_nx, 1) = apr_localx(li_ex,1:li_nx)
                self % opo_elts (li_elt) % opr_pts(1:li_ny, 2) = apr_localy(li_ey,1:li_ny)
                
                DO li_j = 1, li_ny

                    DO li_i = 1, li_nx

                        li_index = li_index + 1

                        self % opo_elts (li_elt) % opr_w(li_index) = apr_localwx(li_ex,li_i) * apr_localwy(li_ey,li_j)

                    END DO

                END DO

                if ( self % opo_elts (li_elt) % oi_npts /= li_index ) then
                    print *, 'Error create_elts_t_from_tensor_product_2D: given a wrong number of points'
                end if
                li_index = 0

            END DO

        END DO
#ifdef _TRACE
        CALL printlog("create_elts_t_from_tensor_product_2D : End", ai_dtllevel = mi_dtllevel_base + 1)
#endif

    end subroutine create_elts_t_from_tensor_product_2D
    !----------------------------------------------------------------------------------------------
    subroutine create_elts_t_from_tensor_product_3D(self, apr_localx, apr_localy, apr_localz, apr_localwx, apr_localwy, apr_localwz)
        implicit none
        TYPE(GRID_DATA) :: self
        !> apr_localx contains for each element a local grid
        !> access : apr_localx ( li_elt, li_pt )
        !> in apr_localx ( li_elt, 0 ) we must put the number of pts inside the current element
        REAL(WP), DIMENSION(:,0:) :: apr_localx
        REAL(WP), DIMENSION(:,0:) :: apr_localy
        REAL(WP), DIMENSION(:,0:) :: apr_localz
        REAL(WP), DIMENSION(:,:) :: apr_localwx
        REAL(WP), DIMENSION(:,:) :: apr_localwy
        REAL(WP), DIMENSION(:,:) :: apr_localwz
        ! LOCAL
        INTEGER :: li_elt
        INTEGER :: li_ex
        INTEGER :: li_ey
        INTEGER :: li_ez
        INTEGER :: li_nex
        INTEGER :: li_ney
        INTEGER :: li_nez
        INTEGER :: li_i
        INTEGER :: li_j
        INTEGER :: li_k
        INTEGER :: li_nx
        INTEGER :: li_ny
        INTEGER :: li_nz
        INTEGER :: li_index
#ifdef _TRACE
        CALL printlog("create_elts_t_from_tensor_product_3D : Start", ai_dtllevel = mi_dtllevel_base + 1)
#endif

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

                    self % opo_elts (li_elt) % opi_npts (1) = li_nx
                    self % opo_elts (li_elt) % opi_npts (2) = li_ny
                    self % opo_elts (li_elt) % opi_npts (3) = li_nz

                    self % opo_elts (li_elt) % oi_npts = li_nx * li_ny * li_nz

                    self % opo_elts (li_elt) % opr_pts(1:li_nx, 1) = apr_localx(li_ex,1:li_nx)
                    self % opo_elts (li_elt) % opr_pts(1:li_ny, 2) = apr_localy(li_ey,1:li_ny)
                    self % opo_elts (li_elt) % opr_pts(1:li_nz, 3) = apr_localz(li_ez,1:li_nz)


                    DO li_k = 1, li_nz

                        DO li_j = 1, li_ny

                            DO li_i = 1, li_nx

                                li_index = li_index + 1

                                self % opo_elts (li_elt) % opr_w(li_index) = apr_localwx(li_ex,li_i) &
                                * apr_localwy(li_ey,li_j) &
                                * apr_localwz(li_ez,li_k)

                            END DO

                        END DO

                    END DO

                    if ( self % opo_elts (li_elt) % oi_npts /= li_index ) then
                        print *, 'Error create_elts_t_from_tensor_product_3D: given a wrong number of points'
                    end if
                    li_index = 0

                END DO

            END DO

        END DO
#ifdef _TRACE
        CALL printlog("create_elts_t_from_tensor_product_3D : End", ai_dtllevel = mi_dtllevel_base + 1)
#endif

    end subroutine create_elts_t_from_tensor_product_3D
    !----------------------------------------------------------------------------------------------
    subroutine set_elts_t_unidirection_tensor_2d( self, ai_direction   &
        , apr_dirpts, apr_dirw  &
        , ar_fixpt, ar_fixw, ai_elt)
        implicit none
        TYPE(GRID_DATA) :: self
        !> apr_localx contains for each element a local grid
        !> access : apr_localx ( li_elt, li_pt )
        !> in apr_localx ( li_elt, 0 ) we must put the number of pts inside the current element
        INTEGER :: ai_direction
        real(8), DIMENSION(:,0:) :: apr_dirpts
        real(8), DIMENSION(:,:) :: apr_dirw
        real(8) :: ar_fixpt
        real(8) :: ar_fixw
        !> ai_elt is the id of the last created element
        INTEGER, intent(inout) :: ai_elt
        ! LOCAL 
        INTEGER :: li_elt
        INTEGER :: li_ne
        INTEGER :: li_e
        INTEGER :: li_npts
#ifdef _TRACE
        CALL printlog("set_elts_t_unidirection_tensor_2d : Start", ai_dtllevel = mi_dtllevel_base + 1)
#endif
        li_ne = self % oi_nel 
        DO li_e = 1, li_ne

            li_elt = ai_elt + li_e
            ! getting the number of points
            li_npts = int(apr_dirpts(li_e,0))

            self % opo_elts (li_elt) % opi_npts (ai_direction+1) = li_npts
            self % opo_elts (li_elt) % opi_npts (2-ai_direction) = 1

            self % opo_elts (li_elt) % opr_pts(1:li_npts, ai_direction+1) = apr_dirpts(li_e,1:li_npts)
            self % opo_elts (li_elt) % opr_pts(1, 2-ai_direction) = ar_fixpt

            self % opo_elts (li_elt) % opr_w(1:li_npts) = ar_fixw * apr_dirw(li_e,1:li_npts)

        END DO

        ai_elt = li_elt

#ifdef _TRACE
        CALL printlog("set_elts_t_unidirection_tensor_2d : End", ai_dtllevel = mi_dtllevel_base + 1)
#endif

    end subroutine set_elts_t_unidirection_tensor_2d
    !----------------------------------------------------------------------------------------------
    !> this routine gives the id of the element where belongs the points pt
    !> \todo to implement
    subroutine find_element (self,apr_pt)
        implicit none
        TYPE(GRID_DATA) :: self
        REAL(WP), DIMENSION(:) :: apr_pt
        print *, 'find_element : NOT YET IMPLEMENTED'
        return
    end subroutine find_element
    !----------------------------------------------------------------------------------------------
    !> this routine gives the id of the element where belongs the points pt
    !> \todo to implement
    subroutine insert_point_elts_t (self,apr_pt)
        implicit none
        TYPE(GRID_DATA) :: self
        REAL(WP), DIMENSION(:) :: apr_pt
        print *, 'insert_point_elts_t : NOT YET IMPLEMENTED'
        return
    end subroutine insert_point_elts_t
    !----------------------------------------------------------------------------------------------
    subroutine print_grid_elts_t(self)
        implicit none
        TYPE(GRID_DATA) :: self
        ! LOCAL VARIABLES
        INTEGER :: li_elt
        INTEGER :: li_npts
        INTEGER :: li_dim
        INTEGER :: li_d
#ifdef _TRACE
        CALL printlog("print_grid_elts_t : Start", ai_dtllevel = mi_dtllevel_base + 1)
#endif

        print *, '----------------'
        print *, 'nel = ', self % oi_nel
        print *, 'dim = ', self % oi_dim
        print *, 'ngrid = ', self % oi_ngrid
        print *, 'nfields = ', self % oi_nfields
        print *, '----------------'


        do li_elt = 1, self % oi_nel
            
            print *, '**********************'
            print *, 'element = ', li_elt
            print *, '**********************'            
            print *, 'points = '
            do li_d = 1, self % oi_dim
                li_npts = self % opo_elts (li_elt) % opi_npts (li_d)
                print *, 'd = ', li_d
                print *, 'npts = ', li_npts
                print *, self % opo_elts (li_elt) % opr_pts(1:li_npts, li_d)
            end do
            print *, 'weights = ', self % opo_elts (li_elt) % opr_w!(1:li_npts)
!            print *, 'fields values = '
!            do li_d = 0, self % oi_nfields-1
!                print *, self % opo_elts (li_elt) % opr_values_f(li_d, :, :)
!            end do
!            print *, 'matrices values = '
!            do li_d = 0, self % oi_nmatrices-1
!                print *, self % opo_elts (li_elt) % opr_values_m(li_d, :, :)
!            end do

        end do
#ifdef _TRACE
        CALL printlog("print_grid_elts_t : End", ai_dtllevel = mi_dtllevel_base + 1)
#endif

    end subroutine print_grid_elts_t
    !---------------------------------------------------------------
    subroutine assembly_points_elts_t(self, ao_geos, ai_id, ai_elt, ai_realelt, ai_npts, ao_bbox, apr_points)
        implicit none
        TYPE(GRID_DATA) :: self
        TYPE(GEOMETRIES) :: ao_geos
        INTEGER :: ai_id !id of the geometry (in multi-patchs)
        INTEGER :: ai_elt ! the id of the element
        INTEGER :: ai_realelt ! the id of the real element
        INTEGER :: ai_npts
        TYPE(BBOX), intent(in) :: ao_bbox
        REAL(WP), DIMENSION(0:,:,:), intent(out) :: apr_points !(ai_nderiv,ai_Rd,ai_npts)
        ! LOCAL
        INTEGER :: li_npts
#ifdef _TRACE
        CALL printlog("assembly_points_elts_t : Start", ai_dtllevel = mi_dtllevel_base + 1)
#endif

        li_npts = MAXVAL ( self % opo_elts (ai_elt) % opi_npts (:) )

        CALL assembly_deriv_point_vect_new(ao_geos % opo_geo(ai_id) &
        , ai_realelt    &
        , ai_npts   &
        , ao_bbox % opr_dBasis   &
        , ao_bbox % opi_leftmk  &
        , apr_points)
#ifdef _TRACE
        CALL printlog("assembly_points_elts_t : Start", ai_dtllevel = mi_dtllevel_base + 1)
#endif

    end subroutine assembly_points_elts_t
    !---------------------------------------------------------------
    subroutine assembly_basis_elts_t(self, ao_geos, ai_nderiv, ai_id, ai_elt, ai_realelt, ao_bbox, ai_nen)
        implicit none
        TYPE(GRID_DATA) :: self
        TYPE(GEOMETRIES) :: ao_geos
        INTEGER :: ai_nderiv
        INTEGER :: ai_id !id of the geometry (in multi-patchs)
        INTEGER :: ai_elt ! the id of the element
        INTEGER :: ai_realelt ! the id of the real element
        TYPE(BBOX), intent(inout) :: ao_bbox
        INTEGER, optional :: ai_nen
        ! LOCAL
        INTEGER :: li_dim
        INTEGER :: li_Rd
        INTEGER :: li_npts
        INTEGER :: li_dirnpts
        INTEGER :: li_b
        INTEGER :: li_pt
        INTEGER :: li_d
        INTEGER :: li_d_code
        INTEGER :: li_i
        INTEGER :: li_p
        REAL(WP), DIMENSION(:,:), pointer :: lpr_sites
        REAL(WP), DIMENSION(self % oi_dim) :: lpr_works
        INTEGER, DIMENSION(self % oi_dim) :: lpi_ni
        INTEGER, DIMENSION(self % oi_dim) :: lpi_ndiffum1
        INTEGER, DIMENSION(self % oi_dim) :: lpi_Pp1
        INTEGER, DIMENSION(self % oi_dim) :: lpi_bi
        INTEGER, DIMENSION(self % oi_dim) :: lpi_pti
        INTEGER, DIMENSION(self % oi_dim) :: lpi_npts
        INTEGER, DIMENSION(self % oi_dim) :: lpi_deriv
        TYPE(GEOMETRY), POINTER :: lp_geo

#ifdef _TRACE
        CALL printlog("assembly_basis_elts_t : Start", ai_dtllevel = mi_dtllevel_base + 1)
#endif
        lp_geo => ao_geos % opo_geo(ai_id)

        li_dim = self % oi_dim
        li_dirnpts = MAXVAL ( self % opo_elts (ai_elt) % opi_npts (:) )
        li_npts = self % opo_elts (ai_elt) % oi_npts
        li_Rd = lp_geo % oi_Rd

        ALLOCATE(lpr_sites(self%oi_dim,li_dirnpts))

        lpr_sites (1:li_dim,1:li_dirnpts) = TRANSPOSE (self % opo_elts (ai_elt) % opr_pts(1:li_dirnpts, 1:li_dim))

        CALL assembly_deriv_basis_vect(lp_geo, ai_nderiv &
        , self % opo_elts (ai_elt) % opi_npts (:), lpr_sites   &
        , ao_bbox % opr_dBasis, ao_bbox % opi_leftmk)

        DEALLOCATE(lpr_sites)
#ifdef _TRACE
        CALL printlog("assembly_basis_elts_t : End", ai_dtllevel = mi_dtllevel_base + 1)
#endif
    end subroutine assembly_basis_elts_t
    !----------------------------------------------------------------------------------------------
    subroutine create_stored_data_elts_t(self, api_nderiv, api_N, api_P)
        implicit none
        TYPE(GRID_DATA) :: self
        INTEGER, DIMENSION(:) :: api_nderiv
        INTEGER, DIMENSION(:) :: api_N
        INTEGER, DIMENSION(:) :: api_P
        ! LOCAL VARIABLES
        INTEGER :: li_d
#ifdef _TRACE
        CALL printlog("create_stored_data_elts_t : Start", ai_dtllevel = mi_dtllevel_base + 1)
#endif

        self % ol_usestored_data = .TRUE.
        ALLOCATE(self % opo_stored(self % oi_dim))

        do li_d = 1, self % oi_dim

            self % opo_stored(li_d) % oi_nderiv = api_nderiv(li_d)

            ALLOCATE ( self % opo_stored(li_d) % opr_dB (0:api_nderiv(li_d), api_N(li_d)    &
            , 0:api_P(li_d), self % oi_maxdir_npts) )

        end do
#ifdef _TRACE
        CALL printlog("create_stored_data_elts_t : End", ai_dtllevel = mi_dtllevel_base + 1)
#endif

    end subroutine create_stored_data_elts_t
    !----------------------------------------------------------------------------------------------
    subroutine free_stored_data_elts_t(self)
        implicit none
        TYPE(GRID_DATA) :: self
        ! LOCAL VARIABLES
        INTEGER :: li_d
#ifdef _TRACE
        CALL printlog("free_stored_data_elts_t : Start", ai_dtllevel = mi_dtllevel_base + 1)
#endif

        do li_d = 1, self % oi_dim

            DEALLOCATE ( self % opo_stored(li_d) % opr_dB )

        end do

        DEALLOCATE(self % opo_stored)
#ifdef _TRACE
        CALL printlog("free_stored_data_elts_t : End", ai_dtllevel = mi_dtllevel_base + 1)
#endif

    end subroutine free_stored_data_elts_t
    !----------------------------------------------------------------------------------------------
    subroutine get_weights_elts_t(self, apr_w)
        implicit none
        TYPE(GRID_DATA) :: self
        REAL(WP), DIMENSION(:,:), intent(inout) :: apr_w
        ! LOCAL VARIABLES
        INTEGER :: li_elt
#ifdef _TRACE
        CALL printlog("get_weights_elts_t : Start", ai_dtllevel = mi_dtllevel_base + 1)
#endif

        DO li_elt = 1, self % oi_nel
            apr_w(li_elt,:) = self % opo_elts (li_elt) % opr_w(:)
        END DO
#ifdef _TRACE
        CALL printlog("get_weights_elts_t : End", ai_dtllevel = mi_dtllevel_base + 1)
#endif

    end subroutine get_weights_elts_t
    !---------------------------------------------------------------
    subroutine set_elts_t_tensor_basis ( self, ai_d, apr_dbasisatx, ai_nderiv, ai_p, ai_ni, ai_dirnpts &
        , ai_maxnderiv, ai_maxp, ai_maxni, ai_maxdirnpts )
        implicit none
        TYPE(GRID_DATA) :: self
        INTEGER :: ai_d
        INTEGER :: ai_nderiv
        INTEGER :: ai_p
        INTEGER :: ai_ni
        INTEGER :: ai_dirnpts
        REAL(wp), DIMENSION(0:ai_nderiv, 0:ai_p, 1:ai_dirnpts, 1:ai_ni) :: apr_dbasisatx
        INTEGER :: ai_maxnderiv
        INTEGER :: ai_maxp
        INTEGER :: ai_maxni
        INTEGER :: ai_maxdirnpts
        ! LOCAL
#ifdef _TRACE
        CALL printlog("set_elts_t_tensor_basis : Start", ai_dtllevel = 1)
#endif

        IF (.NOT. ALLOCATED(self % opr_dBasisatx)) THEN
            ALLOCATE(self % opr_dBasisatx (self % oi_dim, 0:ai_maxnderiv, 0:ai_maxp, 1:ai_maxdirnpts, 1:ai_maxni))
        END IF

        self % opr_dBasisatx (ai_d, 0:ai_nderiv, 0:ai_p, 1:ai_dirnpts, 1:ai_ni) = &
        apr_dbasisatx (0:ai_nderiv, 0:ai_p, 1:ai_dirnpts, 1:ai_ni)
#ifdef _TRACE
        CALL printlog("set_elts_t_tensor_basis : End", ai_dtllevel = 1)
#endif

    end subroutine set_elts_t_tensor_basis
    !---------------------------------------------------------------
    subroutine get_elts_t_tensor_basis ( self, ai_d, ai_nderiv, ai_p, ai_ielt, ai_dirnpts, ao_bbox )
        implicit none
        TYPE(GRID_DATA) :: self
        INTEGER :: ai_d
        INTEGER :: ai_nderiv
        INTEGER :: ai_p
        INTEGER :: ai_ielt
        INTEGER :: ai_dirnpts
        TYPE(BBOX), intent(inout) :: ao_bbox
        ! LOCAL
        INTEGER :: li_pt, li_d

#ifdef _TRACE
        CALL printlog("get_elts_t_tensor_basis : Start", ai_dtllevel = 1)
        print *, 'get_elts_t_tensor_basis'
#endif

!#ifdef _DEBUG       
!        print *, "====== ielt ", ai_ielt, " ======"
!        print *, 'ai_d          = ', ai_d
!        print *, 'ai_nderiv     = ', ai_nderiv
!        print *, 'ai_p          = ', ai_p
!        print *, 'ai_dirnpts    = ', ai_dirnpts
!#endif

        ao_bbox % opr_dB(ai_d, 0:ai_nderiv, 0:ai_p, 1:ai_dirnpts) = &
        self % opr_dBasisatx (ai_d, 0:ai_nderiv, 0:ai_p, 1:ai_dirnpts, ai_ielt)

!        print *, "====== ielt ", ai_ielt, " ======"
!        DO li_pt = 1,ai_dirnpts
!            print *, 'pt = ', li_pt
!            DO li_d = 0,ai_nderiv
!                print *, '----> derivative = ', li_d
!                print *, ao_bbox % opr_dB(ai_d, li_d, 0:ai_p, li_pt)
!            END DO
!        END DO
!        print *, "============================="
#ifdef _TRACE
        CALL printlog("get_elts_t_tensor_basis : End", ai_dtllevel = 1)
#endif

    end subroutine get_elts_t_tensor_basis

end module elements_tensor_module




