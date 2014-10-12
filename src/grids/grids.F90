module grids
    use tracelog_module
    use bbox_def
    use geometries_def
!    use geometry_module
!    use geometries_module
!    use assl_module
    use elements_module
    use elements_tensor_module
    use grids_def
    implicit none

#ifdef _DEBUG
    integer, parameter, private  :: mi_dtllevel_base = 0
#else
    integer, parameter, private  :: mi_dtllevel_base = 2
#endif
    
    contains
    !---------------------------------------------------------------
    subroutine create_elements( self, ai_dim, ai_nel, ai_dof    &
        , ai_nfields &
        , ai_nmatrices &
        , ai_nnorms &
!        , api_npts &
        , ai_npts &
        , ai_dirnpts, ai_tensorlevel &
        , ai_maxnparams_matrices &
        , ai_maxnparams_fields &
        , ai_maxnparams_norms )
        implicit none
        type(GRID_DATA) :: self
        INTEGER :: ai_dim
        INTEGER :: ai_nel
        INTEGER, optional, intent(in)  :: ai_dof
        INTEGER, optional, intent(in)  :: ai_nfields
        INTEGER, optional, intent(in)  :: ai_nmatrices
        INTEGER, optional, intent(in)  :: ai_nnorms
!        INTEGER, dimension(ai_nel), optional, intent(in) :: api_npts
        INTEGER, optional, intent(in)  :: ai_npts
        INTEGER, optional, intent(in)  :: ai_dirnpts
        INTEGER, optional, intent(in)  :: ai_tensorlevel
        INTEGER, optional, intent(in)  :: ai_maxnparams_matrices
        INTEGER, optional, intent(in)  :: ai_maxnparams_fields
        INTEGER, optional, intent(in)  :: ai_maxnparams_norms
        ! LOCAL
        INTEGER, dimension(ai_nel) :: lpi_npts
        INTEGER :: li_totalnpts
#ifdef _TRACE
        CALL printlog("create_elements : Start", ai_dtllevel = mi_dtllevel_base + 1)
#endif
        
        self % oi_tensorlevel = 1
        if ( present(ai_tensorlevel)) then
            self % oi_tensorlevel = ai_tensorlevel
        end if
        
            if ( self % oi_tensorlevel /= 0) then
                if ( present(ai_dirnpts)) then
                    lpi_npts = ai_dirnpts
!                    print*,'ddddd'
                else
                    print *, 'SERIOUS ERROR in create_elements: you must give either ai_dirnpts'
                end if
                li_totalnpts = 0
                if ( present(ai_npts)) then
                    li_totalnpts = ai_npts
!                    print*,'xxxx'
                else
                    print *, 'SERIOUS ERROR in create_elements: you must give either ai_npts needed for li_totalnpts'
                end if
            else
                if ( present(ai_npts)) then
                    lpi_npts = ai_npts
                else
                    print *, 'SERIOUS ERROR in create_elements: you must give either ai_npts'
                end if
            end if
!        end if
        
        self % ol_usestored_data = .FALSE. ! default treatment is to not store data
        
        self % oi_nel = ai_nel
        self % oi_dim = ai_dim
        self % oi_dof = ai_dof

        ! ...
        self % oi_nfields = 0
        if (present(ai_nfields)) then
            self % oi_nfields = ai_nfields
        end if
        
        self % ol_use_values_f = .FALSE.
        if ( self % oi_nfields > 0 ) then
            self % ol_use_values_f = .TRUE.
        end if
        ! ...

        ! ...
        self % oi_nmatrices = 0
        if (present(ai_nmatrices)) then
            self % oi_nmatrices = ai_nmatrices
        endif

        self % ol_use_values_m = .FALSE.
        if ( self % oi_nmatrices > 0 ) then
            self % ol_use_values_m = .TRUE.
        end if
        ! ...

        ! ...
        self % oi_nnorms = 0
        if (present(ai_nnorms)) then
            self % oi_nnorms = ai_nnorms
        endif

        self % ol_use_values_n = .FALSE.
        if ( self % oi_nnorms > 0 ) then
            self % ol_use_values_n = .TRUE.
        end if
        ! ...

        self % oi_ngrid = SUM(lpi_npts)
!        print*,'li_totalnpts=', li_totalnpts
!        print*,'lpi_npts=', lpi_npts
        if (self % oi_tensorlevel /= 0) then
            CALL create_elts_t_elements(self, lpi_npts, li_totalnpts &
            , ai_maxnparams_matrices, ai_maxnparams_fields, ai_maxnparams_norms)
        else
            CALL create_elts_elements  (self, lpi_npts)
        end if
#ifdef _TRACE
        CALL printlog("create_elements : Start", ai_dtllevel = mi_dtllevel_base + 1)
#endif

    end subroutine create_elements
    !---------------------------------------------------------------
    subroutine free_elements( self )
        implicit none
        type(GRID_DATA) :: self
        ! LOCAL
#ifdef _TRACE
        CALL printlog("free_elements : Start", ai_dtllevel = mi_dtllevel_base + 1)
#endif

        if (self % oi_tensorlevel /= 0) then
            CALL free_elements_elts_t(self)
        else
            CALL free_elements_elts(self)
        end if
#ifdef _TRACE
        CALL printlog("free_elements : End", ai_dtllevel = mi_dtllevel_base + 1)
#endif

    end subroutine free_elements
    !---------------------------------------------------------------
    subroutine print_grid( self )
        implicit none
        type(GRID_DATA) :: self
        ! LOCAL
#ifdef _TRACE
        CALL printlog("print_grid : Start", ai_dtllevel = mi_dtllevel_base + 1)
#endif

        if (self % oi_tensorlevel /= 0) then
            CALL print_grid_elts_t(self)
        else
            CALL print_grid_elts(self)
        end if
#ifdef _TRACE
        CALL printlog("print_grid : End", ai_dtllevel = mi_dtllevel_base + 1)
#endif

    end subroutine print_grid
    !---------------------------------------------------------------
    subroutine get_maxnpoints_grid( self, ai_maxnpts )
        implicit none
        type(GRID_DATA) :: self
        INTEGER, INTENT(INOUT) :: ai_maxnpts
        ! LOCAL
        INTEGER :: li_elt
#ifdef _TRACE
        CALL printlog("get_maxnpoints_grid : Start", ai_dtllevel = mi_dtllevel_base + 2)
#endif

        ai_maxnpts = self % opo_elts (1) % oi_npts
        DO li_elt = 2, self % oi_nel
            IF (ai_maxnpts < self % opo_elts (li_elt) % oi_npts) THEN
                ai_maxnpts = self % opo_elts (li_elt) % oi_npts
            END IF
        END DO
#ifdef _TRACE
        CALL printlog("get_maxnpoints_grid : End", ai_dtllevel = mi_dtllevel_base + 2)
#endif
        
    end subroutine get_maxnpoints_grid
    !----------------------------------------------------------------------------------------------
    subroutine create_from_tensor_product_1D( self, apr_localx, apr_localw)
        implicit none
        type(GRID_DATA) :: self
        !> apr_localx contains for each element a local grid
        !> access : apr_localx ( li_elt, li_pt )
        !> in apr_localx ( li_elt, 0 ) we must put the number of pts inside the current element
        real(wp), dimension(:,0:) :: apr_localx
        !> local weight
        real(wp), dimension(:,:) :: apr_localw
#ifdef _TRACE
        CALL printlog("create_from_tensor_product_1D : Start", ai_dtllevel = mi_dtllevel_base + 1)
#endif

        if (self % oi_tensorlevel/=0) then
            CALL create_elts_t_from_tensor_product_1D(self,apr_localx, apr_localw)
        else
            CALL create_elts_from_tensor_product_1D(self,apr_localx, apr_localw)
        end if
#ifdef _TRACE
        CALL printlog("create_from_tensor_product_1D : End", ai_dtllevel = mi_dtllevel_base + 1)
#endif

    end subroutine create_from_tensor_product_1D
    !----------------------------------------------------------------------------------------------
    subroutine create_from_tensor_product_2D( self, apr_localx, apr_localy, apr_localwx, apr_localwy)
        implicit none
        type(GRID_DATA) :: self
        !> apr_localx contains for each element a local grid
        !> access : apr_localx ( li_elt, li_pt )
        !> in apr_localx ( li_elt, 0 ) we must put the number of pts inside the current element
        real(wp), dimension(:,0:) :: apr_localx
        real(wp), dimension(:,0:) :: apr_localy
        real(wp), dimension(:,:) :: apr_localwx
        real(wp), dimension(:,:) :: apr_localwy
#ifdef _TRACE
        CALL printlog("create_from_tensor_product_2D : Start", ai_dtllevel = mi_dtllevel_base + 1)
#endif

        if (self % oi_tensorlevel/=0) then
            CALL create_elts_t_from_tensor_product_2D(self,apr_localx, apr_localy    &
            , apr_localwx, apr_localwy)
        else
            CALL create_elts_from_tensor_product_2D(self,apr_localx, apr_localy  &
            , apr_localwx, apr_localwy)
        end if
#ifdef _TRACE
        CALL printlog("create_from_tensor_product_2D : End", ai_dtllevel = mi_dtllevel_base + 1)
#endif

    end subroutine create_from_tensor_product_2D
    !----------------------------------------------------------------------------------------------
    subroutine create_from_tensor_product_3D( self, apr_localx, apr_localy, apr_localz, apr_localwx, apr_localwy, apr_localwz)
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
#ifdef _TRACE
        CALL printlog("create_from_tensor_product_3D : Start", ai_dtllevel = mi_dtllevel_base + 1)
#endif

        if (self % oi_tensorlevel/=0) then
            CALL create_elts_t_from_tensor_product_3D(self,apr_localx, apr_localy, apr_localz    &
            , apr_localwx, apr_localwy, apr_localwz)
        else
            CALL create_elts_from_tensor_product_3D(self,apr_localx, apr_localy, apr_localz  &
            , apr_localwx, apr_localwy, apr_localwz)
        end if
#ifdef _TRACE
        CALL printlog("create_from_tensor_product_3D : End", ai_dtllevel = mi_dtllevel_base + 1)
#endif

    end subroutine create_from_tensor_product_3D
    !----------------------------------------------------------------------------------------------
    subroutine set_unidirection_tensor_2d( self, ai_direction   &
        , apr_dirpts, apr_dirw  &
        , ar_fixpt, ar_fixw, ai_elt)
        implicit none
        type(GRID_DATA) :: self
        !> apr_localx contains for each element a local grid
        !> access : apr_localx ( li_elt, li_pt )
        !> in apr_localx ( li_elt, 0 ) we must put the number of pts inside the current element
        integer :: ai_direction
        real(8), dimension(:,0:) :: apr_dirpts
        real(8), dimension(:,:) :: apr_dirw
        real(8) :: ar_fixpt
        real(8) :: ar_fixw
        !> ai_elt is the id of the last created element
        integer, intent(inout) :: ai_elt
#ifdef _TRACE
        CALL printlog("set_unidirection_tensor_2d : Start", ai_dtllevel = mi_dtllevel_base + 1)
#endif

        if (self % oi_tensorlevel==1) then
            CALL set_elts_t_unidirection_tensor_2d( self, ai_direction   &
            , apr_dirpts, apr_dirw  &
            , ar_fixpt, ar_fixw, ai_elt)
        else
            print *, 'Not done yet'
        end if
#ifdef _TRACE
        CALL printlog("set_unidirection_tensor_2d : End", ai_dtllevel = mi_dtllevel_base + 1)
#endif

    end subroutine set_unidirection_tensor_2d
    !---------------------------------------------------------------
!    subroutine assembly_elements( self, ao_geos, ai_nderiv, ai_id, ai_elt, apr_points)
!        implicit none
!        type(GRID_DATA) :: self
!        type(GEOMETRIES) :: ao_geos
!        integer :: ai_nderiv
!        integer :: ai_id !id of the geometry (in multi-patchs)
!        integer :: ai_elt ! the id of the element
!        real(wp), dimension(:,:,:), intent(out) :: apr_points !(ai_nderiv,ai_Rd,ai_npts)
!        ! LOCAL
!
!        CALL printlog("assembly_elements : Start", ai_dtllevel = mi_dtllevel_base + 1)
!
!        if (self % oi_tensorlevel) then
!            CALL assembly_elts_t(self, ao_geos, ai_nderiv, ai_id, ai_elt, apr_points)
!        else
!            CALL assembly_elts(self, ao_geos, ai_nderiv, ai_id, ai_elt, apr_points)
!        end if
!
!        CALL printlog("assembly_elements : End", ai_dtllevel = mi_dtllevel_base + 1)
!
!    end subroutine assembly_elements
    !---------------------------------------------------------------
    subroutine assembly_points_elements( self, ao_geos, ao_bbox, ai_id, ai_elt, ai_realelt, ai_npts, apr_points, ai_ptw_evaluation)
        implicit none
        type(GRID_DATA) :: self
        type(GEOMETRIES) :: ao_geos
        type(BBOX), intent(in) :: ao_bbox
        integer :: ai_id !id of the geometry (in multi-patchs)
        integer :: ai_elt ! the id of the element
        integer :: ai_realelt ! the id of the real element
        integer :: ai_npts
        real(wp), dimension(0:,:,:), intent(out) :: apr_points !(ai_nderiv,ai_Rd,ai_npts)
        integer, optional :: ai_ptw_evaluation
        ! LOCAL
        integer :: li_ptw_evaluation
        integer :: li_d
        integer :: li_der
        integer :: li_dcode
        integer :: li_dim
        integer :: li_nderiv_code
        integer :: li_local_dercode
        integer :: li_nnz
        integer :: li_annz
        integer :: li_rnnz
        integer :: li_blocn
        INTEGER, DIMENSION(:), POINTER  :: lpi_deriv
        real(wp), dimension(:), pointer :: lpr_b
        real(wp), dimension(:), pointer :: lpr_r
        real(wp), dimension(:), pointer :: lpr_x
        real(wp), dimension(:), pointer :: lpr_a

#ifdef _TRACE
        CALL printlog("assembly_points_elements : Start", ai_dtllevel = mi_dtllevel_base + 1)
#endif

        li_ptw_evaluation = ai_ptw_evaluation
        
        IF (self % oi_tensorlevel==1) THEN

        CALL assembly_points_elts_t(self, ao_geos, ai_id, ai_elt, ai_realelt, ai_npts &
        , ao_bbox, apr_points)
            
        ELSE
            PRINT *, 'assembly_points_elements : Not yet done'
        END IF

#ifdef _TRACE
        CALL printlog("assembly_points_elements : End", ai_dtllevel = mi_dtllevel_base + 1)
#endif

    end subroutine assembly_points_elements
    !---------------------------------------------------------------
    subroutine assembly_logical_basis( self, ao_geos, ai_nderiv, ai_id, ai_elt, ai_realelt, ao_bbox &
        , ai_nen &
        , ai_ptw_evaluation)
        implicit none
        type(GRID_DATA) :: self
        type(GEOMETRIES) :: ao_geos
        integer :: ai_nderiv
        integer :: ai_id !id of the geometry (in multi-patchs)
        integer :: ai_elt ! the id of the element
        integer :: ai_realelt ! the id of the real element
        type(BBOX), intent(inout) :: ao_bbox
        integer, optional :: ai_nen
        integer, optional :: ai_ptw_evaluation
        ! LOCAL
        real(wp), dimension(:,:), pointer :: lpr_sites
        integer :: li_Rd
        integer :: li_npts
        integer :: li_dim
        integer :: li_ptw_evaluation
#ifdef _TRACE
        CALL printlog("assembly_logical_basis : Start", ai_dtllevel = mi_dtllevel_base + 1)
#endif

        li_ptw_evaluation = ai_ptw_evaluation

        IF (self % oi_tensorlevel==1) THEN

            CALL assembly_basis_elts_t(self, ao_geos, ai_nderiv, ai_id, ai_elt, ai_realelt, ao_bbox, ai_nen)

        ELSE
            PRINT *, 'NOT YET IMPLEMENTED'
        END IF
#ifdef _TRACE
        CALL printlog("assembly_logical_basis : End", ai_dtllevel = mi_dtllevel_base + 1)
#endif

    end subroutine assembly_logical_basis
    !---------------------------------------------------------------
    subroutine set_values_field( self, ai_field, apr_values)
        implicit none
        type(GRID_DATA) :: self
        integer :: ai_field
        !> \param[inout] apr_values dimension(ai_ndof, ai_nparam, ai_npts)
        real(wp), dimension(:,:,:) :: apr_values
        ! LOCAL
        integer :: li_elt
        integer :: li_npts
        integer :: li_ndof
        integer :: li_nparam

        li_nparam = SIZE (apr_values,1)
#ifdef _TRACE
        CALL printlog("set_values_field : Start", ai_dtllevel = mi_dtllevel_base + 1)
#endif

#ifdef _DEBUG
print *, 'SIZE(apr_values,1)=', SIZE(apr_values,1)
print *, 'SIZE(apr_values,2)=', SIZE(apr_values,2)
print *, 'SIZE(apr_values,3)=', SIZE(apr_values,3)

print *, 'SIZE(self % opo_elts (1) % opr_values_f,1)=', SIZE(self % opo_elts (1) % opr_values_f,1)
print *, 'SIZE(self % opo_elts (1) % opr_values_f,2)=', SIZE(self % opo_elts (1) % opr_values_f,2)
print *, 'SIZE(self % opo_elts (1) % opr_values_f,3)=', SIZE(self % opo_elts (1) % opr_values_f,3)
#endif

!        PRINT *, 'FIELD  ', ai_field
!        PRINT *, 'values ', apr_values (1, :, :)

        li_ndof = SIZE(apr_values,1)
        DO li_elt = 1, self % oi_nel

            li_npts = self % opo_elts (li_elt) % oi_npts

            self % opo_elts (li_elt) % opr_values_f(ai_field, 1 : li_nparam, 1 : li_npts) =  &
            apr_values (1 : li_ndof, li_elt, 1 : li_npts )

        END DO
#ifdef _TRACE
        CALL printlog("set_values_field : End", ai_dtllevel = mi_dtllevel_base + 1)
#endif

    end subroutine set_values_field
    !---------------------------------------------------------------
    subroutine set_values_matrix( self, ai_matrix, apr_values)
        implicit none
        type(GRID_DATA) :: self
        integer :: ai_matrix
        real(wp), dimension(:,:,:) :: apr_values
        ! LOCAL
        integer :: li_elt
        integer :: li_npts
        integer :: li_nparam

        li_nparam = SIZE (apr_values,1)
#ifdef _TRACE
        CALL printlog("set_values_matrix : Start", ai_dtllevel = mi_dtllevel_base + 1)
#endif

#ifdef _DEBUG
call concatmsg("SIZE(apr_values,1) = ", ai_dtllevel = mi_dtllevel_base + 4)
call concatmsg(SIZE(apr_values,1), ai_dtllevel = mi_dtllevel_base + 4)

call concatmsg("SIZE(apr_values,2) = ", ai_dtllevel = mi_dtllevel_base + 4)
call concatmsg(SIZE(apr_values,2), ai_dtllevel = mi_dtllevel_base + 4)

call concatmsg("SIZE(apr_values,3) = ", ai_dtllevel = mi_dtllevel_base + 4)
call concatmsg(SIZE(apr_values,3), ai_dtllevel = mi_dtllevel_base + 4)

call concatmsg("SIZE(self % opo_elts (1) % opr_values_m,1) = ", ai_dtllevel = mi_dtllevel_base + 4)
call concatmsg(SIZE(self % opo_elts (1) % opr_values_m,1), ai_dtllevel = mi_dtllevel_base + 4)

call concatmsg("SIZE(self % opo_elts (1) % opr_values_m,2) = ", ai_dtllevel = mi_dtllevel_base + 4)
call concatmsg(SIZE(self % opo_elts (1) % opr_values_m,2), ai_dtllevel = mi_dtllevel_base + 4)

call concatmsg("SIZE(self % opo_elts (1) % opr_values_m,3) = ", ai_dtllevel = mi_dtllevel_base + 4)
call concatmsg(SIZE(self % opo_elts (1) % opr_values_m,3), ai_dtllevel = mi_dtllevel_base + 4)

call printmsg(ai_dtllevel = mi_dtllevel_base + 4)
#endif

        DO li_elt = 1, self % oi_nel

            li_npts = self % opo_elts (li_elt) % oi_npts

            self % opo_elts (li_elt) % opr_values_m(ai_matrix, 1 : li_nparam, 1 : li_npts) =  &
            apr_values (1 : li_nparam, li_elt, 1 : li_npts )

#ifdef _DEBUG
call concatmsg("li_elt = ", ai_dtllevel = mi_dtllevel_base + 4)
call concatmsg(li_elt, ai_dtllevel = mi_dtllevel_base + 4)

call printmsg(ai_dtllevel = mi_dtllevel_base + 4)

call concatmsg("SIZE(self % opo_elts (li_elt) % opr_values_m,1) = ", ai_dtllevel = mi_dtllevel_base + 4)
call concatmsg(SIZE(self % opo_elts (li_elt) % opr_values_m,1), ai_dtllevel = mi_dtllevel_base + 4)

call concatmsg("SIZE(self % opo_elts (li_elt) % opr_values_m,2) = ", ai_dtllevel = mi_dtllevel_base + 4)
call concatmsg(SIZE(self % opo_elts (li_elt) % opr_values_m,2), ai_dtllevel = mi_dtllevel_base + 4)

call concatmsg("SIZE(self % opo_elts (li_elt) % opr_values_m,3) = ", ai_dtllevel = mi_dtllevel_base + 4)
call concatmsg(SIZE(self % opo_elts (li_elt) % opr_values_m,3), ai_dtllevel = mi_dtllevel_base + 4)

call printmsg(ai_dtllevel = mi_dtllevel_base + 4)
#endif

!            print *, "self % opo_elts (li_elt) % opr_values_m(ai_matrix, 1 : li_nparam, 1 : li_npts)=", &
!            self % opo_elts (li_elt) % opr_values_m(ai_matrix, 1 : li_nparam, 1 : li_npts)

        END DO
#ifdef _TRACE
        CALL printlog("set_values_matrix : End", ai_dtllevel = mi_dtllevel_base + 1)
#endif

    end subroutine set_values_matrix
    !---------------------------------------------------------------
    subroutine set_values_norm( self, ai_norm, apr_values)
        implicit none
        type(GRID_DATA) :: self
        integer :: ai_norm
        real(wp), dimension(:,:,:) :: apr_values
        ! LOCAL
        integer :: li_elt
        integer :: li_npts
        integer :: li_nparam

        li_nparam = SIZE (apr_values,1)
#ifdef _TRACE
        CALL printlog("set_values_norm : Start", ai_dtllevel = mi_dtllevel_base + 1)
#endif

#ifdef _DEBUG
call concatmsg("SIZE(apr_values,1) = ", ai_dtllevel = mi_dtllevel_base + 4)
call concatmsg(SIZE(apr_values,1), ai_dtllevel = mi_dtllevel_base + 4)

call concatmsg("SIZE(apr_values,2) = ", ai_dtllevel = mi_dtllevel_base + 4)
call concatmsg(SIZE(apr_values,2), ai_dtllevel = mi_dtllevel_base + 4)

call concatmsg("SIZE(apr_values,3) = ", ai_dtllevel = mi_dtllevel_base + 4)
call concatmsg(SIZE(apr_values,3), ai_dtllevel = mi_dtllevel_base + 4)

call concatmsg("SIZE(self % opo_elts (1) % opr_values_m,1) = ", ai_dtllevel = mi_dtllevel_base + 4)
call concatmsg(SIZE(self % opo_elts (1) % opr_values_m,1), ai_dtllevel = mi_dtllevel_base + 4)

call concatmsg("SIZE(self % opo_elts (1) % opr_values_m,2) = ", ai_dtllevel = mi_dtllevel_base + 4)
call concatmsg(SIZE(self % opo_elts (1) % opr_values_m,2), ai_dtllevel = mi_dtllevel_base + 4)

call concatmsg("SIZE(self % opo_elts (1) % opr_values_m,3) = ", ai_dtllevel = mi_dtllevel_base + 4)
call concatmsg(SIZE(self % opo_elts (1) % opr_values_m,3), ai_dtllevel = mi_dtllevel_base + 4)

call printmsg(ai_dtllevel = mi_dtllevel_base + 4)
#endif

        DO li_elt = 1, self % oi_nel

            li_npts = self % opo_elts (li_elt) % oi_npts

            self % opo_elts (li_elt) % opr_values_n(ai_norm, 1 : li_nparam, 1 : li_npts) =  &
            apr_values (1 : li_nparam, li_elt, 1 : li_npts )

#ifdef _DEBUG
call concatmsg("li_elt = ", ai_dtllevel = mi_dtllevel_base + 4)
call concatmsg(li_elt, ai_dtllevel = mi_dtllevel_base + 4)

call printmsg(ai_dtllevel = mi_dtllevel_base + 4)

call concatmsg("SIZE(self % opo_elts (li_elt) % opr_values_m,1) = ", ai_dtllevel = mi_dtllevel_base + 4)
call concatmsg(SIZE(self % opo_elts (li_elt) % opr_values_m,1), ai_dtllevel = mi_dtllevel_base + 4)

call concatmsg("SIZE(self % opo_elts (li_elt) % opr_values_m,2) = ", ai_dtllevel = mi_dtllevel_base + 4)
call concatmsg(SIZE(self % opo_elts (li_elt) % opr_values_m,2), ai_dtllevel = mi_dtllevel_base + 4)

call concatmsg("SIZE(self % opo_elts (li_elt) % opr_values_m,3) = ", ai_dtllevel = mi_dtllevel_base + 4)
call concatmsg(SIZE(self % opo_elts (li_elt) % opr_values_m,3), ai_dtllevel = mi_dtllevel_base + 4)

call printmsg(ai_dtllevel = mi_dtllevel_base + 4)
#endif

!            print *, "self % opo_elts (li_elt) % opr_values_m(ai_norm, 1 : li_nparam, 1 : li_npts)=", &
!            self % opo_elts (li_elt) % opr_values_m(ai_norm, 1 : li_nparam, 1 : li_npts)

        END DO
#ifdef _TRACE
        CALL printlog("set_values_norm : End", ai_dtllevel = mi_dtllevel_base + 1)
#endif

    end subroutine set_values_norm
    !----------------------------------------------------------------------------------------------
    subroutine create_stored_data(self, ao_geos, ai_id, api_nderiv)
        implicit none
        type(GRID_DATA) :: self
        type(GEOMETRIES) :: ao_geos
        integer :: ai_id
        integer, dimension(:) :: api_nderiv
        ! LOCAL VARIABLES
        integer :: li_d
#ifdef _TRACE
        CALL printlog("create_stored_data : Start", ai_dtllevel = mi_dtllevel_base + 1)
#endif

        if (self % oi_tensorlevel==1) then

            CALL create_stored_data_elts_t(self, api_nderiv &
            , ao_geos % opo_geo(ai_id) % opi_N  &
            , ao_geos % opo_geo(ai_id) % opi_P )

        else
            print *, 'NOT YET IMPLEMENTED'
        end if
#ifdef _TRACE
        CALL printlog("create_stored_data : End", ai_dtllevel = mi_dtllevel_base + 1)
#endif

    end subroutine create_stored_data
    !----------------------------------------------------------------------------------------------
    subroutine free_stored_data(self)
        implicit none
        type(GRID_DATA) :: self
        ! LOCAL VARIABLES
        integer :: li_d
#ifdef _TRACE
        CALL printlog("free_stored_data : Start", ai_dtllevel = mi_dtllevel_base + 1)
#endif

        if (self % oi_tensorlevel==1) then

            CALL free_stored_data_elts_t(self)

        else
            print *, 'NOT YET IMPLEMENTED'
        end if
#ifdef _TRACE
        CALL printlog("free_stored_data : End", ai_dtllevel = mi_dtllevel_base + 1)
#endif

    end subroutine free_stored_data
    !---------------------------------------------------------------
    subroutine get_weights(self, apr_w)
        implicit none
        type(GRID_DATA) :: self
        real(wp), dimension(:,:), intent(inout) :: apr_w
        ! LOCAL
#ifdef _TRACE
        CALL printlog("get_weights : Start", ai_dtllevel = mi_dtllevel_base + 1)
#endif

        if (self % oi_tensorlevel==1) then

            CALL get_weights_elts_t(self, apr_w)

        else
            print *, 'NOT YET IMPLEMENTED'
        end if
#ifdef _TRACE
        CALL printlog("get_weights : End", ai_dtllevel = mi_dtllevel_base + 1)
#endif

    end subroutine get_weights
    !---------------------------------------------------------------
    subroutine set_tensor_basis ( self, ai_d, apr_dbasisatx, ai_nderiv, ai_p, ai_ni, ai_dirnpts &
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
#ifdef _TRACE
        CALL printlog("set_tensor_basis : Start", ai_dtllevel = 1)
#endif

        IF (self % oi_tensorlevel==1) THEN
            CALL set_elts_t_tensor_basis ( self, ai_d, apr_dbasisatx, ai_nderiv, ai_p, ai_ni, ai_dirnpts &
        , ai_maxnderiv, ai_maxp, ai_maxni, ai_maxdirnpts )
        END IF
#ifdef _TRACE
        CALL printlog("set_tensor_basis : End", ai_dtllevel = 1)
#endif

    end subroutine set_tensor_basis
end module grids
