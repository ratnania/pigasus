!
! File:   grids_def.F90
! Author: root
!
! Created on January 2, 2012, 2:14 PM
!

module grids_def
    use used_precision
    implicit none

    !******************************************************************************************
    !	STORED structure
    !   WE STORE IN THIS STRUCTURE THE DERIVATION OF ALL B-SPLINES IN EACH DIRECTION
    !   THE GOLBAL MEMORY COST (FOR 3D) IS (Nx + Ny + Nz) * ( nderiv + 1 )
    !******************************************************************************************
    type, public :: STORED_DATA

        INTEGER :: oi_nderiv

        ! opr_dB(li_der, li_N, 0:li_p, li_i)
        real(wp), dimension(:,:,:,:), pointer :: opr_dB

    end type STORED_DATA
    !******************************************************************************************

    !******************************************************************************************
    !	element structure
    !******************************************************************************************
    type, public :: ELEMENT_DATA

        !> default False
        !> True if we have computed and stored all desired values for the current
        logical :: ol_stored

        !> number of points inside the element
        integer :: oi_npts

        !> number of points inside each 1-direction element
        integer, dimension(:), pointer :: opi_npts
        
        !> evaluation points where we need to evaluate all splines : quadrature points.
        !> opr_pts ( i , : ) are the coordinates of the i^th point
        real(wp), dimension (:,:), pointer :: opr_pts
        !> weights for quadrature points.
        real(wp), dimension (:), pointer :: opr_w
        !> default = TRUE
        !> will be set on False, if the pts belongs to many elements,
        !> we then will activate only one
        logical, dimension (:), pointer :: opl_activ_pts
        !> will be used for local refinement
        logical :: ol_activ

        !> this array is used to store any values of fields defined on the current grid
        !> the size of this array must be equal to the number of fields defined
        !> on the current space
        !> opr_values(li_field, li_dof, li_pt) ; dof to handle scalar/vector fields
        real(wp), dimension (:,:,:), pointer :: opr_values_f

        !> this array is used to store any values of matrices defined on the current grid
        !> the size of this array must be equal to the number of matrices defined
        !> on the current space
        !> opr_values(li_matrix, li_param, li_pt) ; 
        real(wp), dimension (:,:,:), pointer :: opr_values_m

        !> this array is used to store any values of norms defined on the current grid
        !> the size of this array must be equal to the number of norms defined
        !> on the current space
        !> opr_values(li_norm, li_param, li_pt) ; 
        real(wp), dimension (:,:,:), pointer :: opr_values_n
    end type ELEMENT_DATA
    !******************************************************************************************

    !******************************************************************************************
    !	GRID structure
    !******************************************************************************************
    type, public :: GRID_DATA

        !> used only for tensor elements, the maximum number of elements per direction
        integer :: oi_maxdir_npts

        !> ol_tensor = .TRUE. if we want to save data as a tensor,
        !> we will only save the values for each direction,
        !> this is very important for fast evaluations
        !> default value = .FALSE.
        integer :: oi_tensorlevel
        !> IF WE USE FIELDS VALUES
        logical :: ol_use_values_f
        !> IF WE USE MATRICES VALUES
        logical :: ol_use_values_m
        !> IF WE USE NORMS VALUES
        logical :: ol_use_values_n

        logical :: ol_elts_alloc
        INTEGER :: oi_ngrid
        INTEGER :: oi_nel
        INTEGER :: oi_dim
        !> dimension of field : scalar/vector
        INTEGER :: oi_dof
        INTEGER :: oi_nfields
        INTEGER :: oi_nmatrices
        INTEGER :: oi_nnorms

        type(ELEMENT_DATA), dimension(:), pointer :: opo_elts

        !> IF WE NEED TO STORE SOME COMPUTATIONS
        !> THIS IS DONE ONLY IS THE CASE WHEN TENSOR PRODUCT IS ACTIVATED
        logical :: ol_usestored_data
        
        type(STORED_DATA), dimension(:), pointer :: opo_stored

        !> used to store the values of bsplines and their derivatives for each direction
        ! opr_dBasisatx (ai_d, 0:ai_nderiv, 0:ai_p, 1:ai_dirnpts, ai_ielt)
        real(wp), dimension(:,:,:,:,:), allocatable :: opr_dBasisatx

    end type GRID_DATA
    !******************************************************************************************

    !******************************************************************************************
    type, public :: GRIDS_DATA
        !> the number of patchs
        integer :: oi_npatchs
        !> for each patch we must have the appriate grid-data
        type(GRID_DATA), dimension(:), pointer :: opo_grid
    end type GRIDS_DATA
    !******************************************************************************************

    contains
        !---------------------------------------------------------------
        subroutine from_code_to_nderiv(ai_dim, ai_nderiv_code, api_deriv)
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

        end subroutine from_code_to_nderiv
        !---------------------------------------------------------------
        subroutine from_nderiv_to_code(ai_dim, api_deriv, ai_nderiv_code)
            !> compute the corresponding deriv indices
            !> to put in a new file grids_tools
            implicit none
            integer, intent(in) :: ai_dim
            integer, dimension (:), intent(in) :: api_deriv
            integer, intent(out) :: ai_nderiv_code
            ! LOCAL

            if (ai_dim == 1) then
                ai_nderiv_code = api_deriv ( 1 )
                return
            end if

            if (ai_dim == 2) then
                IF ( (api_deriv ( 1 ) == 0) .AND. (api_deriv ( 2 ) == 0) ) THEN
                    ai_nderiv_code = 0
                END IF
                IF ( (api_deriv ( 1 ) == 1) .AND. (api_deriv ( 2 ) == 0) ) THEN
                    ai_nderiv_code = 1
                END IF
                IF ( (api_deriv ( 1 ) == 0) .AND. (api_deriv ( 2 ) == 1) ) THEN
                    ai_nderiv_code = 2
                END IF
                IF ( (api_deriv ( 1 ) == 2) .AND. (api_deriv ( 2 ) == 0))  THEN
                    ai_nderiv_code = 3
                END IF
                IF ( (api_deriv ( 1 ) == 1) .AND. (api_deriv ( 2 ) == 1) ) THEN
                    ai_nderiv_code = 4
                END IF
                IF ( (api_deriv ( 1 ) == 0) .AND. (api_deriv ( 2 ) == 2) ) THEN
                    ai_nderiv_code = 5
                END IF
            end if

            if (ai_dim == 3) then
                IF ( (api_deriv ( 1 ) == 0) .AND. (api_deriv ( 2 ) == 0) .AND. (api_deriv ( 3 ) == 0))  THEN
                    ai_nderiv_code = 0
                END IF
                IF ( (api_deriv ( 1 ) == 1) .AND. (api_deriv ( 2 ) == 0) .AND. (api_deriv ( 3 ) == 0))  THEN
                    ai_nderiv_code = 1
                END IF
                IF ( (api_deriv ( 1 ) == 0) .AND. (api_deriv ( 2 ) == 1) .AND. (api_deriv ( 3 ) == 0))  THEN
                    ai_nderiv_code = 2
                END IF
                IF ( (api_deriv ( 1 ) == 0) .AND. (api_deriv ( 2 ) == 0) .AND. (api_deriv ( 3 ) == 1))  THEN
                    ai_nderiv_code = 3
                END IF
                IF ( (api_deriv ( 1 ) == 2) .AND. (api_deriv ( 2 ) == 0) .AND. (api_deriv ( 3 ) == 0))  THEN
                    ai_nderiv_code = 4
                END IF
                IF ( (api_deriv ( 1 ) == 1) .AND. (api_deriv ( 2 ) == 1) .AND. (api_deriv ( 3 ) == 0))  THEN
                    ai_nderiv_code = 5
                END IF
                IF ( (api_deriv ( 1 ) == 1) .AND. (api_deriv ( 2 ) == 0) .AND. (api_deriv ( 3 ) == 1) ) THEN
                    ai_nderiv_code = 6
                END IF
                IF ( (api_deriv ( 1 ) == 0) .AND. (api_deriv ( 2 ) == 2) .AND. (api_deriv ( 3 ) == 0))  THEN
                    ai_nderiv_code = 7
                END IF
                IF ( (api_deriv ( 1 ) == 0) .AND. (api_deriv ( 2 ) == 1) .AND. (api_deriv ( 3 ) == 1))  THEN
                    ai_nderiv_code = 8
                END IF
                IF ( (api_deriv ( 1 ) == 0) .AND. (api_deriv ( 2 ) == 0) .AND. (api_deriv ( 3 ) == 2))  THEN
                    ai_nderiv_code = 9
                END IF
            end if

        end subroutine from_nderiv_to_code
end module grids_def