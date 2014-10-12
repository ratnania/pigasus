!     
! File:   geometries_def.F90
! Author: ratnani
!
! Created on December 5, 2011, 2:13 PM
!

module geometries_def
    use used_precision
    implicit none

    private

    type, public :: GEOMETRY
        !> dimension of the geometry : 1D, 2D, 3D, .. nD
        integer :: oi_dim
        !> 1 if we use nurbs 0 otherwise
        integer :: oi_rational
        !> dimension of the Points : 1D, 2D, 3D, .. nD
        integer :: oi_Rd
        !> the number of the Points
        integer :: oi_npts
        !> list of the number of control points for all directions
        !> this array is of dimension oi_dim
        integer, dimension(:), pointer :: opi_N
        integer, dimension(:), pointer :: opi_Nw
        !> list of the splines degrees for all directions
        !> this array is of dimension oi_dim
        integer, dimension(:), pointer :: opi_P
        !> list of knot vectors for all dimensions
        !> the first index is allocated with the maximum of N+P+1
        !> the second index is for the direction
        real(wp), dimension(:,:), pointer :: opr_u
        !> list of control points ordered in 1D :
        !> if oi_dim = 3 then the index of the (i,j,k) control point is
        !> index = ( (k-1)N[2] + j - 1 ) N[1] + i
        real(wp), dimension(:,:), pointer :: opr_P
        !> weights
        real(wp), dimension(:), pointer :: opr_W
        real(wp), dimension(:,:,:), pointer :: opr_Ww
        !> the number of different knots for each direction
        integer, dimension(:), pointer :: opi_ndiffu
        !> the correspondant index without repeating knots
        real(wp), dimension(:,:), pointer :: opi_Real_Index
        !> indices defining the element
        real(wp), dimension(:,:), pointer :: opi_ELT_INDEX
        !>
        real(wp), dimension(:), pointer :: opi_elt
        !> TRUE if we have allocated memory for ELt_INDEX and ELT
        logical :: ol_elt_allocated
        !> in opi_bi we store for each patch, element and local basis id, the local 1D direction id,
        !> this works only if we create space with the tensor property activated
        integer, dimension(:,:,:), pointer :: opi_bi
        !> in opi_ni we store for each patch and element, the local 1D direction id of the element,
        !> this works only if we create space with the tensor property activated
        integer, dimension(:,:), pointer :: opi_ni
        !> additional weights used to normalize splines (Maxwell's equations for ex)
        !> indices are like for the real knots
        real(wp), dimension(:,:), pointer :: opr_L
        !*************
        ! tensor-level 2
        !*************
        !> for each direction Rd and each element of the hyperplan
        !> we store for each element in the 1st direction
        !> the relevant control points for the band-line element
        real(8), dimension(:,:,:), pointer :: opr_X
        !*************
    end type GEOMETRY

    type, public :: GEOMETRY_INFO
        integer :: oi_nen
    end type GEOMETRY_INFO
    
    type, public :: GEOMETRIES
        integer :: oi_ngeo
        type(GEOMETRY), dimension(:), pointer :: opo_geo
    end type GEOMETRIES

end module geometries_def

