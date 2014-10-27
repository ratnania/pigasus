!     
! File:   bbox_def.F90
! Author: ratnani
!
! Created on December 5, 2011, 2:13 PM
!

module bbox_def
    use used_precision
!    use assl_module
    implicit none

    private

    type, public :: BBOX
        integer :: oi_dim
        integer :: oi_nderiv
        integer :: oi_nderiv_code
        integer :: oi_nen
        integer :: oi_npts
        integer :: oi_tensorlevel ! 0 -> no tensor, 1-> element-tensor (Bezier surfaces), 2-> patch-tensor (IGA approach)

        integer :: oi_maxnpts
        integer :: oi_maxp
    ! ***************************************
    !       TENSOR LEVEL = 1
    ! ***************************************
    ! ...
        !> \todo peut etre a supprimer d'ici, et les mettre en local si besoin
        integer, dimension(:), pointer :: opi_leftmk
        real(wp), dimension(:,:,:,:), pointer :: opr_dB
    ! ...        
        real(wp), dimension(:,:,:), pointer :: opr_dBasis
    ! ***************************************

!    ! ***************************************
!    !       TENSOR LEVEL = 2
!    ! ***************************************
!      type(ASSL_DATA) :: oo_assl  ! needed to assemble the basis
!      type(ASSL_DATA) :: oo_asslx ! needed to assemble the geometry      
!      real(8), dimension(:,:), pointer :: opr_r
!      ! lbasis(der,ij, elt) est la valeur de la i eme spline sur le j eme point de quadrature
!      ! dans l'element elt: elt est un element appartenant a l'orthogonal de la direction 1
!      ! nel = n2 ou n2 * n3
!      real(8), dimension(:,:,:), pointer :: opr_lBasis
!      ! lbasis_t(der,d,iiprimej, elt) est PRODUCT(N_i(gl_j) * N_iprime(gl_j))
!      ! dans l'element elt: elt est un element appartenant a l'orthogonal de la direction 1
!      ! nel = n2 ou n2 * n3
!      real(8), dimension(:,:,:), pointer :: opr_lBasis_t
!      ! gbasis(der,ij) est la valeur de la i eme spline sur le j eme point de quadrature
!      ! en parcourant tous les elements de la direction 1
!      real(8), dimension(:,:), pointer :: opr_gBasis
!      ! gbasis_t(der,ij) est N_i(gl_j) * N_iprime(gl_j)
!      ! en parcourant tous les elements de la direction 1
!      real(8), dimension(:,:), pointer :: opr_gBasis_t
!    ! ***************************************
    end type BBOX

    type, public :: BBOX_DIAGS
        integer :: oi_dim
        integer :: oi_nderiv
        integer :: oi_nderiv_code
        integer :: oi_nen
        integer :: oi_npts
        integer :: oi_tensorlevel

        integer :: oi_maxnpts
        integer :: oi_maxp
    ! ...
        !> \todo peut etre a supprimer d'ici, et les mettre en local si besoin
        integer, dimension(:,:), pointer :: opi_leftmk
        real(wp), dimension(:,:,:,:), pointer :: opr_dB
    ! ...
        real(wp), dimension(:,:,:), pointer :: opr_dBasis
    end type BBOX_DIAGS

end module bbox_def

