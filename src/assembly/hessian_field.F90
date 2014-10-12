!     
! File:   second_deriv_s.F90
! Author: root
!
! Created on October 11, 2012, 9:05 AM
!


module hessian_field_module
use used_precision
use tracelog_module
use grids_def
use geometries_def
use geometry_tools
use bbox_def
use bbox_module
use connectivities_def
use fem_def
use assembly_def
implicit none

private

public :: build_hessian_field_local

#ifdef _DEBUG
integer, parameter, private  :: mi_dtllevel_base = 0
#else
integer, parameter, private  :: mi_dtllevel_base = 2
#endif
!> \todo to remove from here, and use an adequate variable
integer, parameter, private :: mi_MAXDIM = 3
integer, parameter, private :: mi_MAXDIM_invH = 3
integer, parameter, private :: mi_NPARAM = 3

contains
!----------------------------------------------------------------------------------------------
subroutine build_hessian_field_local(ao_ASS, ao_FEM, ai_grids_id, ai_id, ai_elt, ai_field, apr_fieldh)
implicit none
TYPE(ASSEMBLY) :: ao_ASS
type(FEM) :: ao_FEM
!> grids id
integer :: ai_grids_id
!> patch id
integer :: ai_id
!> current element
integer :: ai_elt
!> current field-id
integer :: ai_field
real*8, dimension(0:,:,:,:), intent(inout) :: apr_fieldh
! LOCAL VARIABLES
INTEGER :: li_dim

#ifdef _DEBUG
call printlog("build_hessian_field_local : Begin", ai_dtllevel = mi_dtllevel_base + 1)
#endif

li_dim = ao_FEM % opi_dim(ai_grids_id)

IF (li_dim == 1) THEN
        CALL printlog("build_hessian_field_local : Not Implemented in 1D", ai_dtllevel = 0)
ELSE IF (li_dim == 2) THEN
        CALL build_hessian_field_local_2D(ao_ASS, ao_FEM, ai_grids_id, ai_id, ai_elt, ai_field, apr_fieldh)
END IF

#ifdef _DEBUG
call printlog("build_hessian_field_local : End", ai_dtllevel = mi_dtllevel_base + 1)
#endif

end subroutine build_hessian_field_local
!----------------------------------------------------------------------------------------------
subroutine build_hessian_field_local_2D(ao_ASS, ao_FEM, ai_grids_id, ai_id, ai_elt, ai_field, apr_fieldh)
implicit none
TYPE(ASSEMBLY) :: ao_ASS
type(FEM) :: ao_FEM
!> grids id
integer :: ai_grids_id
!> patch id
integer :: ai_id
!> current element
integer :: ai_elt
!> current field-id
integer :: ai_field
real*8, dimension(0:,:,:,:), intent(inout) :: apr_fieldh
! LOCAL VARIABLES
integer :: li_npts
integer :: li_b1
integer :: li_b2
integer :: li_i
integer :: li_iprime
integer :: li_A1
integer :: li_A2
integer :: li_sp
integer :: li_bloc
integer :: li_ndof
integer :: li_map
INTEGER :: li_loc_id
INTEGER :: li_nparam
INTEGER :: li_field_id
INTEGER :: li_dim
INTEGER :: li_conelt
INTEGER :: li_pt
INTEGER :: li_iparam
real(wp) :: lr_contribution
INTEGER, PARAMETER :: MAXDIM = 3
real(wp), dimension(ao_ASS % oi_maxndof,MAXDIM,ao_ASS % oi_maxnpts) :: lpr_Basis
real(wp), dimension(ao_ASS % oi_maxnpts) :: lpr_jacobian
real(wp), DIMENSION(ao_ASS % oi_maxnpts) :: lpr_contribution
real(wp), dimension(ao_ASS % oi_maxnpts) :: lpr_w
REAL(wp), DIMENSION (mi_NPARAM, ao_ASS % oi_maxnpts) :: lpr_A
REAL(wp), DIMENSION (mi_MAXDIM, ao_ASS % oi_maxnpts) :: lpr_alpha
REAL(wp), DIMENSION (mi_MAXDIM, ao_ASS % oi_maxnpts) :: lpr_beta
REAL(wp), DIMENSION (mi_MAXDIM_invH, mi_MAXDIM_invH, ao_ASS % oi_maxnpts) :: lpr_invH
TYPE(CONNECTIVITY), pointer :: lp_con
TYPE(CONNECTIVITY), pointer :: lp_conprime
TYPE(GRID_DATA), pointer :: lp_grid

#ifdef _DEBUG
call printlog("build_hessian_field_local_2D : Begin", ai_dtllevel = mi_dtllevel_base + 1)
#endif
!print *, 'ai_elt=', ai_elt

li_dim = ao_FEM % opi_dim(ai_grids_id)
lpr_Basis = 0.0_wp

!\todo prblm du id => id+1
lp_grid => ao_FEM % opo_grids ( ai_grids_id ) % opo_grid ( ai_id-1 )
li_npts = lp_grid % opo_elts ( ai_elt ) % oi_npts

li_sp       = ao_FEM % opi_InfoField(ai_field, INFOFIELD_SPACE )
li_ndof     = ao_FEM % opi_InfoField(ai_field, INFOFIELD_NDOF )
li_loc_id   = ao_FEM % opi_InfoField(ai_field, INFOFIELD_LOCID )
li_nparam   = ao_FEM % opi_InfoField(ai_field, INFOFIELD_NPARAM )

li_field_id = ao_FEM % opi_InfoField(ai_field, INFOFIELD_OPERANDE )

lp_con          => ao_FEM % opo_spaces (li_sp) % oo_con

li_conelt = lp_con % opi_real_elts (ai_id-1,ai_elt)

lpr_jacobian = ao_ASS % opo_info_gr (ai_grids_id) % opr_jacobians
lpr_jacobian = dabs(lpr_jacobian)

!... we begin by computing the matrix H^{-1}
!        lpr_alpha (1:li_dim, 1:li_npts) = ao_ASS % opo_info_sp (li_sp) % opr_points (1, 1:li_dim, 1:li_npts)
!        lpr_beta  (1:li_dim, 1:li_npts) = ao_ASS % opo_info_sp (li_sp) % opr_points (2, 1:li_dim, 1:li_npts)
!
!        lpr_invH(1, 1, 1:li_npts) = lpr_beta(2,1:li_npts)**2
!        lpr_invH(1, 2, 1:li_npts) = -2.0_wp * lpr_beta(1,1:li_npts) * lpr_beta(2,1:li_npts)
!        lpr_invH(1, 3, 1:li_npts) = lpr_alpha(2,1:li_npts)**2
!
!        lpr_invH(2, 1, 1:li_npts) = -lpr_alpha(2,1:li_npts) * lpr_beta(2,1:li_npts)
!        lpr_invH(2, 2, 1:li_npts) = lpr_alpha(1,1:li_npts) * lpr_beta(2,1:li_npts) + lpr_alpha(2,1:li_npts) * lpr_beta(1,1:li_npts)
!        lpr_invH(2, 3, 1:li_npts) = -lpr_alpha(1,1:li_npts) * lpr_beta(1,1:li_npts)
!
!        lpr_invH(3, 1, 1:li_npts) = lpr_beta(1,1:li_npts)**2
!        lpr_invH(3, 2, 1:li_npts) = -2.0_wp * lpr_alpha(1,1:li_npts) * lpr_alpha(2,1:li_npts)
!        lpr_invH(3, 3, 1:li_npts) = lpr_alpha(1,1:li_npts)**2
!
!        DO li_pt = 1, li_npts
!!            lpr_invH (:,:,li_pt) = TRANSPOSE(lpr_invH (:,:,li_pt)) / (lpr_jacobian(li_pt)**2)
!            lpr_invH (:,:,li_pt) = lpr_invH (:,:,li_pt) / (lpr_jacobian(li_pt)**2)
!        END DO

lpr_alpha(1, 1:li_npts) = ao_ASS % opo_info_gr (ai_grids_id) % opr_invJacobian(2,2,1:li_npts)
lpr_alpha(2, 1:li_npts) = ao_ASS % opo_info_gr (ai_grids_id) % opr_invJacobian(2,1,1:li_npts)
lpr_beta (1, 1:li_npts) = ao_ASS % opo_info_gr (ai_grids_id) % opr_invJacobian(1,2,1:li_npts)
lpr_beta (2, 1:li_npts) = ao_ASS % opo_info_gr (ai_grids_id) % opr_invJacobian(1,1,1:li_npts)

!        print *, 'alpha_1 = ', lpr_alpha(1, 1:li_npts)
!        print *, 'alpha_2 = ', lpr_alpha(2, 1:li_npts)
!        print *, 'beta_1  = ', lpr_beta(1, 1:li_npts)
!        print *, 'beta_2  = ', lpr_beta(2, 1:li_npts)

lpr_invH(1, 1, 1:li_npts) = lpr_beta(2,1:li_npts)**2
lpr_invH(1, 2, 1:li_npts) = 2.0_wp * lpr_beta(1,1:li_npts) * lpr_beta(2,1:li_npts)
lpr_invH(1, 3, 1:li_npts) = lpr_alpha(2,1:li_npts)**2

lpr_invH(2, 1, 1:li_npts) = lpr_alpha(2,1:li_npts) * lpr_beta(2,1:li_npts)
lpr_invH(2, 2, 1:li_npts) = lpr_alpha(1,1:li_npts) * lpr_beta(2,1:li_npts) + lpr_alpha(2,1:li_npts) * lpr_beta(1,1:li_npts)
lpr_invH(2, 3, 1:li_npts) = lpr_alpha(1,1:li_npts) * lpr_beta(1,1:li_npts)

lpr_invH(3, 1, 1:li_npts) = lpr_beta(1,1:li_npts)**2
lpr_invH(3, 2, 1:li_npts) = 2.0_wp * lpr_alpha(1,1:li_npts) * lpr_alpha(2,1:li_npts)
lpr_invH(3, 3, 1:li_npts) = lpr_alpha(1,1:li_npts)**2

!        print *, 'invH = ', lpr_invH(:, :, 1)
!...

!...
! compute the param matrix
!...
lpr_A (1:li_nparam, 1:li_npts) = lp_grid % opo_elts (ai_elt) % opr_values_f(li_loc_id, 1:li_nparam, 1:li_npts)
!...

!> \todo attention +1
do li_b1 = 1, lp_con % opi_nen(ai_id)
do li_b2 = 1, lp_con % opi_nen(ai_id)

li_A1 = lp_con % opi_LM(ai_id, li_b1, li_conelt)
li_A2 = lp_con % opi_LM(ai_id, li_b2, li_conelt)

if (li_A1 == 0) then
        cycle
end if

if (li_A2 == 0) then
        cycle
end if

CALL af_Basis_2D(ao_ASS, ao_FEM, li_sp, li_ndof, ai_grids_id, ai_id, ai_elt, li_b1, li_b2, lpr_invH, lpr_A, lpr_Basis )

!> \todo traiter le cas ou ndof > 1
! apr_Basis (1:li_ndof,1:li_npts)
!> \todo on doit enlever 1 a cause des indices
!> bbox % dBasis commence a 0

apr_fieldh(ai_field,1,1, 1:li_npts) = apr_fieldh(ai_field,1,1, 1:li_npts) &
+ ao_FEM % opo_F(li_field_id) % opr_c(li_A1) &
* ao_FEM % opo_F(li_field_id) % opr_c(li_A2) &
* lpr_Basis (1,1,:)
!            print *, "base =", lpr_Basis (1,1,:)

end do
end do
!        print *, 'uh =', apr_fieldh(ai_field,1,1, 1:li_npts)

#ifdef _DEBUG
call printlog("build_hessian_field_local_2D : End", ai_dtllevel = mi_dtllevel_base + 1)
#endif
end subroutine build_hessian_field_local_2D
!---------------------------------------------------------------
subroutine af_Basis_2D(ao_ASS, ao_FEM, ai_space, ai_ndof, ai_grids_id, ai_id, ai_elt, ai_b1, ai_b2, apr_invH, apr_A, apr_Basis)
implicit none
TYPE(ASSEMBLY) :: ao_ASS
TYPE(FEM) :: ao_FEM
INTEGER :: ai_space
INTEGER :: ai_ndof
INTEGER :: ai_grids_id
INTEGER :: ai_id
INTEGER :: ai_elt
INTEGER :: ai_b1
INTEGER :: ai_b2
real(wp), DIMENSION(:,:,:) :: apr_invH
real(wp), DIMENSION(:,:) :: apr_A
real(wp), DIMENSION(:,:,:) :: apr_Basis
! LOCAL
TYPE(BBOX), POINTER :: lp_bbox
TYPE(GRID_DATA), pointer :: lp_grid
TYPE(GEOMETRIES), pointer :: lp_geos
TYPE(METRIC_INFO), pointer :: lp_info
real(wp), dimension(:,:,:), pointer :: lpr_invJacobian
REAL(wp), DIMENSION (ao_ASS % oi_maxndof, mi_MAXDIM, ao_ASS % oi_maxnpts) :: lpr_Basis
real(wp), dimension ( mi_MAXDIM_invH) :: lpr_rhs1
real(wp), dimension ( mi_MAXDIM_invH) :: lpr_rhs2
real(wp), dimension ( mi_MAXDIM_invH) :: lpr_d2_1
real(wp), dimension ( mi_MAXDIM_invH) :: lpr_d2_2
INTEGER :: li_d
INTEGER :: li_npts
INTEGER :: li_dim
INTEGER :: li_pt
REAL(wp), dimension (ao_ASS % oi_maxnpts) :: lpr_phi_s
REAL(wp), dimension (ao_ASS % oi_maxnpts) :: lpr_phi_t
REAL(wp), dimension (ao_ASS % oi_maxnpts) :: lpr_phi_ss
REAL(wp), dimension (ao_ASS % oi_maxnpts) :: lpr_phi_st
REAL(wp), dimension (ao_ASS % oi_maxnpts) :: lpr_phi_tt
REAL(wp), dimension (ao_ASS % oi_maxnpts) :: lpr_s_x
REAL(wp), dimension (ao_ASS % oi_maxnpts) :: lpr_s_y
REAL(wp), dimension (ao_ASS % oi_maxnpts) :: lpr_t_x
REAL(wp), dimension (ao_ASS % oi_maxnpts) :: lpr_t_y
REAL(wp), dimension (ao_ASS % oi_maxnpts) :: lpr_s_xx
REAL(wp), dimension (ao_ASS % oi_maxnpts) :: lpr_s_xy
REAL(wp), dimension (ao_ASS % oi_maxnpts) :: lpr_s_yy
REAL(wp), dimension (ao_ASS % oi_maxnpts) :: lpr_t_xx
REAL(wp), dimension (ao_ASS % oi_maxnpts) :: lpr_t_xy
REAL(wp), dimension (ao_ASS % oi_maxnpts) :: lpr_t_yy

#ifdef _DEBUG
call printlog("af_Basis_2D : Begin", ai_dtllevel = mi_dtllevel_base + 3)
#endif

li_dim = ao_FEM % opi_dim(ai_grids_id)

!\todo prblm du id => id+1
lp_grid => ao_FEM % opo_grids (ai_grids_id) % opo_grid ( ai_id-1 )

lpr_invJacobian => ao_ASS % opo_info_gr (ai_grids_id) % opr_invJacobian

li_npts = lp_grid % opo_elts ( ai_elt ) % oi_npts

lp_bbox => ao_ASS % opo_bbox_sp ( ai_space )
lp_geos => ao_FEM % opo_spaces ( ai_space ) % oo_mapping

DO li_pt=1, li_npts
CALL compute_rhs(ai_b1, li_pt, lpr_rhs1)
CALL compute_rhs(ai_b2, li_pt, lpr_rhs2)

!... compute the contribution
lpr_d2_1(1:mi_MAXDIM_invH)  = MATMUL(apr_invH(1:mi_MAXDIM_invH, 1:mi_MAXDIM_invH, li_pt), lpr_rhs1(1:mi_MAXDIM_invH))
lpr_d2_2(1:mi_MAXDIM_invH)  = MATMUL(apr_invH(1:mi_MAXDIM_invH, 1:mi_MAXDIM_invH, li_pt), lpr_rhs2(1:mi_MAXDIM_invH))

apr_Basis (1,1,li_pt) = lpr_d2_1(1) * lpr_d2_2(3) - lpr_d2_1(2) * lpr_d2_2(2)
END DO

#ifdef _DEBUG
call printlog("af_Basis_2D : End", ai_dtllevel = mi_dtllevel_base + 3)
#endif

contains
        subroutine compute_rhs(ai_b, ai_pt, apr_rhs)
        implicit none
        integer :: ai_b
        integer :: ai_pt
        real(wp), dimension (mi_MAXDIM_invH) :: apr_rhs


        lpr_Basis (1,1:li_dim,ai_pt) = MATMUL(lpr_invJacobian(1:li_dim, 1:li_dim, ai_pt) &
        , lp_bbox % opr_dBasis (1:li_dim,ai_b-1,ai_pt))

        !...
        !... COMPUTING THE RIGHT HAND SIDE : rhs = D2u_tilda - C
        !... where C = (alpha_11 phi_x + beta_11 phi_y, alpha_12 phi_x + beta_12 phi_y, alpha_22 phi_x + beta_22 phi_y )^T
        !...
        apr_rhs(1) = lp_bbox % opr_dBasis (3,ai_b-1,ai_pt) &
        - ( ao_ASS % opo_info_gr (ai_grids_id) % opr_points (3, 1, ai_pt) * lpr_Basis(1,1, ai_pt)  &
        +   ao_ASS % opo_info_gr (ai_grids_id) % opr_points (3, 2, ai_pt) * lpr_Basis(1,2, ai_pt) )

        apr_rhs(2) = lp_bbox % opr_dBasis (4,ai_b-1,ai_pt) &
        - ( ao_ASS % opo_info_gr (ai_grids_id) % opr_points (4, 1, ai_pt) * lpr_Basis(1,1, ai_pt)  &
        +   ao_ASS % opo_info_gr (ai_grids_id) % opr_points (4, 2, ai_pt) * lpr_Basis(1,2, ai_pt) )

        apr_rhs(3) = lp_bbox % opr_dBasis (5,ai_b-1,ai_pt) &
        - ( ao_ASS % opo_info_gr (ai_grids_id) % opr_points (5, 1, ai_pt) * lpr_Basis(1,1, ai_pt)  &
        +   ao_ASS % opo_info_gr (ai_grids_id) % opr_points (5, 2, ai_pt) * lpr_Basis(1,2, ai_pt) )

        !            print *, 'lpr_rhs = ', lpr_rhs

        !...

        end subroutine compute_rhs
end subroutine af_Basis_2D

end module hessian_field_module



