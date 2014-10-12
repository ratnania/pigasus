!     
! File:   rotational_scalar_scalar.F90
! Author: root
!
! Created on April 8, 2012, 5:14 PM
!

module rotational_scalar_module
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

    public :: build_rotational_scalar_Local

#ifdef _DEBUG
    integer, parameter, private  :: mi_dtllevel_base = 0
#else
    integer, parameter, private  :: mi_dtllevel_base = 2
#endif
    !> \todo to remove from here, and use an adequate variable
    integer, parameter, private :: mi_MAXDIM = 3
contains

    !----------------------------------------------------------------------------------------------
    subroutine build_rotational_scalar_Local(ao_ASS, ao_FEM, ai_grids_id, ai_id, ai_elt, ai_matrix, apr_Matrix_elt)
        implicit none
        TYPE(ASSEMBLY) :: ao_ASS
        type(FEM) :: ao_FEM
        !> grids id
        integer :: ai_grids_id
        !> patch id
        integer :: ai_id
        !> current element
        integer :: ai_elt
        !> current matrix-id
        integer :: ai_matrix
        real(wp), dimension(0:,:,:) :: apr_Matrix_elt
        ! LOCAL VARIABLES
        integer :: li_npts
        integer :: li_b, li_bprime
        integer :: li_i, li_iprime
        integer :: li_A, li_Aprime
        integer :: li_sp, li_spprime
        integer :: li_bloc, li_bprimeloc
        integer :: li_ndof, li_ndofprime
        INTEGER :: li_cur_sp, li_cur_spprime
        INTEGER :: li_cur_b, li_cur_bprime
        INTEGER :: li_cur_nen, li_cur_nenprime
        integer :: li_cur_sp_id, li_cur_spprime_id
        integer :: li_map
        INTEGER :: li_loc_id
        INTEGER :: li_nparam
        INTEGER :: li_dim
        INTEGER :: li_d
        INTEGER :: li_conelt
        INTEGER :: li_conprimeelt
        real(wp) :: lr_contribution
        real(wp) :: lr_dotproduct
        real(wp), dimension(ao_ASS % oi_maxndof,ao_ASS % oi_maxnpts) :: lpr_Basis
        real(wp), dimension (ao_ASS % oi_maxndof, mi_MAXDIM, ao_ASS % oi_maxnpts) :: lpr_Basisprime
        real(wp), dimension(ao_ASS % oi_maxnpts) :: lpr_jacobian
        real(wp), DIMENSION(ao_ASS % oi_maxnpts) :: lpr_contribution
        real(wp), dimension(ao_ASS % oi_maxnpts) :: lpr_w
        real(wp), dimension(ao_ASS % oi_maxndof, ao_ASS % oi_maxnpts) :: lpr_vprime ! the unit vector
        TYPE(CONNECTIVITY), pointer :: lp_con
        TYPE(CONNECTIVITY), pointer :: lp_conprime
        TYPE(CONNECTIVITY), pointer :: lp_cur_con
        TYPE(CONNECTIVITY), pointer :: lp_cur_conprime
        TYPE(GRID_DATA), pointer :: lp_grid

#ifdef _DEBUG
        call printlog("build_rotational_scalar_Local : Begin", ai_dtllevel = mi_dtllevel_base + 1)
#endif
!print *, 'ai_elt=', ai_elt
        print *, 'build_rotational_scalar_Local: Not yet implemented'
        STOP

        lpr_Basis = 0.0_wp
        lpr_Basisprime = 0.0_wp

        li_dim = ao_FEM % opi_dim(ai_grids_id)

        !\todo prblm du id => id+1
        lp_grid => ao_FEM % opo_grids ( ai_grids_id ) % opo_grid ( ai_id-1 )
        li_npts = lp_grid % opo_elts ( ai_elt ) % oi_npts

        li_sp       = ao_FEM % opi_InfoMatrix (ai_matrix, INFOMATRIX_SPACE_1)
        li_spprime  = ao_FEM % opi_InfoMatrix (ai_matrix, INFOMATRIX_SPACE_2)

        li_ndof         = ao_FEM % opi_InfoSpace(li_sp     , INFOSPACE_NDOF )
        li_ndofprime    = ao_FEM % opi_InfoSpace(li_spprime, INFOSPACE_NDOF )

!        print*,'current spaces are ', li_sp, li_spprime
        lp_con          => ao_FEM % opo_spaces (li_sp) % oo_con
        lp_conprime     => ao_FEM % opo_spaces (li_spprime) % oo_con

        li_loc_id = ao_FEM % opi_InfoMatrix(ai_matrix, INFOMATRIX_LOCID)
        li_nparam = ao_FEM % opi_InfoMatrix(ai_matrix, INFOMATRIX_NPARAM)
!        lpr_b (1 : li_nparam, 1 : li_npts) =    &
!        lp_grid % opo_elts (ai_elt) % opr_values_m(li_loc_id, 1 : li_nparam, 1 : li_npts)

        ! if we don't use an exterior mapping, then we must use the transformation defined
        ! by the current geometry
!        if ( ao_FEM % opi_InfoSpace ( li_sp , INFOSPACE_EXTMAPPING ) == 1 ) then
!            li_map = ao_FEM % opi_InfoSpace ( li_sp , INFOSPACE_MAPPING )
!            lpr_jacobian = ao_ASS % opo_info_mp (li_map) % opr_jacobians
!        else
!            lpr_jacobian = ao_ASS % opo_info_sp (li_sp) % opr_jacobians
!        end if
!        lpr_jacobian = dabs(lpr_jacobian)
!        !> \todo a supprimer apres
        lpr_jacobian = 1.0_wp

        li_conelt = lp_con % opi_real_elts (ai_id-1,ai_elt)
        li_conprimeelt = lp_conprime % opi_real_elts (ai_id-1,ai_elt)

        do li_b = 1, lp_con % opi_nen (ai_id, li_conelt)

            li_A = lp_con % opi_LM ( ai_id, li_b, li_conelt )

            if (li_A == 0) then
                cycle
            end if

!            print *, 'li_b=', li_b, 'li_cur_b=', li_cur_b

            CALL af_Basis(ao_ASS, ao_FEM, li_sp, li_ndof, ai_grids_id, ai_id, ai_elt, li_b, lpr_Basis )

!print *, 'lpr_Basis=', lpr_Basis

            li_cur_spprime = 1
            li_cur_spprime_id = ao_FEM % opi_InfoComposedSpace(li_spprime, li_cur_spprime)
            li_cur_nenprime = ao_FEM % opo_spaces (li_cur_spprime_id) % oo_con % opi_nen (ai_id, ai_elt)
            li_cur_bprime = 1
            lpr_vprime (:,:) = 0.0_wp ; lpr_vprime (li_cur_spprime,1 : li_npts) = 1.0_wp
            do li_bprime = 1, lp_conprime % opi_nen (ai_id, li_conprimeelt)

                li_Aprime = lp_conprime % opi_LM ( ai_id, li_bprime, li_conprimeelt )

                if (li_Aprime == 0) then
                    if ( (li_cur_bprime == li_cur_nenprime) .AND. (li_bprime /= lp_conprime % opi_nen (ai_id, ai_elt)) ) then
                        li_cur_bprime = 0
                        !> \todo this works only for V = V1xV2
                        li_cur_spprime = li_cur_spprime + 1
                        li_cur_spprime_id = ao_FEM % opi_InfoComposedSpace(li_spprime, li_cur_spprime)
                        li_cur_nenprime = ao_FEM % opo_spaces (li_cur_spprime_id) % oo_con % opi_nen (ai_id, ai_elt)
                        lpr_vprime (:,:) = 0.0_wp ; lpr_vprime (li_cur_spprime,1 : li_npts) = 1.0_wp
                    end if
                    li_cur_bprime = li_cur_bprime + 1
                    
                    cycle
                end if

!print *, 'li_bprime=', li_bprime , 'li_cur_bprime=', li_cur_bprime
!print *, 'li_cur_bprime=', li_cur_bprime

                CALL af_Basisprime(ao_ASS, ao_FEM, li_spprime, li_ndofprime, ai_grids_id, ai_id, ai_elt &
                , li_bprime, li_cur_spprime_id, li_cur_bprime, li_cur_nenprime, lpr_Basisprime)

!print *, 'lpr_Basisprime=', lpr_Basisprime

                lr_dotproduct = 0.0_wp
!                DO li_d = 1, li_dim
!                    lr_dotproduct = lr_dotproduct + DOT_PRODUCT(lpr_Basis (1,li_d,1:li_npts),lpr_vprime(li_d,1:li_npts))
!                END DO
                !> \todo only work if ndof=1 (scalar unknowns)
                ! multiply by input data
                lpr_contribution (1:li_npts) = lp_grid % opo_elts (ai_elt) % opr_values_m(li_loc_id, 1, 1 : li_npts)   &                
                ! multiply by the right basis
!                * lpr_Basisprime (1,1:li_npts)   &
                ! multiply by the dot product for each direction
                * lr_dotproduct  &
                ! multiply by quadratures weights
                * lp_grid % opo_elts ( ai_elt ) % opr_w (1:li_npts) &
                ! multiply by jacobians
                * lpr_jacobian (1:li_npts)

                lr_contribution = SUM ( lpr_contribution (1:li_npts) )
!                print *, 'lr_contribution=', lr_contribution

                apr_Matrix_elt(ai_matrix, li_b, li_bprime) = apr_Matrix_elt(ai_matrix, li_b, li_bprime) &
                + lr_contribution


                if ( (li_cur_bprime == li_cur_nenprime) .AND. (li_bprime /= lp_conprime % opi_nen (ai_id, ai_elt)) ) then
                    li_cur_bprime = 0
                    !> \todo this works only for V = V1xV2
                    li_cur_spprime = li_cur_spprime + 1
                    li_cur_spprime_id = ao_FEM % opi_InfoComposedSpace(li_spprime, li_cur_spprime)
                    li_cur_nenprime = ao_FEM % opo_spaces (li_cur_spprime_id) % oo_con % opi_nen (ai_id, ai_elt)
                    lpr_vprime (:,:) = 0.0_wp ; lpr_vprime (li_cur_spprime,1 : li_npts) = 1.0_wp
                end if
                li_cur_bprime = li_cur_bprime + 1
                
            end do

        end do
!print*,'apr_Matrix_elt',apr_Matrix_elt
!STOP

#ifdef _DEBUG
        call printlog("build_rotational_scalar_Local : End", ai_dtllevel = mi_dtllevel_base + 1)
#endif
    end subroutine build_rotational_scalar_Local
    !---------------------------------------------------------------
    subroutine af_Basis(ao_ASS, ao_FEM, ai_space, ai_ndof, ai_grids_id, ai_id, ai_elt, ai_b, apr_Basis)
        implicit none
        TYPE(ASSEMBLY) :: ao_ASS
        TYPE(FEM) :: ao_FEM
        INTEGER :: ai_space
        INTEGER :: ai_ndof
        INTEGER :: ai_grids_id
        INTEGER :: ai_id
        INTEGER :: ai_elt
        INTEGER :: ai_b
        real(wp), DIMENSION(:,:) :: apr_Basis
        ! LOCAL
        TYPE(BBOX), POINTER :: lp_bbox
        TYPE(GRID_DATA), pointer :: lp_grid
        TYPE(GEOMETRIES), pointer :: lp_geos
        INTEGER :: li_d
        INTEGER :: li_npts

#ifdef _DEBUG
        call printlog("af_Basis : Begin", ai_dtllevel = mi_dtllevel_base + 3)
#endif

        !\todo prblm du id => id+1
        lp_grid => ao_FEM % opo_grids (ai_grids_id) % opo_grid ( ai_id-1 )

        li_npts = lp_grid % opo_elts ( ai_elt ) % oi_npts

        lp_bbox => ao_ASS % opo_bbox_sp ( ai_space )
        lp_geos => ao_FEM % opo_spaces ( ai_space ) % oo_mapping

        !> \todo traiter le cas ou ndof > 1
        ! apr_Basis (1:li_ndof,1:li_npts)
        !> \todo on doit enlever 1 a cause des indices
        !> bbox % dBasis commence a 0
        apr_Basis (1,:) = lp_bbox % opr_dBasis (0,ai_b-1,:)

#ifdef _DEBUG
        call printlog("af_Basis : End", ai_dtllevel = mi_dtllevel_base + 3)
#endif
    end subroutine af_Basis
    !---------------------------------------------------------------
        subroutine af_Basisprime(ao_ASS, ao_FEM, ai_space, ai_ndof, ai_grids_id, ai_id, ai_elt, ai_b &
        , ai_cur_sp, ai_cur_b, ai_cur_nen, apr_Basis)
        implicit none
        TYPE(ASSEMBLY) :: ao_ASS
        TYPE(FEM) :: ao_FEM
        INTEGER :: ai_space
        INTEGER :: ai_ndof
        INTEGER :: ai_grids_id
        INTEGER :: ai_id
        INTEGER :: ai_elt
        INTEGER :: ai_b
        INTEGER :: ai_cur_sp
        INTEGER :: ai_cur_b
        INTEGER :: ai_cur_nen
        real(wp), DIMENSION(:,:,:) :: apr_Basis
        ! LOCAL
        TYPE(BBOX), POINTER :: lp_bbox
        TYPE(GRID_DATA), pointer :: lp_grid
        TYPE(GEOMETRIES), pointer :: lp_geos
        TYPE(METRIC_INFO), pointer :: lp_info
        real(wp), dimension(:,:,:), pointer :: lpr_invJacobian
        INTEGER :: li_d
        INTEGER :: li_npts
        INTEGER :: li_dim
        INTEGER :: li_pt

#ifdef _DEBUG
        call printlog("af_Basisprime : Begin", ai_dtllevel = mi_dtllevel_base + 3)
#endif

        li_dim = ao_FEM % opi_dim(ai_grids_id)

        !\todo prblm du id => id+1
        lp_grid => ao_FEM % opo_grids (ai_grids_id) % opo_grid ( ai_id-1 )
        lp_info => ao_ASS % opo_info_sp (ai_id-1)

        li_npts = lp_grid % opo_elts ( ai_elt ) % oi_npts

        lp_bbox => ao_ASS % opo_bbox_sp ( ai_cur_sp )
        lp_geos => ao_FEM % opo_spaces ( ai_space ) % oo_mapping

        lpr_invJacobian => ao_ASS % opo_info_gr (ai_grids_id) % opr_invJacobian

        !> \todo traiter le cas ou ndof > 1
        ! apr_Basis (1:li_ndof,1:li_npts)
        !> \todo on doit enlever 1 a cause des indices
        !> bbox % dBasis commence a 0
!        apr_Basis (1,:) = lp_bbox % opr_dBasis (1,ai_b-1,:)

        !> we multiply by the inverse of the jacobian  matrix
        do li_pt=1, li_npts
            apr_Basis (1,1:li_dim,li_pt) = MATMUL(lpr_invJacobian(1:li_dim, 1:li_dim, li_pt) &
            , lp_bbox % opr_dBasis (1:li_dim,ai_b-1,li_pt))
        end do

#ifdef _DEBUG
        call printlog("af_Basisprime : End", ai_dtllevel = mi_dtllevel_base + 3)
#endif
        end subroutine af_Basisprime

end module rotational_scalar_module





