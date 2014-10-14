!     
! File:   stiffness.F90
! Author: root
!
! Created on January 27, 2012, 4:55 PM
!

module stiffness_module
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

    public :: build_stiffness_Local

#ifdef _DEBUG
    integer, parameter, private  :: mi_dtllevel_base = 0
#else
    integer, parameter, private  :: mi_dtllevel_base = 2
#endif
    !> \todo to remove from here, and use an adequate variable
    integer, parameter, private :: mi_MAXDIM = 3

contains

    !----------------------------------------------------------------------------------------------
    subroutine build_stiffness_Local(ao_ASS, ao_FEM, ai_grids_id, ai_id, ai_elt, ai_ref, apr_Operator_elt)
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
        integer :: ai_ref
        real(wp), dimension(0:,:,:) :: apr_Operator_elt
        ! LOCAL VARIABLES
        integer :: li_npts
        integer :: li_b, li_bprime
        integer :: li_i, li_iprime
        integer :: li_A, li_Aprime
        integer :: li_sp, li_spprime
        integer :: li_bloc, li_bprimeloc
        integer :: li_ndof, li_ndofprime
        integer :: li_map
        integer :: li_d
        integer :: li_dprime
        integer :: li_dim
        integer :: li_iparam
        INTEGER :: li_loc_id
        INTEGER :: li_nparam
        INTEGER :: li_nen
        INTEGER :: li_nenprime         
        real(wp) :: lr_contribution
        real(wp), dimension (ao_ASS % oi_maxndof, ao_ASS % oi_maxnpts, mi_MAXDIM) :: lpr_gradB, lpr_gradBprime
        real(wp), dimension(ao_ASS % oi_maxnpts) :: lpr_jacobian
        real(wp), DIMENSION(ao_ASS % oi_maxnpts) :: lpr_contribution
        real(wp), dimension(ao_ASS % oi_maxnpts) :: lpr_w
        real(wp), dimension(ao_ASS % oi_maxnpts, NPARAM_MATRIX) :: lpr_m
        REAL(WP), dimension(ao_FEM % oi_maxnen, ao_FEM % oi_maxnen) :: lpr_Matrix        
        TYPE(CONNECTIVITY), pointer :: lp_con
        TYPE(CONNECTIVITY), pointer :: lp_conprime
        TYPE(GRID_DATA), pointer :: lp_grid
        TYPE(PHYSICAL_BASIS), pointer :: lp_pBasis
        TYPE(PHYSICAL_BASIS), pointer :: lp_pBasisprime
#ifdef _PERF        
        REAL(8) time_begin, time_end
#endif

#ifdef _DEBUG
        call printlog("build_stiffness_Local : Begin", ai_dtllevel = mi_dtllevel_base + 1)
#endif
!print *, 'ai_elt=', ai_elt

        lpr_gradB = 0.0_wp
        lpr_gradBprime = 0.0_wp

        li_dim = ao_FEM % opi_dim(ai_grids_id)
        
        !\todo prblm du id => id+1
        lp_grid         => ao_FEM % opo_grids ( ai_grids_id ) % opo_grid ( ai_id-1 )
        li_npts         = lp_grid % opo_elts ( ai_elt ) % oi_npts

        li_sp           = ao_FEM % opi_InfoOperator (ai_ref, INFOOPERATOR_SPACE_1)
        li_spprime      = ao_FEM % opi_InfoOperator (ai_ref, INFOOPERATOR_SPACE_2)

        li_ndof         = ao_FEM % opi_InfoSpace(li_sp     , INFOSPACE_NDOF )
        li_ndofprime    = ao_FEM % opi_InfoSpace(li_spprime, INFOSPACE_NDOF )

        lp_con          => ao_FEM % opo_spaces (li_sp) % oo_con
        lp_conprime     => ao_FEM % opo_spaces (li_spprime) % oo_con

        lp_pBasis       => ao_ASS % opo_pBasis (li_sp)  
        lp_pBasisprime  => ao_ASS % opo_pBasis (li_spprime) 

        li_loc_id       = ao_FEM % opi_InfoOperator(ai_ref, INFOOPERATOR_LOCID)
        li_nparam       = ao_FEM % opi_InfoOperator(ai_ref, INFOOPERATOR_NPARAM)

        lpr_jacobian    = ao_ASS % opo_info_gr (ai_grids_id) % opr_jacobians
        lpr_jacobian    = dabs(lpr_jacobian)

        li_nen          = lp_con % opi_nen(ai_id)
        li_nenprime     = lp_conprime % opi_nen(ai_id)

        lpr_m (1:li_npts, 1:li_nparam) = TRANSPOSE(lp_grid % opo_elts (ai_elt) % opr_values_m(li_loc_id, 1:li_nparam, 1:li_npts))
        lpr_w (1 : li_npts)            = lp_grid % opo_elts ( ai_elt ) % opr_w (1:li_npts)

        lpr_Matrix = 0.0_wp
#ifdef _PERF        
        CALL CPU_TIME ( time_begin )
#endif
        do li_b = 1, li_nen

            lpr_gradB (1, 1:li_npts, 1:li_dim) = TRANSPOSE(lp_pBasis % opr_gradB (1:li_dim, li_b, 1:li_npts))            

            do li_bprime = 1, li_nenprime

                lpr_gradBprime (1, 1:li_npts, 1:li_dim) =  &
                TRANSPOSE(lp_pBasisprime % opr_gradB (1:li_dim, li_bprime, 1:li_npts))

                !> \todo only work if ndof=1 (scalar unknowns)
                lpr_contribution = 0.0_wp
                li_iparam = 1
                DO li_dprime = 1, li_dim
                    DO li_d = 1, li_dim
                        ! multiply by the dot product of the two basis
                        lpr_contribution (1:li_npts) = lpr_contribution (1:li_npts) &
                        + lpr_gradB (1, 1:li_npts, li_d) * lpr_gradBprime (1, 1:li_npts, li_dprime) &
                        * lpr_m(1:li_npts, li_iparam)
                        li_iparam = li_iparam + 1
                    END DO
                END DO

                ! multiply by input data
                lpr_contribution (1:li_npts) = lpr_contribution (1:li_npts) &              
                ! multiply by quadratures weights
                * lpr_w (1:li_npts) &
                ! multiply by jacobians
                * lpr_jacobian (1:li_npts)

                lr_contribution = SUM ( lpr_contribution (1:li_npts) )

                lpr_Matrix(li_b, li_bprime) = lpr_Matrix(li_b, li_bprime) &
                + lr_contribution

            end do

        end do

#ifdef _PERF        
        CALL CPU_TIME ( time_end )
        PRINT *, 'Time of operation was ', &
        time_end - time_begin, ' seconds'
#endif

        apr_Operator_elt(ai_ref, 1:li_nen, 1:li_nenprime) = lpr_Matrix (1:li_nen, 1:li_nenprime)

#ifdef _DEBUG        
        call printlog("build_stiffness_Local : End", ai_dtllevel = mi_dtllevel_base + 1)
#endif
    end subroutine build_stiffness_Local

end module stiffness_module

