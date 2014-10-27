!     
! File:   second_deriv.F90
! Author: root
!
! Created on June 4, 2012, 9:46 AM
!

module second_deriv_module
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

    public :: build_second_deriv_Local

#ifdef _DEBUG
    integer, parameter, private  :: mi_dtllevel_base = 0
#else
    integer, parameter, private  :: mi_dtllevel_base = 2
#endif
    !> \todo to remove from here, and use an adequate variable
    integer, parameter, private :: mi_MAXDIM = 3
    integer, parameter, private :: mi_NSIZE  = 3 ! 1 => 1D, 3 => 2D, 6 => 3D 
    integer, parameter, private :: mi_NPARAM = 9 ! 1 => 1D, 9 => 2D, 36 => 3D 

contains

    !----------------------------------------------------------------------------------------------
    subroutine build_second_deriv_Local(ao_ASS, ao_FEM, ai_grids_id, ai_id, ai_elt, ai_ref, apr_Operator_elt)
        implicit none
        TYPE(ASSEMBLY) :: ao_ASS
        TYPE(FEM) :: ao_FEM
        !> grids id
        INTEGER :: ai_grids_id
        !> patch id
        INTEGER :: ai_id
        !> current element
        INTEGER :: ai_elt
        !> current matrix-id
        INTEGER :: ai_ref
        REAL(wp), DIMENSION(0:,:,:) :: apr_Operator_elt
        ! LOCAL VARIABLES
        INTEGER :: li_npts
        INTEGER :: li_pt
        INTEGER :: li_b, li_bprime
        INTEGER :: li_i, li_iprime
        INTEGER :: li_A, li_Aprime
        INTEGER :: li_d, li_dprime
        INTEGER :: li_sp, li_spprime
        INTEGER :: li_bloc, li_bprimeloc
        INTEGER :: li_ndof, li_ndofprime
        INTEGER :: li_map
        INTEGER :: li_dim
        INTEGER :: li_iparam
        INTEGER :: li_loc_id
        INTEGER :: li_nparam
        INTEGER :: li_nen
        INTEGER :: li_nenprime          
        INTEGER :: li_nsize
        REAL(wp) :: lr_contribution
        real(wp), dimension (ao_ASS % oi_maxndof, ao_ASS % oi_maxnpts, mi_NSIZE) :: lpr_HessianB, lpr_HessianBprime 
        REAL(wp), DIMENSION (ao_ASS % oi_maxnpts) :: lpr_jacobian
        REAL(wp), DIMENSION (ao_ASS % oi_maxnpts) :: lpr_contribution
        real(wp), dimension(ao_ASS % oi_maxnpts) :: lpr_w
        real(wp), dimension(ao_ASS % oi_maxnpts, NPARAM_MATRIX) :: lpr_m
        REAL(WP), dimension(ao_FEM % oi_maxnen, ao_FEM % oi_maxnen) :: lpr_Matrix        
        TYPE(CONNECTIVITY), POINTER :: lp_con
        TYPE(CONNECTIVITY), POINTER :: lp_conprime
        TYPE(GRID_DATA), POINTER :: lp_grid
        TYPE(PHYSICAL_BASIS), pointer :: lp_pBasis
        TYPE(PHYSICAL_BASIS), pointer :: lp_pBasisprime
#ifdef _PERF        
        REAL(8) time_begin, time_end
#endif

#ifdef _DEBUG
        call printlog("build_second_deriv_Local : Begin", ai_dtllevel = mi_dtllevel_base + 1)
#endif

        lpr_HessianB = 0.0_wp
        lpr_HessianBprime = 0.0_wp

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

        if (li_dim==1) then
                li_nsize = 1
        end if
        if (li_dim==2) then
                li_nsize = 3 
        end if
        if (li_dim==3) then
                li_nsize = 6
        end if

        lpr_m (1:li_npts, 1:li_nparam) = TRANSPOSE(lp_grid % opo_elts (ai_elt) % opr_values_m(li_loc_id, 1:li_nparam, 1:li_npts))
        lpr_w (1 : li_npts)            = lp_grid % opo_elts ( ai_elt ) % opr_w (1:li_npts)

        lpr_Matrix = 0.0_wp

#ifdef _PERF        
        CALL CPU_TIME ( time_begin )
#endif  
        do li_b = 1, li_nen

            lpr_HessianB (1, 1:li_npts, 1:li_nsize) = TRANSPOSE(lp_pBasis % opr_HessianB (1:li_nsize, li_b, 1:li_npts))

            do li_bprime = 1, li_nenprime

                lpr_HessianBprime (1, 1:li_npts, 1:li_nsize) = &
                TRANSPOSE(lp_pBasisprime % opr_HessianB (1:li_nsize, li_bprime, 1:li_npts))

                lpr_contribution = 0.0_wp
                li_iparam = 1
                DO li_dprime = 1, li_nsize
!                PRINT *, '>>>>>>> dprime ', li_dprime
                    DO li_d = 1, li_nsize
!                        PRINT *, '>> d, iprime ', li_d, li_iparam
                        ! multiply by the dot product of the two basis
                        lpr_contribution (1:li_npts) = lpr_contribution (1:li_npts) &
                        + lpr_HessianB (1, 1:li_npts, li_d) &
                        * lpr_HessianBprime (1, 1:li_npts, li_dprime) &
                        * lpr_m(1:li_npts, li_iparam)
!                        print *, '------'
!                        print *, 'B  ', lpr_HessianB (1, 1:li_npts, li_d)    
!                        print *, 'BB ', lpr_HessianBprime (1, 1:li_npts, li_dprime) 
!                        print *, 'm  ', lpr_m(1:li_npts, li_iparam)
!                        print *, '------'
                        li_iparam = li_iparam + 1

                    END DO
                END DO                
!                STOP

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
        call printlog("build_second_deriv_Local : End", ai_dtllevel = mi_dtllevel_base + 1)
#endif
    end subroutine build_second_deriv_Local
    !---------------------------------------------------------------

end module second_deriv_module


