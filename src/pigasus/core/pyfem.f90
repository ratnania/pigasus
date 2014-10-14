!     
! File:   pyfem.F90
! Author: root
!
! Created on January 11, 2012, 2:34 PM
!

!
! File:   pyfem.F90
! Author: ratnani
!
! Created on december 8, 2011, 9:28 AM
!

module pyfem
    use used_precision
    use tracelog_module
    use grids, only : set_values_field, set_values_matrix, set_values_norm, get_weights &
                         ,  create_from_tensor_product_1d,  create_from_tensor_product_2d,  create_from_tensor_product_3d &
                         , set_unidirection_tensor_2d, print_grid
!    use assl_module
    use fem_def
    use grids_module
    use fem_module
    use assembly_def
    use assembly_module
    use fields_tools_module
!#ifdef _FIGA
!    use fast_iga_def
!    use fast_iga_module
!#endif
    implicit none

    type(FEM)        , save, private :: mo_fem
    type(ASSEMBLY)   , save, private :: mo_ass
!#ifdef _FIGA
!    type(FAST_IGA)   , save, private :: mo_figa
!#endif

!#ifdef _DEBUG
!    integer, parameter, private  :: mi_dtllevel_base = 0
!#else
!    integer, parameter, private  :: mi_dtllevel_base = 2
!#endif
    integer, parameter, private  :: mi_dtllevel_base = 2
    
contains
!#ifdef _FIGA
!include "pyfem_figa.F90"
!#endif
    !---------------------------------------------------------------
    subroutine init(ai_stdoutput, ai_detail)
        implicit none
        integer :: ai_stdoutput
        integer :: ai_detail
        ! LOCAL
        logical :: ll_stdoutput

        ll_stdoutput = .TRUE.
        if (ai_stdoutput == 1) then
            ll_stdoutput = .FALSE.
        end if

        call opentracelog(al_stdoutput = ll_stdoutput, ai_dtllevel = ai_detail)

    end subroutine init
    !---------------------------------------------------------------
    subroutine pyfem_get_maxnderiv_field (ai_val)
        implicit none
        INTEGER, INTENT(OUT) :: ai_val

        ai_val = 1
    end subroutine pyfem_get_maxnderiv_field
    !---------------------------------------------------------------
    subroutine set_nspaces (ai_nspaces)
        implicit none
        INTEGER :: ai_nspaces

        CALL set_nspaces_fem (mo_fem, ai_nspaces)
    end subroutine set_nspaces
    !---------------------------------------------------------------
    subroutine set_nmappings (ai_nmappings)
        implicit none
        INTEGER :: ai_nmappings

        CALL set_nmappings_fem (mo_fem, ai_nmappings)
    end subroutine set_nmappings
    !---------------------------------------------------------------
    subroutine set_nmetrics (ai_val)
        implicit none
        INTEGER :: ai_val

        CALL set_nmetrics_fem (mo_fem, ai_val)
    end subroutine set_nmetrics
    !---------------------------------------------------------------
    subroutine set_nsolvers (ai_n)
        implicit none
        INTEGER :: ai_n

        CALL set_nsolvers_fem (mo_fem, ai_n)
    end subroutine set_nsolvers     
    !---------------------------------------------------------------
    subroutine set_ngraphs (ai_n)
        implicit none
        INTEGER :: ai_n

        CALL set_ngraphs_fem (mo_fem, ai_n)
    end subroutine set_ngraphs    
    !---------------------------------------------------------------
    subroutine set_noperators (ai_n)
        implicit none
        INTEGER :: ai_n

        CALL set_noperators_fem (mo_fem, ai_n)
    end subroutine set_noperators    
    !---------------------------------------------------------------
    subroutine set_nmatrices (ai_nMatrices)
        implicit none
        INTEGER :: ai_nMatrices

        CALL set_nmatrices_fem (mo_fem, ai_nMatrices)
    end subroutine set_nmatrices
    !---------------------------------------------------------------
    subroutine set_nfields (ai_nFields)
        implicit none
        INTEGER :: ai_nFields

        CALL set_nfields_fem (mo_fem, ai_nFields)
    end subroutine set_nfields
    !---------------------------------------------------------------
    subroutine set_nnorms (ai_n)
        implicit none
        INTEGER :: ai_n

        CALL set_nnorms_fem (mo_fem, ai_n)
    end subroutine set_nnorms
    !---------------------------------------------------------------
    subroutine set_ngrids (ai_nGrids)
        implicit none
        INTEGER :: ai_nGrids

        CALL set_ngrids_fem (mo_fem, ai_nGrids)
    end subroutine set_ngrids
    !---------------------------------------------------------------
    subroutine set_maxnpatchs (ai_val)
        implicit none
        INTEGER :: ai_val

        CALL set_maxnpatchs_fem (mo_fem, ai_val)
    end subroutine set_maxnpatchs
    !---------------------------------------------------------------
    subroutine set_maxnparams_operators (ai_val)
        implicit none
        INTEGER :: ai_val

        CALL set_maxnparams_operators_fem (mo_fem, ai_val)
    end subroutine set_maxnparams_operators
    !---------------------------------------------------------------
    subroutine set_maxnaddto (ai_val)
        implicit none
        INTEGER :: ai_val

        CALL set_maxnaddto_fem (mo_fem, ai_val)
    end subroutine set_maxnaddto    
    !---------------------------------------------------------------
    subroutine set_maxnparams_fields (ai_val)
        implicit none
        INTEGER :: ai_val

        CALL set_maxnparams_fields_fem (mo_fem, ai_val)
    end subroutine set_maxnparams_fields
    !---------------------------------------------------------------
    subroutine set_maxnparams_norms (ai_val)
        implicit none
        INTEGER :: ai_val

        CALL set_maxnparams_norms_fem (mo_fem, ai_val)
    end subroutine set_maxnparams_norms
    !---------------------------------------------------------------
    subroutine set_maxder_field_op (ai_val)
        implicit none
        INTEGER :: ai_val

        CALL set_maxder_field_op_fem (mo_fem, ai_val)
    end subroutine set_maxder_field_op
!----------------------------------------------------------------------------------------------
    subroutine set_solver ( ai_id, ai_matrix, ai_residual, ai_solver )
        implicit none
        integer  :: ai_id
        integer  :: ai_matrix
        integer  :: ai_residual
        integer  :: ai_solver

        mo_fem % opi_InfoSolver ( ai_id , : ) = -1

        mo_fem % opi_InfoSolver ( ai_id , INFOSOLVER_MATRIX ) = ai_matrix
        mo_fem % opi_InfoSolver ( ai_id , INFOSOLVER_RESIDUAL  ) = ai_residual
        mo_fem % opi_InfoSolver ( ai_id , INFOSOLVER_SOLVER  ) = ai_solver

    end subroutine set_solver    
!----------------------------------------------------------------------------------------------
    subroutine set_field ( ai_id, ai_ndof, ai_space, ai_size, ai_loc_id, ai_type &
        , ai_operator, ai_operande, ai_parameval, ai_nparam)
        implicit none
        integer  :: ai_id
        integer  :: ai_ndof
        integer  :: ai_size
        integer  :: ai_space
        integer  :: ai_loc_id
        integer  :: ai_type
        integer  :: ai_operator
        integer  :: ai_operande
        integer  :: ai_parameval
        integer  :: ai_nparam

        mo_fem % opi_InfoField ( ai_id , : ) = -1

        mo_fem % opi_InfoField ( ai_id , INFOFIELD_NDOF         ) = ai_ndof
        mo_fem % opi_InfoField ( ai_id , INFOFIELD_SIZE         ) = ai_size
        mo_fem % opi_InfoField ( ai_id , INFOFIELD_SPACE        ) = ai_space
        mo_fem % opi_InfoField ( ai_id , INFOFIELD_LOCID        ) = ai_loc_id
        mo_fem % opi_InfoField ( ai_id , INFOFIELD_TYPE         ) = ai_type
        mo_fem % opi_InfoField ( ai_id , INFOFIELD_OPERATOR     ) = ai_operator
        mo_fem % opi_InfoField ( ai_id , INFOFIELD_OPERANDE     ) = ai_operande
        mo_fem % opi_InfoField ( ai_id , INFOFIELD_PARAMEVAL    ) = ai_parameval
        mo_fem % opi_InfoField ( ai_id , INFOFIELD_NPARAM    ) = ai_nparam

    end subroutine set_field
!----------------------------------------------------------------------------------------------
    subroutine set_norm ( ai_id, ai_field, ai_type, ai_locid )
        implicit none
        integer  :: ai_id
        integer  :: ai_field
        integer  :: ai_type
        integer  :: ai_locid

        mo_fem % opi_InfoNorm ( ai_id , : ) = -1

        mo_fem % opi_InfoNorm ( ai_id , INFONORM_FIELD ) = ai_field
        mo_fem % opi_InfoNorm ( ai_id , INFONORM_TYPE  ) = ai_type
        mo_fem % opi_InfoNorm ( ai_id , INFONORM_LOCID  ) = ai_locid

    end subroutine set_norm
!--------------------------------------------------------------------
    subroutine set_matrix ( ai_id, ai_type, ai_graph_id, ai_ijvassembly )
        implicit none
        integer  :: ai_id
        integer  :: ai_type
        integer  :: ai_graph_id
        integer  :: ai_ijvassembly

        CALL printlog("set_matrix : Start", ai_dtllevel = 1)

        mo_fem % opi_InfoMatrix ( ai_id , : ) = -1

        mo_fem % opi_InfoMatrix ( ai_id , INFOMATRIX_TYPE ) = ai_type
        mo_fem % opi_InfoMatrix ( ai_id , INFOMATRIX_GRAPH  ) = ai_graph_id
        mo_fem % opi_InfoMatrix ( ai_id , INFOMATRIX_IJVASSEMBLY ) = ai_ijvassembly 

        CALL printlog("set_matrix : End", ai_dtllevel = 1)

    end subroutine set_matrix
!--------------------------------------------------------------------
    subroutine set_operator( ai_id, ai_type, ai_space_1, ai_space_2, ai_loc_id   &
        , ai_nparam, ai_parameval, ai_transpose)
        implicit none
        integer  :: ai_id
        integer  :: ai_type
        integer  :: ai_space_1
        integer  :: ai_space_2
        integer  :: ai_loc_id
        integer  :: ai_nparam
        integer  :: ai_parameval
        integer  :: ai_transpose

        CALL printlog("set_operator: Start", ai_dtllevel = 1)

        mo_fem % opi_InfoOperator( ai_id , : ) = -1

        mo_fem % opi_InfoOperator ( ai_id , INFOOPERATOR_TYPE ) = ai_type
        mo_fem % opi_InfoOperator ( ai_id , INFOOPERATOR_SPACE_1  ) = ai_space_1
        mo_fem % opi_InfoOperator ( ai_id , INFOOPERATOR_SPACE_2  ) = ai_space_2
        mo_fem % opi_InfoOperator ( ai_id , INFOOPERATOR_LOCID  ) = ai_loc_id
        mo_fem % opi_InfoOperator ( ai_id , INFOOPERATOR_NPARAM  ) = ai_nparam
        mo_fem % opi_InfoOperator ( ai_id , INFOOPERATOR_PARAMEVAL  ) = ai_parameval
        mo_fem % opi_InfoOperator ( ai_id , INFOOPERATOR_TRANSPOSE  ) = ai_transpose
        
        CALL printlog("set_operator : End", ai_dtllevel = 1)

    end subroutine set_operator    
!----------------------------------------------------------------------------------------------
    subroutine set_space ( ai_id, ai_extmapping, ai_mapping, ai_tensor, ai_grids_id &
        , ai_storeddata, ai_ndof, ai_size, ai_maxnen, ai_composed_space, ai_composed_nspaces)
        implicit none
        integer  :: ai_id
        integer  :: ai_extmapping
        integer  :: ai_mapping
        integer  :: ai_tensor
        integer  :: ai_grids_id
        integer  :: ai_storeddata
        integer  :: ai_ndof
        integer  :: ai_size
        integer  :: ai_maxnen
        integer  :: ai_composed_space
        integer  :: ai_composed_nspaces

        mo_fem % opi_InfoSpace ( ai_id , : ) = -1

        mo_fem % opi_InfoSpace ( ai_id , INFOSPACE_EXTMAPPING   ) = ai_extmapping
        mo_fem % opi_InfoSpace ( ai_id , INFOSPACE_MAPPING      ) = ai_mapping
        mo_fem % opi_InfoSpace ( ai_id , INFOSPACE_TENSOR       ) = ai_tensor
        mo_fem % opi_InfoSpace ( ai_id , INFOSPACE_GRIDS        ) = ai_grids_id
        mo_fem % opi_InfoSpace ( ai_id , INFOSPACE_STOREDDATA   ) = ai_storeddata
        mo_fem % opi_InfoSpace ( ai_id , INFOSPACE_NDOF         ) = ai_ndof
        mo_fem % opi_InfoSpace ( ai_id , INFOSPACE_SIZE         ) = ai_size
        mo_fem % opi_InfoSpace ( ai_id , INFOSPACE_MAXNEN       ) = ai_maxnen
        mo_fem % opi_InfoSpace ( ai_id , INFOSPACE_COMPOSED     ) = ai_composed_space
        mo_fem % opi_InfoSpace ( ai_id , INFOSPACE_COMPOSED_NSP ) = ai_composed_nspaces


    end subroutine set_space
!----------------------------------------------------------------------------------------------
    subroutine set_graph ( ai_id, ai_space_1, ai_space_2 )
        implicit none
        integer  :: ai_id
        integer  :: ai_space_1
        integer  :: ai_space_2

        mo_fem % opi_InfoGraph ( ai_id , : ) = -1

        mo_fem % opi_InfoGraph ( ai_id , INFOGRAPH_SPACE_1 ) = ai_space_1
        mo_fem % opi_InfoGraph ( ai_id , INFOGRAPH_SPACE_2 ) = ai_space_2

    end subroutine set_graph    
!----------------------------------------------------------------------------------------------
    subroutine set_mapping ( ai_id, ai_tensor, ai_storeddata, ai_space )
        implicit none
        integer  :: ai_id
        integer  :: ai_tensor
        integer  :: ai_storeddata
        integer  :: ai_space

        CALL printlog("set_mapping : Start", ai_dtllevel = 1)

        mo_fem % opi_InfoMapping ( ai_id , : ) = -1

        mo_fem % opi_InfoMapping ( ai_id , INFOMAPPING_TENSOR     ) = ai_tensor
        mo_fem % opi_InfoMapping ( ai_id , INFOMAPPING_STOREDDATA ) = ai_storeddata
        mo_fem % opi_InfoMapping ( ai_id , INFOMAPPING_SPACE      ) = ai_space

        CALL printlog("set_mapping : End", ai_dtllevel = 1)

    end subroutine set_mapping
!----------------------------------------------------------------------------------------------
    subroutine set_patch ( ai_grids_id, ai_id, ai_nel, ai_maxnpts, ai_tensor, ai_dirmaxnpts )
        implicit none
        integer  :: ai_grids_id
        integer  :: ai_id
        integer  :: ai_nel
        integer  :: ai_maxnpts
        integer  :: ai_tensor
        integer  :: ai_dirmaxnpts

        CALL printlog("set_patch : Start", ai_dtllevel = 1)

        mo_fem % opi_InfoPatch ( ai_grids_id, ai_id , : ) = -1

        mo_fem % opi_InfoPatch ( ai_grids_id, ai_id , INFOPATCH_NEL      ) = ai_nel
        mo_fem % opi_InfoPatch ( ai_grids_id, ai_id , INFOPATCH_MAXNPTS  ) = ai_maxnpts
        mo_fem % opi_InfoPatch ( ai_grids_id, ai_id , INFOPATCH_TENSOR   ) = ai_tensor
        mo_fem % opi_InfoPatch ( ai_grids_id, ai_id , INFOPATCH_DIRMAXNPTS  ) = ai_dirmaxnpts

        CALL printlog("set_patch : End", ai_dtllevel = 1)

    end subroutine set_patch
!----------------------------------------------------------------------------------------------
    subroutine set_grids ( ai_id, ai_nPatchs, ai_dim, ai_Rd, ai_dof, ai_nfields &
        , ai_nmatrices, ai_usemetric, ai_metric_id, ai_nnorms )
        implicit none
        integer  :: ai_id
        integer  :: ai_nPatchs
        integer  :: ai_dim
        integer  :: ai_Rd
        integer  :: ai_dof
        integer  :: ai_nfields
        integer  :: ai_nmatrices
        integer  :: ai_usemetric
        integer  :: ai_metric_id
        integer  :: ai_nnorms

        CALL printlog("set_grids : Start", ai_dtllevel = 1)

        mo_fem % opi_InfoGrids ( ai_id , : ) = -1

        mo_fem % opi_InfoGrids ( ai_id , INFOGRIDS_NPATCHS ) = ai_npatchs
        mo_fem % opi_InfoGrids ( ai_id , INFOGRIDS_DIM ) = ai_dim
        mo_fem % opi_InfoGrids ( ai_id , INFOGRIDS_RD ) = ai_Rd
        mo_fem % opi_InfoGrids ( ai_id , INFOGRIDS_DOF ) = ai_dof
        mo_fem % opi_InfoGrids ( ai_id , INFOGRIDS_NFIELDS ) = ai_nfields
        mo_fem % opi_InfoGrids ( ai_id , INFOGRIDS_NMATRICES ) = ai_nmatrices
        mo_fem % opi_InfoGrids ( ai_id , INFOGRIDS_USEMETRIC ) = ai_usemetric
        mo_fem % opi_InfoGrids ( ai_id , INFOGRIDS_METRIC_ID ) = ai_metric_id
        mo_fem % opi_InfoGrids ( ai_id , INFOGRIDS_NNORMS ) = ai_nnorms

        CALL printlog("set_grids : End", ai_dtllevel = 1)

    end subroutine set_grids
!----------------------------------------------------------------------------------------------
    subroutine pyfem_reset_terms_toassembly ( ai_grids_id )
        implicit none
        integer  :: ai_grids_id

        CALL reset_patchs_toassembly(ai_grids_id)
        CALL reset_operators_toassembly()
        CALL reset_matrices_toassembly()
        CALL reset_fields_toassembly()
        CALL reset_norms_toassembly()
        CALL reset_spaces_toassembly()

    end subroutine pyfem_reset_terms_toassembly
!----------------------------------------------------------------------------------------------
    subroutine pyfem_save_terms_toassembly ( )
        implicit none

        CALL save_terms_toassembly(mo_ASS)

    end subroutine pyfem_save_terms_toassembly
!----------------------------------------------------------------------------------------------
    subroutine pyfem_load_terms_toassembly ( )
        implicit none

        CALL load_terms_toassembly(mo_ASS)

    end subroutine pyfem_load_terms_toassembly
!----------------------------------------------------------------------------------------------
    subroutine set_space_toassembly ( ai_id )
        implicit none
        integer  :: ai_id

        mo_fem % opi_Infospace ( ai_id , INFOSPACE_TOASSEMBLY ) = 1

    end subroutine set_space_toassembly    
!----------------------------------------------------------------------------------------------
    subroutine reset_spaces_toassembly ( )
        implicit none

        mo_fem % opi_InfoSpace ( : , INFOSPACE_TOASSEMBLY ) = 0

    end subroutine reset_spaces_toassembly    
!----------------------------------------------------------------------------------------------
    subroutine set_patch_toassembly ( ai_grids_id, ai_id )
        implicit none
        integer  :: ai_grids_id
        integer  :: ai_id

        mo_fem % opi_Infopatch ( ai_grids_id, ai_id , INFOPATCH_TOASSEMBLY ) = 1

    end subroutine set_patch_toassembly
!----------------------------------------------------------------------------------------------
    subroutine reset_patchs_toassembly ( ai_grids_id )
        implicit none
        INTEGER :: ai_grids_id

        mo_fem % opi_Infopatch ( ai_grids_id, : , INFOPATCH_TOASSEMBLY ) = 0

    end subroutine reset_patchs_toassembly
!----------------------------------------------------------------------------------------------
    subroutine pyfem_set_patchs_toassembly ( ai_grids_id )
        implicit none
        INTEGER, intent(in)  :: ai_grids_id

        CALL set_patchs_toassembly(mo_ass, mo_fem, ai_grids_id)

    end subroutine pyfem_set_patchs_toassembly
!----------------------------------------------------------------------------------------------
    subroutine pyfem_set_elts_toassembly ( ai_grids_id, api_values, ai_size )
        implicit none
        INTEGER, intent(in)  :: ai_grids_id
        INTEGER, intent(in)  :: ai_size
        INTEGER , dimension(ai_size), intent(in)  :: api_values

        CALL set_elts_toassembly(mo_ass, mo_fem, ai_grids_id, api_values, ai_size )
        
    end subroutine pyfem_set_elts_toassembly
!----------------------------------------------------------------------------------------------
    subroutine pyfem_set_operator_matrices_toassembly ( ai_operator, api_values, ai_size )
    !> for eahc operator ,we specify the matrices to assembly
        implicit none
        INTEGER, intent(in)  :: ai_operator
        INTEGER, intent(in)  :: ai_size
        INTEGER , dimension(ai_size), intent(in)  :: api_values

        mo_fem % opi_InfoOperator ( ai_operator , INFOOPERATOR_TOASSEMBLY ) = 1
        CALL set_operator_matrices_toassembly(mo_fem, ai_operator, api_values, ai_size )
        
    end subroutine pyfem_set_operator_matrices_toassembly
!----------------------------------------------------------------------------------------------
    subroutine set_matrix_toassembly (ai_id)
        implicit none
        INTEGER :: ai_id

        mo_fem % opi_InfoMatrix ( ai_id , INFOMATRIX_TOASSEMBLY ) = 1 

    end subroutine set_matrix_toassembly
!----------------------------------------------------------------------------------------------
    subroutine reset_matrices_toassembly ( )
        implicit none

        mo_fem % opi_InfoMatrix ( : , INFOMATRIX_TOASSEMBLY ) = 0   

    end subroutine reset_matrices_toassembly    
!----------------------------------------------------------------------------------------------
    subroutine reset_operators_toassembly ( )
        implicit none

        mo_fem % opi_InfoOperator ( : , INFOOPERATOR_TOASSEMBLY ) = 0

    end subroutine reset_operators_toassembly
!----------------------------------------------------------------------------------------------
    subroutine pyfem_reset_matrix ( ai_id )
        implicit none
        integer :: ai_id
        ! LOCAL

        CALL reset_matrix(mo_fem, ai_id)

    end subroutine pyfem_reset_matrix
!---------------------------------------------------------------
    subroutine set_field_toassembly ( ai_id )
        implicit none
        integer  :: ai_id

        mo_fem % opi_InfoField ( ai_id , INFOFIELD_TOASSEMBLY ) = 1

    end subroutine set_field_toassembly
!----------------------------------------------------------------------------------------------
    subroutine reset_fields_toassembly ( )
        implicit none

        mo_fem % opi_InfoField ( : , INFOFIELD_TOASSEMBLY ) = 0

    end subroutine reset_fields_toassembly
!----------------------------------------------------------------------------------------------
    subroutine set_norm_toassembly ( ai_id )
        implicit none
        integer  :: ai_id

        mo_fem % opi_InfoNorm ( ai_id , INFONORM_TOASSEMBLY ) = 1

    end subroutine set_norm_toassembly
!----------------------------------------------------------------------------------------------
    subroutine reset_norms_toassembly ( )
        implicit none

        mo_fem % opi_InfoNorm ( : , INFONORM_TOASSEMBLY ) = 0

    end subroutine reset_norms_toassembly
    !---------------------------------------------------------------
    subroutine SET_SOLVER_MAXITER(ai_id, ai_val)
        implicit none
        INTEGER :: ai_id
        INTEGER :: ai_val

        CALL SET_SOLVER_MAXITER_FEM(mo_fem, ai_id, ai_val)

    end subroutine SET_SOLVER_MAXITER
    !---------------------------------------------------------------
    subroutine SET_SOLVER_RTOL(ai_id, ar_val)
        implicit none
        INTEGER :: ai_id
        REAL(8) :: ar_val

        CALL SET_SOLVER_RTOL_FEM(mo_fem, ai_id, ar_val)

    end subroutine SET_SOLVER_RTOL
    !---------------------------------------------------------------
    subroutine SET_SOLVER_ATOL(ai_id, ar_val)
        implicit none
        INTEGER :: ai_id
        REAL(8) :: ar_val

        CALL SET_SOLVER_ATOL_FEM(mo_fem, ai_id, ar_val)

    end subroutine SET_SOLVER_ATOL
    !---------------------------------------------------------------
    subroutine SET_SOLVER_EPS(ai_id, ar_val)
        implicit none
        INTEGER :: ai_id
        REAL(8) :: ar_val

        CALL SET_SOLVER_EPS_FEM(mo_fem, ai_id, ar_val)

    end subroutine SET_SOLVER_EPS
    !---------------------------------------------------------------
    subroutine GET_SOLVER_NITER(ai_id, ai_val)
        implicit none
        INTEGER :: ai_id
        INTEGER, INTENT(OUT) :: ai_val

        CALL GET_SOLVER_NITER_FEM(mo_fem, ai_id, ai_val)

    end subroutine GET_SOLVER_NITER
    !---------------------------------------------------------------
    subroutine GET_SOLVER_ERR(ai_id, apr_val, ai_size)
        implicit none
        INTEGER, INTENT(IN) :: ai_id
        INTEGER, INTENT(IN) :: ai_size
        REAL(8), DIMENSION(ai_size), INTENT(OUT) :: apr_val

        CALL GET_SOLVER_ERR_FEM(mo_fem, ai_id, apr_val, ai_size)

    end subroutine GET_SOLVER_ERR
    !---------------------------------------------------------------
    subroutine GET_SOLVER_RESNORM(ai_id, ar_val)
        implicit none
        INTEGER, INTENT(IN) :: ai_id
        REAL(8), INTENT(OUT) :: ar_val

        CALL GET_SOLVER_RESNORM_FEM(mo_fem, ai_id, ar_val)

    end subroutine GET_SOLVER_RESNORM    
    !---------------------------------------------------------------
    subroutine create_partone()
        implicit none

        CALL create_fem_partone (mo_fem)
        
    end subroutine create_partone
    !---------------------------------------------------------------
    subroutine create_parttwo()
        implicit none
        ! LOCAL

        CALL create_fem_parttwo(mo_fem)

    end subroutine create_parttwo
    !---------------------------------------------------------------
    subroutine create_partthree()
        implicit none
        ! LOCAL

        CALL create_fem_partthree(mo_fem)

    end subroutine create_partthree
    !---------------------------------------------------------------
    subroutine create_partfour()
        implicit none
        ! LOCAL

        CALL Assembly_initialize(mo_ass,mo_fem)

    end subroutine create_partfour
    !---------------------------------------------------------------
    subroutine pyfem_free()
        implicit none

        CALL printlog("pyfem_free : Start", ai_dtllevel = 1)

        CALL Assembly_finalize(mo_ass,mo_fem)
        CALL free_fem(mo_fem)
        CALL closetracelog()

        CALL printlog("pyfem_free : End", ai_dtllevel = 1)

    end subroutine pyfem_free
    !---------------------------------------------------------------
    subroutine pyfem_print_grid(ai_grids_id, ai_id)
        implicit none
        INTEGER :: ai_grids_id
        INTEGER :: ai_id

        CALL print_grid(mo_fem % opo_grids (ai_grids_id) % opo_grid (ai_id))

    end subroutine pyfem_print_grid
    !----------------------------------------------------------------------------------------------
    subroutine pyfem_get_weights_grid(ai_grids_id, ai_id, apr_w, ai_nel, ai_size)
        implicit none
        INTEGER :: ai_grids_id
        INTEGER :: ai_id
        INTEGER :: ai_nel
        integer, intent ( in) :: ai_size
        real(8), dimension (ai_nel,ai_size), intent ( out) :: apr_w
        ! LOCAL

        CALL get_weights(mo_fem % opo_grids (ai_grids_id) % opo_grid (ai_id), apr_w)

    end subroutine pyfem_get_weights_grid
    !---------------------------------------------------------------
    subroutine field_setval_coefs (ai_id, ar_val)
        implicit none
        INTEGER :: ai_id
        REAL(8) :: ar_val

        CALL field_setval_coefs_fem(mo_fem, ai_id, ar_val)

    end subroutine field_setval_coefs
    !----------------------------------------------------------------------------------------------
    subroutine field_setarray_coefs(ai_id, apr_val, ai_size)
        implicit none
        integer, intent(in) :: ai_id
        integer, intent (in) :: ai_size
        real(8), dimension (ai_size), intent ( in) :: apr_val

        CALL field_setarray_coefs_fem(mo_fem, ai_id, apr_val, ai_size)
        
    end subroutine field_setarray_coefs
    !---------------------------------------------------------------
    subroutine field_setval_values (ai_id, ar_val)
        implicit none
        INTEGER :: ai_id
        REAL(8) :: ar_val

        CALL field_setval_values_fem(mo_fem, ai_id, ar_val)

    end subroutine field_setval_values
    !----------------------------------------------------------------------------------------------
    subroutine field_setarray_values(ai_id, apr_val, ai_ndof, ai_maxnderiv, ai_maxnel, ai_maxnpts)
        implicit none
        integer, intent(in) :: ai_id
        integer, intent (in) :: ai_ndof
        integer, intent (in) :: ai_maxnderiv
        integer, intent (in) :: ai_maxnel
        integer, intent (in) :: ai_maxnpts
        real(8), dimension (ai_ndof, ai_maxnderiv, ai_maxnel, ai_maxnpts), intent ( in) :: apr_val

        CALL field_setarray_values_fem(mo_fem, ai_id, apr_val, ai_ndof, ai_maxnderiv, ai_maxnel, ai_maxnpts)

    end subroutine field_setarray_values
    !----------------------------------------------------------------------------------------------
    subroutine getfield(ai_field, apr_result, ai_size)
        implicit none
        integer, intent(in) :: ai_field
        integer, intent ( in) :: ai_size
        real(8), dimension ( ai_size), intent ( out) :: apr_result

        CALL getfield_fem(mo_fem,ai_field, apr_result, ai_size)

    end subroutine getfield
    !----------------------------------------------------------------------------------------------
    subroutine field_addition2(ai_i, ai_j, apr_result, ai_size)
        implicit none
        integer, intent(in) :: ai_i
        integer, intent(in) :: ai_j
        integer, intent ( in) :: ai_size
        real(8), dimension ( ai_size), intent ( out) :: apr_result

        CALL addop_field2(mo_fem,ai_i, ai_j, apr_result, ai_size)

    end subroutine field_addition2
    !----------------------------------------------------------------------------------------------
    subroutine field_multiplication2(ai_i, ai_j, apr_result, ai_size)
        implicit none
        integer, intent(in) :: ai_i
        integer, intent(in) :: ai_j
        integer, intent ( in) :: ai_size
        real(8), dimension ( ai_size), intent ( out) :: apr_result

        CALL multop_field2(mo_fem,ai_i, ai_j, apr_result, ai_size)

    end subroutine field_multiplication2
    !----------------------------------------------------------------------------------------------
    subroutine field_scaladdition2(ai_i, ar_val, apr_result, ai_size)
        implicit none
        integer, intent(in) :: ai_i
        real(8), intent(in) :: ar_val
        integer, intent ( in) :: ai_size
        real(8), dimension ( ai_size), intent ( out) :: apr_result

        CALL addscalop_field2(mo_fem,ai_i, ar_val, apr_result, ai_size)

    end subroutine field_scaladdition2
    !----------------------------------------------------------------------------------------------
    subroutine field_scalmultiplication2(ai_i, ar_val, apr_result, ai_size)
        implicit none
        integer, intent(in) :: ai_i
        real(8), intent(in) :: ar_val
        integer, intent ( in) :: ai_size
        real(8), dimension ( ai_size), intent ( out) :: apr_result

        CALL multscalop_field2(mo_fem,ai_i, ar_val, apr_result, ai_size)

    end subroutine field_scalmultiplication2
    !---------------------------------------------------------------
    subroutine field_addition (ai_i, ai_j, ai_k)
        implicit none
        INTEGER :: ai_i
        INTEGER :: ai_j
        INTEGER :: ai_k

        CALL addop_field (mo_fem, ai_i, ai_j, ai_k)

    end subroutine field_addition
    !---------------------------------------------------------------
    subroutine field_mult_field (ai_i, ai_j, ai_k)
        implicit none
        INTEGER :: ai_i
        INTEGER :: ai_j
        INTEGER :: ai_k

        CALL multop_field (mo_fem, ai_i, ai_j, ai_k)

    end subroutine field_mult_field
    !---------------------------------------------------------------
    subroutine field_mult_array (ai_field_id, apr_val_in, ai_size_in)
        implicit none
        INTEGER :: ai_field_id
        INTEGER :: ai_size_in
        real(8), dimension (ai_size_in), intent (in) :: apr_val_in

        CALL multop_field_array  (mo_fem, ai_field_id, apr_val_in, ai_size_in)

    end subroutine field_mult_array
    !---------------------------------------------------------------
    subroutine field_add_scal (ai_i, ar_val, ai_k)
        implicit none
        INTEGER :: ai_i
        REAL(8) :: ar_val
        INTEGER :: ai_k

        CALL addop_field_scal (mo_fem, ai_i, ar_val, ai_k)

    end subroutine field_add_scal
    !---------------------------------------------------------------
    subroutine field_mult_scal (ai_i, ar_val, ai_k)
        implicit none
        INTEGER :: ai_i
        REAL(8) :: ar_val
        INTEGER :: ai_k

        CALL multop_field_scal (mo_fem, ai_i, ar_val, ai_k)

    end subroutine field_mult_scal
    !---------------------------------------------------------------
    subroutine matrix_mult_field (ai_i, ai_j, ai_k)
        implicit none
        INTEGER :: ai_i
        INTEGER :: ai_j
        INTEGER :: ai_k

        CALL multop_matrix_field  (mo_fem, ai_i, ai_j, ai_k)

    end subroutine matrix_mult_field
    !---------------------------------------------------------------
    subroutine matrix_mult_array (ai_matrix_id, ai_field_id, apr_val, ai_size)
        implicit none
        INTEGER :: ai_matrix_id
        INTEGER :: ai_field_id
        INTEGER :: ai_size
        real(8), dimension (ai_size), intent (out) :: apr_val

        CALL multop_matrix_array  (mo_fem, ai_matrix_id, ai_field_id, apr_val, ai_size)

    end subroutine matrix_mult_array
    !---------------------------------------------------------------
    subroutine matrix_mult_array2 (ai_matrix_id, apr_val_in, ai_size_in, apr_val_out, ai_size_out)
        implicit none
        INTEGER :: ai_matrix_id
        INTEGER :: ai_size_in
        real(8), dimension (ai_size_in), intent (in) :: apr_val_in
        INTEGER :: ai_size_out
        real(8), dimension (ai_size_out), intent (out) :: apr_val_out

        CALL multop_matrix_array2  (mo_fem, ai_matrix_id, apr_val_in, ai_size_in, apr_val_out, ai_size_out)

    end subroutine matrix_mult_array2
    !---------------------------------------------------------------
    subroutine matrix_add_matrix (ai_i, ai_j, ai_k)
        implicit none
        INTEGER, INTENT(IN) :: ai_i
        INTEGER, INTENT(IN) :: ai_j
        INTEGER, INTENT(IN) :: ai_k

        CALL addop_matrix_matrix  (mo_fem, ai_i, ai_j, ai_k)

    end subroutine matrix_add_matrix
    !---------------------------------------------------------------
    subroutine matrix_mult_matrix (ai_i, ai_j, ai_k)
        implicit none
        INTEGER, INTENT(IN) :: ai_i
        INTEGER, INTENT(IN) :: ai_j
        INTEGER, INTENT(IN) :: ai_k

        CALL multop_matrix_matrix  (mo_fem, ai_i, ai_j, ai_k)

    end subroutine matrix_mult_matrix
    !---------------------------------------------------------------
    subroutine matrix_add_scal (ai_i, ar_val, ai_k)
        implicit none
        INTEGER :: ai_i
        REAL(8) :: ar_val
        INTEGER :: ai_k

        CALL addop_matrix_scal  (mo_fem, ai_i, ar_val, ai_k)

    end subroutine matrix_add_scal
    !---------------------------------------------------------------
    subroutine matrix_mult_scal (ai_i, ar_val, ai_k)
        implicit none
        INTEGER :: ai_i
        REAL(8) :: ar_val
        INTEGER :: ai_k

        CALL multop_matrix_scal  (mo_fem, ai_i, ar_val, ai_k)

    end subroutine matrix_mult_scal
    !----------------------------------------------------------------------------------------------
    subroutine solver_setvaluesrhs(ai_matrix_id, apr_val, ai_size)
        implicit none
        integer, intent(in) :: ai_matrix_id
        integer, intent (in) :: ai_size
        real(8), dimension (ai_size), intent ( in) :: apr_val

        CALL solver_setvaluesrhs_fem(mo_fem, ai_matrix_id, apr_val, ai_size)

    end subroutine solver_setvaluesrhs
    !----------------------------------------------------------------------------------------------
    subroutine solver_setfieldrhs(ai_matrix_id, ai_field_id)
        implicit none
        integer, intent(in) :: ai_matrix_id
        integer, intent(in) :: ai_field_id

        CALL solver_setfieldrhs_fem(mo_fem, ai_matrix_id, ai_field_id)

    end subroutine solver_setfieldrhs
    !----------------------------------------------------------------------------------------------
    subroutine solver_getsolution(ai_matrix_id, apr_val, ai_size)
        implicit none
        integer, intent(in) :: ai_matrix_id
        integer, intent ( in) :: ai_size
        real(8), dimension ( ai_size), intent ( out) :: apr_val

        CALL solver_getsolution_fem(mo_fem,ai_matrix_id, apr_val, ai_size)

    end subroutine solver_getsolution
    !----------------------------------------------------------------------------------------------
    subroutine solver_getfieldsolution(ai_matrix_id, ai_field_id)
        implicit none
        integer, intent(in) :: ai_matrix_id
        integer, intent(in) :: ai_field_id

        CALL solver_getfieldsolution_fem(mo_fem, ai_matrix_id, ai_field_id)

    end subroutine solver_getfieldsolution
        !----------------------------------------------------------------------------------------------
    subroutine solver_initialize(ai_solver_id)
        implicit none
        integer, intent(in) :: ai_solver_id

        CALL solver_initialize_fem(mo_fem, ai_solver_id)

    end subroutine solver_initialize    
        !----------------------------------------------------------------------------------------------
    subroutine solver_run(ai_solver_id)
        implicit none
        integer, intent(in) :: ai_solver_id

        CALL solver_run_fem(mo_fem, ai_solver_id)

    end subroutine solver_run
        !----------------------------------------------------------------------------------------------
    subroutine solver_finalize(ai_solver_id)
        implicit none
        integer, intent(in) :: ai_solver_id

        CALL solver_finalize_fem(mo_fem, ai_solver_id)

    end subroutine solver_finalize     
    !---------------------------------------------------------------
    subroutine set_patch_space(ai_space, ai_id, ai_rational, api_N, api_P, apr_u, apr_P, apr_W, ai_dim, ai_Rd, ai_maxnk, ai_npoint)
        !> this routine creates and allocates memory for the object geometry
        implicit none
        integer, intent(in) :: ai_space
        integer, intent(in) :: ai_id
        integer, intent(in) :: ai_rational
        integer :: ai_dim
        integer :: ai_Rd
        integer :: ai_maxnk
        integer :: ai_npoint
        integer, dimension(ai_dim) :: api_N
        integer, dimension(ai_dim) :: api_P
        real(8), dimension(ai_maxnk, ai_dim) :: apr_u
        real(8), dimension(ai_npoint, ai_Rd) :: apr_P
        real(8), dimension(ai_npoint) :: apr_W
        ! LOCAL
        type(GEOMETRY) :: lo_geo ! the current geometry

        CALL printlog("set_patch_space : Start", ai_dtllevel = 1)

        ! TRANSPOSE a enlever, apres le changement des acces tableaux dans geometry pour U et P
        CALL create_geometry(lo_geo, ai_dim, ai_Rd, ai_rational, api_N, api_P, TRANSPOSE(apr_u), TRANSPOSE(apr_P), apr_W)

        ! id =0,npatchs-1, but in geometries it must be =1,npatchs
        CALL add_geo(mo_fem % opo_spaces (ai_space) % oo_mapping, lo_geo, ai_id+1)

        CALL free_geometry(lo_geo)

        CALL printlog("set_patch_space : End", ai_dtllevel = 1)

    end subroutine set_patch_space
    !---------------------------------------------------------------
    subroutine set_patch_mapping(ai_mapping, ai_id, ai_rational, api_N, api_P &
    , apr_u, apr_P, apr_W, ai_dim, ai_Rd, ai_maxnk, ai_npoint)
        !> this routine creates and allocates memory for the object geometry
        implicit none
        integer, intent(in) :: ai_mapping
        integer, intent(in) :: ai_id
        integer, intent(in) :: ai_rational
        integer :: ai_dim
        integer :: ai_Rd
        integer :: ai_maxnk
        integer :: ai_npoint
        integer, dimension(ai_dim) :: api_N
        integer, dimension(ai_dim) :: api_P
        real(8), dimension(ai_maxnk, ai_dim) :: apr_u
        real(8), dimension(ai_npoint, ai_Rd) :: apr_P
        real(8), dimension(ai_npoint) :: apr_W
        ! LOCAL
        type(GEOMETRY) :: lo_geo ! the current geometry

        CALL printlog("set_patch_mapping : Start", ai_dtllevel = 1)

        ! TRANSPOSE a enlever, apres le changement des acces tableaux dans geometry pour U et P
        CALL create_geometry(lo_geo, ai_dim, ai_Rd, ai_rational, api_N, api_P, TRANSPOSE(apr_u), TRANSPOSE(apr_P), apr_W)
        
        CALL add_geo(mo_fem % opo_mappings (ai_mapping) % oo_mapping, lo_geo, ai_id+1)
        
        CALL free_geometry(lo_geo)

        CALL printlog("set_patch_mapping : End", ai_dtllevel = 1)

    end subroutine set_patch_mapping
    !---------------------------------------------------------------
    subroutine update_points_patch_space(ai_space, ai_id, apr_P, apr_W, ai_Rd, ai_npoint)
        !> this routine creates and allocates memory for the object geometry
        implicit none
        integer, intent(in) :: ai_space
        integer, intent(in) :: ai_id
        integer :: ai_Rd
        integer :: ai_npoint
        real(8), dimension(ai_npoint, ai_Rd) :: apr_P
        real(8), dimension(ai_npoint) :: apr_W
        ! LOCAL

        CALL printlog("update_points_patch_space: Start", ai_dtllevel = 1)

        CALL update_points_geometry(&
        mo_fem % opo_spaces (ai_space) % oo_mapping % opo_geo (ai_id+1)&
        , TRANSPOSE(apr_P), apr_W)

        CALL printlog("update_points_patch_space: End", ai_dtllevel = 1)

    end subroutine update_points_patch_space    
    !---------------------------------------------------------------
    subroutine update_points_patch_mapping(ai_mapping, ai_id, apr_P, apr_W, ai_Rd, ai_npoint)
        !> this routine creates and allocates memory for the object geometry
        implicit none
        integer, intent(in) :: ai_mapping
        integer, intent(in) :: ai_id
        integer :: ai_Rd
        integer :: ai_npoint
        real(8), dimension(ai_npoint, ai_Rd) :: apr_P
        real(8), dimension(ai_npoint) :: apr_W
        ! LOCAL

        CALL printlog("update_points_patch_mapping: Start", ai_dtllevel = 1)

        CALL update_points_geometry(&
        mo_fem % opo_mappings (ai_mapping) % oo_mapping % opo_geo (ai_id+1)&
        , TRANSPOSE(apr_P), apr_W)        

        CALL printlog("supdate_points_patch_mapping: End", ai_dtllevel = 1)

    end subroutine update_points_patch_mapping
    !---------------------------------------------------------------
    subroutine pyfem_create_from_tensor_product_1d ( ai_grids_id, ai_id         &
                                                    , apr_localx, apr_localwx   &
                                                    , ai_nx, ai_maxnptsp1x      )
        implicit none
        integer :: ai_grids_id
        integer :: ai_id
        integer :: ai_nx
        integer :: ai_maxnptsp1x
        real(8), dimension(ai_nx,0:ai_maxnptsp1x-1) :: apr_localx
        real(8), dimension(ai_nx,ai_maxnptsp1x-1) :: apr_localwx

        CALL printlog("pyfem_create_from_tensor_product_1d : Start", ai_dtllevel = 1)

        CALL create_from_tensor_product_1D(mo_fem % opo_grids (ai_grids_id) % opo_grid(ai_id)   &
        ,apr_localx, apr_localwx)

        CALL printlog("pyfem_create_from_tensor_product_1d : End", ai_dtllevel = 1)

    end subroutine pyfem_create_from_tensor_product_1d
    !---------------------------------------------------------------
    subroutine pyfem_create_from_tensor_product_2d ( ai_grids_id, ai_id       &
                                                    , apr_localx, apr_localwx &
                                                    , apr_localy, apr_localwy &
                                                    , ai_nx, ai_maxnptsp1x    &
                                                    , ai_ny, ai_maxnptsp1y    )
        implicit none
        integer :: ai_grids_id
        integer :: ai_id
        integer :: ai_nx
        integer :: ai_maxnptsp1x
        integer :: ai_ny
        integer :: ai_maxnptsp1y
        real(8), dimension(ai_nx,0:ai_maxnptsp1x-1) :: apr_localx
        real(8), dimension(ai_nx,ai_maxnptsp1x-1) :: apr_localwx
        real(8), dimension(ai_ny,0:ai_maxnptsp1y-1) :: apr_localy
        real(8), dimension(ai_ny,ai_maxnptsp1y-1) :: apr_localwy

        CALL printlog("pyfem_create_from_tensor_product_2d : Start", ai_dtllevel = 1)

        CALL create_from_tensor_product_2D(mo_fem % opo_grids (ai_grids_id) % opo_grid(ai_id)   &
        , apr_localx, apr_localy, apr_localwx, apr_localwy)

        CALL printlog("pyfem_create_from_tensor_product_2d : End", ai_dtllevel = 1)

    end subroutine pyfem_create_from_tensor_product_2d
    !---------------------------------------------------------------
    subroutine pyfem_create_from_tensor_product_3d ( ai_grids_id, ai_id       &
                                                    , apr_localx, apr_localwx &
                                                    , apr_localy, apr_localwy &
                                                    , apr_localz, apr_localwz &
                                                    , ai_nx, ai_maxnptsp1x    &
                                                    , ai_ny, ai_maxnptsp1y    &
                                                    , ai_nz, ai_maxnptsp1z    )
        implicit none
        integer :: ai_grids_id
        integer :: ai_id
        integer :: ai_nx
        integer :: ai_maxnptsp1x
        integer :: ai_ny
        integer :: ai_maxnptsp1y
        integer :: ai_nz
        integer :: ai_maxnptsp1z
        real(8), dimension(ai_nx,0:ai_maxnptsp1x-1) :: apr_localx
        real(8), dimension(ai_nx,ai_maxnptsp1x-1) :: apr_localwx
        real(8), dimension(ai_ny,0:ai_maxnptsp1y-1) :: apr_localy
        real(8), dimension(ai_ny,ai_maxnptsp1y-1) :: apr_localwy
        real(8), dimension(ai_ny,0:ai_maxnptsp1y-1) :: apr_localz
        real(8), dimension(ai_nz,ai_maxnptsp1z-1) :: apr_localwz

        CALL printlog("pyfem_create_from_tensor_product_3d : Start", ai_dtllevel = 1)

        CALL create_from_tensor_product_3D(mo_fem % opo_grids (ai_grids_id) % opo_grid(ai_id)   &
        , apr_localx, apr_localy, apr_localz  &
        , apr_localwx, apr_localwy, apr_localwz)

        CALL printlog("pyfem_create_from_tensor_product_3d : End", ai_dtllevel = 1)

    end subroutine pyfem_create_from_tensor_product_3d
    !---------------------------------------------------------------
    subroutine pyfem_set_unidirection_tensor_2d ( ai_grids_id, ai_id, ai_direction, ai_curelt  &
        , ar_fixpt, ar_fixw &
        , apr_dirpts, apr_dirw  &
        , ai_newelt    &
        , ai_n, ai_maxnptsp1)
        implicit none
        integer :: ai_grids_id
        integer :: ai_id
        integer :: ai_direction
        !> ai_elt is the id of the last created element
        integer :: ai_curelt
        real(8) :: ar_fixpt
        real(8) :: ar_fixw
        integer :: ai_n
        integer :: ai_maxnptsp1
        real(8), dimension(ai_n,0:ai_maxnptsp1-1) :: apr_dirpts
        real(8), dimension(ai_n,ai_maxnptsp1-1) :: apr_dirw
        !> ai_elt is the id of the last created element
        integer, intent(out) :: ai_newelt
        
        CALL printlog("pyfem_set_unidirection_tensor_2d : Start", ai_dtllevel = 1)

        CALL set_unidirection_tensor_2d(mo_fem % opo_grids (ai_grids_id) % opo_grid(ai_id)   &
        , ai_direction  &
        , apr_dirpts, apr_dirw  &
        , ar_fixpt, ar_fixw, ai_curelt)

        ai_newelt = ai_curelt
        
        CALL printlog("pyfem_set_unidirection_tensor_2d : End", ai_dtllevel = 1)

    end subroutine pyfem_set_unidirection_tensor_2d
    !---------------------------------------------------------------
    subroutine pyfem_set_tensor_basis ( ai_grids_id, ai_id, ai_d, apr_dbasisatx, ai_nderiv, ai_p, ai_ni, ai_dirnpts &
        , ai_maxnderiv, ai_maxp, ai_maxni, ai_maxdirnpts )
        implicit none
        integer :: ai_grids_id
        integer :: ai_id
        INTEGER :: ai_d
        INTEGER :: ai_nderiv
        INTEGER :: ai_p
        INTEGER :: ai_ni
        INTEGER :: ai_dirnpts
        REAL(8), DIMENSION(0:ai_nderiv, 0:ai_p, 1:ai_dirnpts, 1:ai_ni) :: apr_dbasisatx
        INTEGER :: ai_maxnderiv
        INTEGER :: ai_maxp
        INTEGER :: ai_maxni
        INTEGER :: ai_maxdirnpts
        ! LOCAL

        CALL printlog("pyfem_set_tensor_basis : Start", ai_dtllevel = 1)

        CALL set_tensor_basis ( mo_fem % opo_grids (ai_grids_id) % opo_grid(ai_id) &
                                    , ai_d, apr_dbasisatx, ai_nderiv, ai_p, ai_ni, ai_dirnpts &
                                    , ai_maxnderiv, ai_maxp, ai_maxni, ai_maxdirnpts )
        
        CALL printlog("pyfem_set_tensor_basis : End", ai_dtllevel = 1)

    end subroutine pyfem_set_tensor_basis
    !---------------------------------------------------------------
    subroutine pyfem_create_connectivity(ai_id, ai_npatch, ai_maxnen, ai_nelt)
        implicit none
        integer :: ai_id
        integer :: ai_npatch
        integer :: ai_maxnen
        integer :: ai_nelt
        ! LOCAL

        CALL printlog("pyfem_create_connectivity : Start", ai_dtllevel = 1)

        CALL create_connectivity(mo_fem % opo_spaces(ai_id) % oo_con, ai_npatch, ai_maxnen, ai_nelt)

        CALL printlog("pyfem_create_connectivity : End", ai_dtllevel = 1)
    end subroutine pyfem_create_connectivity    
    !---------------------------------------------------------------
    subroutine pyfem_set_connectivity_id(ai_id, api_ID, ai_N)
        implicit none
        integer :: ai_id
        integer :: ai_N
        integer, dimension(ai_N) :: api_ID
        ! LOCAL

        CALL printlog("pyfem_set_connectivity_id : Start", ai_dtllevel = 1)

        CALL set_connectivity_id(mo_fem % opo_spaces(ai_id) % oo_con, api_ID, ai_N  )

        CALL printlog("pyfem_set_connectivity_id : End", ai_dtllevel = 1)

    end subroutine pyfem_set_connectivity_id    
    !---------------------------------------------------------------
    subroutine pyfem_set_connectivity_real_elts(ai_id, api_real_elts, ai_npatch, ai_nelt)
        implicit none
        integer :: ai_id
        integer :: ai_npatch
        integer :: ai_nelt
        integer, dimension(ai_npatch, ai_nelt) :: api_real_elts
        ! LOCAL

        CALL printlog("pyfem_set_connectivity_real_elts : Start", ai_dtllevel = 1)

        CALL set_real_elts(mo_fem % opo_spaces(ai_id) % oo_con, api_real_elts, ai_npatch, ai_nelt)

        CALL printlog("pyfem_set_connectivity_real_elts : End", ai_dtllevel = 1)

    end subroutine pyfem_set_connectivity_real_elts
    !---------------------------------------------------------------
    subroutine pyfem_set_connectivity_lm(ai_id, api_LM, ai_npatch, ai_maxnen, ai_nrealelt)
        implicit none
        integer :: ai_id
        integer :: ai_npatch
        integer :: ai_maxnen
        integer :: ai_nrealelt
        integer, dimension(ai_npatch, ai_maxnen, ai_nrealelt) :: api_LM
        ! LOCAL

        CALL printlog("pyfem_set_connectivity_lm : Start", ai_dtllevel = 1)

        CALL set_connectivity_lm(mo_fem % opo_spaces(ai_id) % oo_con    &
        , api_LM, ai_npatch, ai_maxnen, ai_nrealelt )

        CALL printlog("pyfem_set_connectivity_lm : End", ai_dtllevel = 1)

    end subroutine pyfem_set_connectivity_lm    
    !---------------------------------------------------------------
    subroutine pyfem_set_connectivity_ien(ai_id, api_IEN, ai_npatch, ai_maxnen, ai_nrealelt)
        implicit none
        integer :: ai_id
        integer :: ai_npatch
        integer :: ai_maxnen
        integer :: ai_nrealelt
        integer, dimension(ai_npatch, ai_maxnen, ai_nrealelt) :: api_IEN
        ! LOCAL

        CALL printlog("pyfem_set_connectivity_ien : Start", ai_dtllevel = 1)

        CALL set_connectivity_ien(mo_fem % opo_spaces(ai_id) % oo_con    &
        , api_IEN , ai_npatch, ai_maxnen, ai_nrealelt)

        CALL printlog("pyfem_set_connectivity_ien : End", ai_dtllevel = 1)

    end subroutine pyfem_set_connectivity_ien    
    !---------------------------------------------------------------
    subroutine pyfem_set_connectivity_loc(ai_id, ai_patch, api_id_loc, ai_N)
        implicit none
        integer :: ai_id
        integer :: ai_patch
        integer :: ai_N
        integer, dimension(ai_n) :: api_id_loc
        ! LOCAL

        CALL printlog("pyfem_set_connectivity_loc: Start", ai_dtllevel = 1)

        CALL create_ID_loc(mo_fem % opo_spaces(ai_id) % oo_con, ai_patch, ai_n, api_id_loc)

        CALL printlog("pyfem_set_connectivity_loc: End", ai_dtllevel = 1)

    end subroutine pyfem_set_connectivity_loc    
    !---------------------------------------------------------------
    subroutine pyfem_assembly(ai_grids_id)
        implicit none
        INTEGER, INTENT(IN) :: ai_grids_id

        CALL printlog("pyfem_assembly : Start", ai_dtllevel = 1)

        CALL Assembly_Matrix(mo_ass,mo_fem, ai_grids_id)

        CALL printlog("pyfem_assembly : End", ai_dtllevel = 1)

    end subroutine pyfem_assembly
!----------------------------------------------------------------------------------------------
    subroutine set_metric_points ( ai_ref, ai_patch, apr_F, apr_DF, ai_nel, ai_npts, ai_Rd, ai_Rd2 )
        implicit none
        integer, intent(in)  :: ai_ref  !THE METRIC ID
        integer, intent(in)  :: ai_patch  !THE ai_patch ID
        integer, intent(in)  :: ai_nel   !NUMBER OF ELEMENTS
        integer, intent(in)  :: ai_npts   !NUMBER OF POINTS
        integer, intent(in)  :: ai_Rd   !THE DOMAIN DIMENSION
        integer, intent(in)  :: ai_Rd2  !GENERALLY THIS ai_Rd**2
        real*8 , dimension(ai_nel, ai_npts,ai_Rd), intent(in)  :: apr_F
        real*8 , dimension(ai_nel, ai_npts,ai_Rd2), intent(in)  :: apr_DF
        ! LOCAL VARIABLES

        CALL printlog("set_metric_points : Begin", ai_dtllevel = 1)

        CALL set_metric_points_ass ( mo_ass, mo_fem, ai_ref, ai_patch, apr_F, apr_DF, ai_nel, ai_npts, ai_Rd, ai_Rd2 )

        CALL printlog("set_metric_points : End", ai_dtllevel = 1)

    end subroutine set_metric_points
!----------------------------------------------------------------------------------------------
    subroutine set_metric_points_advanced ( ai_ref, ai_patch, apr_points, ai_nel, ai_npts, ai_Rd, ai_nderivp1 )
        implicit none
        integer, intent(in)  :: ai_ref  !THE METRIC ID
        integer, intent(in)  :: ai_patch  !THE ai_patch ID
        integer, intent(in)  :: ai_nel   !NUMBER OF ELEMENTS
        integer, intent(in)  :: ai_npts   !NUMBER OF POINTS
        integer, intent(in)  :: ai_Rd   !THE DOMAIN DIMENSION
        integer, intent(in)  :: ai_nderivp1  !GENERALLY THIS ai_Rd
        real*8 , dimension(ai_nel, ai_npts,ai_Rd, 0:ai_nderivp1-1), intent(in)  :: apr_points
        ! LOCAL VARIABLES

        CALL printlog("set_metric_points_advanced : Begin", ai_dtllevel = 1)

        CALL set_metric_points_advanced_ass ( mo_ass, mo_FEM, ai_ref, ai_patch, apr_points, ai_nel, ai_npts, ai_Rd, ai_nderivp1 )

        CALL printlog("set_metric_points_advanced : End", ai_dtllevel = 1)

    end subroutine set_metric_points_advanced
!----------------------------------------------------------------------------------------------
    subroutine get_space_points_advanced ( ai_space, ai_patch, apr_points, ai_nelts, ai_npts, ai_dim, ai_nderiv )
        implicit none
        integer, intent(in)  :: ai_space
        integer, intent(in)  :: ai_patch
        integer, intent(in)  :: ai_nelts
        integer, intent(in)  :: ai_npts
        integer, intent(in)  :: ai_dim
        integer, intent(in)  :: ai_nderiv
        real*8 , dimension(ai_nelts,ai_npts,ai_dim,0:ai_nderiv), intent(out)  :: apr_points
        ! LOCAL VARIABLES

        CALL printlog("get_space_points_advanced : Begin", ai_dtllevel = 1)

        CALL assembly_space_points_advanced ( mo_ass, mo_fem, ai_space, ai_patch, ai_nderiv, apr_points )
        
        CALL printlog("get_space_points_advanced : End", ai_dtllevel = 1)

    end subroutine get_space_points_advanced
!----------------------------------------------------------------------------------------------
    subroutine get_space_parametricpoints ( ai_space, ai_patch, apr_points, ai_nelts, ai_npts, ai_dim )
        implicit none
        integer, intent(in)  :: ai_space
        integer, intent(in)  :: ai_patch
        integer, intent(in)  :: ai_nelts
        integer, intent(in)  :: ai_npts
        integer, intent(in)  :: ai_dim
        real*8 , dimension(ai_nelts,ai_npts,ai_dim), intent(out)  :: apr_points
        ! LOCAL VARIABLES

        CALL printlog("get_space_parametricpoints : Begin", ai_dtllevel = 1)

        CALL get_space_parametric_points ( mo_ass, mo_fem, ai_space, ai_patch, apr_points, ai_nelts, ai_npts, ai_dim )

        CALL printlog("get_space_parametricpoints : End", ai_dtllevel = 1)

    end subroutine get_space_parametricpoints
!----------------------------------------------------------------------------------------------
    subroutine set_field_on_grids ( ai_field, ai_patch, apr_values, ai_nparam, ai_nel, ai_npts )
        implicit none
        integer, intent(in)  :: ai_field
        integer, intent(in)  :: ai_patch
        integer, intent(in)  :: ai_nparam
        integer, intent(in)  :: ai_nel
        integer, intent(in)  :: ai_npts
        real*8 , dimension(ai_nparam, ai_nel, ai_npts), intent(in)  :: apr_values
        ! LOCAL VARIABLES
        integer :: li_space
        integer :: li_grids
        integer :: li_locid

        CALL printlog("set_field_on_grids : Start", ai_dtllevel = 1)

!#ifdef _DEBUG
        call concatmsg("current field is = ", ai_dtllevel = mi_dtllevel_base + 1)
        call concatmsg(ai_field, ai_dtllevel = mi_dtllevel_base + 1)
        call printmsg(ai_dtllevel = mi_dtllevel_base + 1)
!#endif

        li_space = mo_fem % opi_InfoField ( ai_field, INFOFIELD_SPACE )
        li_grids = mo_fem % opi_InfoSpace ( li_space, INFOSPACE_GRIDS )
        li_locid = mo_fem % opi_InfoField ( ai_field, INFOFIELD_LOCID )

        CALL set_values_field( mo_fem % opo_grids(li_grids) % opo_grid ( ai_patch ), li_locid, apr_values)

        CALL printlog("set_field_on_grids : End", ai_dtllevel = 1)

    end subroutine set_field_on_grids
!----------------------------------------------------------------------------------------------
    subroutine set_operator_on_grids ( ai_operator, ai_patch, apr_values, ai_nparams, ai_nel, ai_npts)
        implicit none
        integer, intent(in)  :: ai_operator
        integer, intent(in)  :: ai_patch
        integer, intent(in)  :: ai_nparams
        integer, intent(in)  :: ai_nel
        integer, intent(in)  :: ai_npts
        real*8 , dimension(ai_nparams, ai_nel, ai_npts) :: apr_values
        ! LOCAL VARIABLES
        integer :: li_space
        integer :: li_grids
        integer :: li_locid

        CALL printlog("set_operator_on_grids : Start", ai_dtllevel = 1)
        
        li_space = mo_fem % opi_InfoOperator ( ai_operator, INFOOPERATOR_SPACE_1 )
        li_grids = mo_fem % opi_InfoSpace ( li_space, INFOSPACE_GRIDS )
        li_locid = mo_fem % opi_InfoOperator ( ai_operator, INFOOPERATOR_LOCID )

        CALL set_values_matrix( mo_fem % opo_grids(li_grids) % opo_grid ( ai_patch ), li_locid, apr_values)

        CALL printlog("set_operator_on_grids : End", ai_dtllevel = 1)

    end subroutine set_operator_on_grids
!----------------------------------------------------------------------------------------------
    subroutine set_norm_on_grids ( ai_norm, ai_patch, apr_values, ai_nparams, ai_nel, ai_npts)
        implicit none
        integer, intent(in)  :: ai_norm
        integer, intent(in)  :: ai_patch
        integer, intent(in)  :: ai_nparams
        integer, intent(in)  :: ai_nel
        integer, intent(in)  :: ai_npts
        real*8 , dimension(ai_nparams, ai_nel, ai_npts) :: apr_values
        ! LOCAL VARIABLES
        integer :: li_field
        integer :: li_space
        integer :: li_grids
        integer :: li_locid

        CALL printlog("set_norm_on_grids : Start", ai_dtllevel = 1)

        li_field = mo_fem % opi_InfoNorm ( ai_norm, INFONORM_FIELD )
        li_space = mo_fem % opi_InfoField ( li_field, INFOFIELD_SPACE )
        li_grids = mo_fem % opi_InfoSpace ( li_space, INFOSPACE_GRIDS )
        li_locid = mo_fem % opi_InfoNorm ( ai_norm, INFONORM_LOCID )

        CALL set_values_norm( mo_fem % opo_grids(li_grids) % opo_grid ( ai_patch ), li_locid, apr_values)

        CALL printlog("set_norm_on_grids : End", ai_dtllevel = 1)

    end subroutine set_norm_on_grids
!----------------------------------------------------------------------------------------------
    subroutine set_space_weights ( ai_space, ai_patch, ai_direction, apr_values, ai_dim)
        implicit none
        integer, intent(in)  :: ai_space
        integer, intent(in)  :: ai_patch
        integer, intent(in)  :: ai_direction
        integer, intent(in)  :: ai_dim
        real*8 , dimension(ai_dim), intent(in) :: apr_values
        ! LOCAL VARIABLES

        CALL printlog("set_space_weights : Start", ai_dtllevel = 1)

        CALL set_space_weights_fem ( mo_fem, ai_space, ai_patch, ai_direction, apr_values, ai_dim)

        CALL printlog("set_space_weights : End", ai_dtllevel = 1)

    end subroutine set_space_weights
!----------------------------------------------------------------------------------------------
    subroutine eval_field_on_grids ( ai_field, ai_patch, apr_fieldh, ai_ndof, ai_maxder, ai_nel, ai_maxnpts )
        implicit none
        integer, intent(in)  :: ai_field
        integer, intent(in)  :: ai_patch
        integer, intent(in)  :: ai_ndof
        integer, intent(in)  :: ai_maxder
        integer, intent(in)  :: ai_nel
        integer, intent(in)  :: ai_maxnpts
        real*8 , dimension(ai_ndof, ai_maxder, ai_nel, ai_maxnpts), intent(out)  :: apr_fieldh
        ! LOCAL VARIABLES

        CALL printlog("eval_field_on_grids : Start", ai_dtllevel = 1)

        CALL assembly_field_on_grids ( mo_ass, mo_fem, ai_field, ai_patch, apr_fieldh   &
        , ai_ndof, ai_maxder, ai_nel, ai_maxnpts )

        CALL printlog("eval_field_on_grids : End", ai_dtllevel = 1)

    end subroutine eval_field_on_grids
!----------------------------------------------------------------------------------------------
    subroutine pyfem_get_fieldh ( ai_field, ai_patch, apr_fieldh, ai_k, ai_ndof, ai_maxder, ai_nel, ai_maxnpts )
        !> ai_k must be equal to ai_ndof * ai_maxder
        implicit none
        integer, intent(in)  :: ai_field
        integer, intent(in)  :: ai_patch
        integer, intent(in)  :: ai_k
        integer, intent(in)  :: ai_ndof
        integer, intent(in)  :: ai_maxder
        integer, intent(in)  :: ai_nel
        integer, intent(in)  :: ai_maxnpts
        real*8 , dimension( ai_k, ai_nel, ai_maxnpts), intent(out)  :: apr_fieldh
        ! LOCAL VARIABLES

        CALL printlog("pyfem_get_fieldh : Start", ai_dtllevel = 1)

        CALL get_fieldh ( mo_fem, ai_field, ai_patch, apr_fieldh, ai_k, ai_ndof, ai_maxder, ai_nel, ai_maxnpts )

        CALL printlog("pyfem_get_fieldh : End", ai_dtllevel = 1)

    end subroutine pyfem_get_fieldh
    !----------------------------------------------------------------------------------------------
    subroutine getGlobalnorm(ai_norm, ar_result)
        implicit none
        integer, intent(in) :: ai_norm
        real(8), intent ( out) :: ar_result

        CALL getGlobalNorm_fem(mo_fem, ai_norm, ar_result)

    end subroutine getGlobalnorm
    !----------------------------------------------------------------------------------------------
    subroutine getPatchNorm(ai_norm, apr_result, ai_npatchs)
        implicit none
        integer, intent(in) :: ai_norm
        integer, intent ( in) :: ai_npatchs
        real(8), dimension ( ai_npatchs), intent ( out) :: apr_result

        CALL getPatchNorm_fem(mo_fem, ai_norm, apr_result, ai_npatchs)

    end subroutine getPatchNorm
    !----------------------------------------------------------------------------------------------
    subroutine getElementNorm(ai_norm, ai_patch, apr_result, ai_nel)
        implicit none
        integer, intent(in) :: ai_norm
        integer, intent(in) :: ai_patch
        integer, intent(in) :: ai_nel
        real(8), dimension (ai_nel), intent ( out) :: apr_result

        CALL getElementNorm_fem(mo_fem, ai_norm, ai_patch, apr_result, ai_nel)

    end subroutine getElementNorm
!----------------------------------------------------------------------------------------------
    subroutine eval_field_on_grids_2d ( ai_field, ai_patch, apr_fieldh, ai_ndof, ai_nel, ai_maxdirnpts )
        implicit none
        integer, intent(in)  :: ai_field
        integer, intent(in)  :: ai_patch
        integer, intent(in)  :: ai_ndof
        integer, intent(in)  :: ai_nel
        integer, intent(in)  :: ai_maxdirnpts
        real*8 , dimension(ai_ndof, ai_nel, ai_maxdirnpts, ai_maxdirnpts), intent(out)  :: apr_fieldh
        ! LOCAL VARIABLES

        CALL printlog("eval_field_on_grids_2d : Start", ai_dtllevel = 1)

        CALL assembly_field_on_grids_2d ( mo_ass, mo_fem, ai_field, ai_patch, apr_fieldh    &
        , ai_ndof, ai_nel, ai_maxdirnpts )

        CALL printlog("eval_field_on_grids_2d : End", ai_dtllevel = 1)

    end subroutine eval_field_on_grids_2d
!----------------------------------------------------------------------------------------------
    subroutine pyfem_reset_field(ai_field)
        implicit none
        integer, intent(in)  :: ai_field
        ! LOCAL VARIABLES

        CALL printlog("pyfem_reset_field: Start", ai_dtllevel = 1)

        CALL reset_field( mo_fem, ai_field )

        CALL printlog("pyfem_reset_field: End", ai_dtllevel = 1)

    end subroutine pyfem_reset_field
!----------------------------------------------------------------------------------------------
    subroutine pyfem_set_initialize_field( ai_field)
        implicit none
        integer, intent(in)  :: ai_field
        ! LOCAL VARIABLES

        CALL printlog("pyfem_set_initialize_field: Start", ai_dtllevel = 1)

        CALL set_initialize_field( mo_fem, ai_field )

        CALL printlog("pyfem_set_initialize_field: End", ai_dtllevel = 1)

    end subroutine pyfem_set_initialize_field
!----------------------------------------------------------------------------------------------
    subroutine pyfem_set_finalize_field( ai_field)
        implicit none
        integer, intent(in)  :: ai_field
        ! LOCAL VARIABLES

        CALL printlog("pyfem_set_finalize_field: Start", ai_dtllevel = 1)

        CALL set_finalize_field( mo_fem, ai_field )

        CALL printlog("pyfem_set_finalize_field: End", ai_dtllevel = 1)

    end subroutine pyfem_set_finalize_field
!----------------------------------------------------------------------------------------------
    subroutine field_to_matrix_1d ( ai_field, ai_patch, apr_matrix, ai_n_1 )
        implicit none
        integer, intent(in)  :: ai_field
        integer, intent(in)  :: ai_patch
        integer, intent(in)  :: ai_n_1
        real*8 , dimension(ai_n_1), intent(out)  :: apr_matrix
        ! LOCAL VARIABLES

        CALL printlog("field_to_matrix_1d : Start", ai_dtllevel = 1)

        CALL field_to_matrix_1D_fem ( mo_fem, ai_field, ai_patch, apr_matrix, ai_n_1 )

        CALL printlog("field_to_matrix_1d : End", ai_dtllevel = 1)

    end subroutine field_to_matrix_1d
!----------------------------------------------------------------------------------------------
    subroutine field_from_matrix_1d( ai_field, ai_patch, ar_alpha, apr_matrix, ai_n_1 )
        implicit none
        integer, intent(in)  :: ai_field
        integer, intent(in)  :: ai_patch
        REAL*8 , INTENT(IN)  :: ar_alpha
        integer, intent(in)  :: ai_n_1
        real*8 , dimension(ai_n_1), intent(in)  :: apr_matrix
        ! LOCAL VARIABLES

        CALL printlog("field_from_matrix_1d: Start", ai_dtllevel = 1)

        CALL field_from_matrix_1D_fem ( mo_fem, ai_field, ai_patch, ar_alpha, apr_matrix, ai_n_1 )

        CALL printlog("field_from_matrix_1d: End", ai_dtllevel = 1)

    end subroutine field_from_matrix_1d   
!----------------------------------------------------------------------------------------------
    subroutine field_to_matrix_2d ( ai_field, ai_patch, apr_matrix, ai_n_1, ai_n_2 )
        implicit none
        integer, intent(in)  :: ai_field
        integer, intent(in)  :: ai_patch
        integer, intent(in)  :: ai_n_1
        integer, intent(in)  :: ai_n_2
        real*8 , dimension(ai_n_1, ai_n_2), intent(out)  :: apr_matrix
        ! LOCAL VARIABLES

        CALL printlog("field_to_matrix_2d : Start", ai_dtllevel = 1)

        CALL field_to_matrix_2D_fem ( mo_fem, ai_field, ai_patch, apr_matrix, ai_n_1, ai_n_2 )

        CALL printlog("field_to_matrix_2d : End", ai_dtllevel = 1)

    end subroutine field_to_matrix_2d
!----------------------------------------------------------------------------------------------
    subroutine field_from_matrix_2d ( ai_field, ai_patch, ar_alpha, apr_matrix, ai_n_1, ai_n_2 )
        implicit none
        integer, intent(in)  :: ai_field
        integer, intent(in)  :: ai_patch
        REAL*8 , INTENT(IN)  :: ar_alpha
        integer, intent(in)  :: ai_n_1
        integer, intent(in)  :: ai_n_2
        real*8 , dimension(ai_n_1, ai_n_2), intent(in)  :: apr_matrix
        ! LOCAL VARIABLES

        CALL printlog("field_from_matrix_2d: Start", ai_dtllevel = 1)

        CALL field_from_matrix_2D_fem ( mo_fem, ai_field, ai_patch, ar_alpha, apr_matrix, ai_n_1, ai_n_2 )

        CALL printlog("field_from_matrix_2d: End", ai_dtllevel = 1)

    end subroutine field_from_matrix_2d
!----------------------------------------------------------------------------------------------
    subroutine field_to_matrix_3d ( ai_field, ai_patch, apr_matrix, ai_n_1, ai_n_2 , ai_n_3 )
        implicit none
        integer, intent(in)  :: ai_field
        integer, intent(in)  :: ai_patch
        integer, intent(in)  :: ai_n_1
        integer, intent(in)  :: ai_n_2
        integer, intent(in)  :: ai_n_3
        real*8 , dimension(ai_n_1, ai_n_2, ai_n_3), intent(out)  :: apr_matrix
        ! LOCAL VARIABLES

        CALL printlog("field_to_matrix_3d : Start", ai_dtllevel = 1)

        CALL field_to_matrix_3D_fem ( mo_fem, ai_field, ai_patch, apr_matrix, ai_n_1, ai_n_2 , ai_n_3 )

        CALL printlog("field_to_matrix_3d : End", ai_dtllevel = 1)

    end subroutine field_to_matrix_3d
!----------------------------------------------------------------------------------------------
    subroutine field_from_matrix_3d ( ai_field, ai_patch, ar_alpha, apr_matrix, ai_n_1, ai_n_2 , ai_n_3 )
        implicit none
        integer, intent(in)  :: ai_field
        integer, intent(in)  :: ai_patch
        REAL*8 , INTENT(IN)  :: ar_alpha
        integer, intent(in)  :: ai_n_1
        integer, intent(in)  :: ai_n_2
        integer, intent(in)  :: ai_n_3
        real*8 , dimension(ai_n_1, ai_n_2, ai_n_3), intent(in)  :: apr_matrix
        ! LOCAL VARIABLES

        CALL printlog("field_from_matrix_3d: Start", ai_dtllevel = 1)

        CALL field_from_matrix_3D_fem ( mo_fem, ai_field, ai_patch, ar_alpha, apr_matrix, ai_n_1, ai_n_2 , ai_n_3 )

        CALL printlog("field_from_matrix_3d: End", ai_dtllevel = 1)

    end subroutine field_from_matrix_3d
!!----------------------------------------------------------------------------------------------
!    subroutine pyfem_set_assl_position ( ai_space, apr_a, api_ia, api_ja, ai_nel, ai_nR )
!        implicit none
!        integer, intent(in)  :: ai_space
!        integer, intent(in)  :: ai_nR   !NUMBER OF ROWS
!        integer, intent(in)  :: ai_nel  !NUMBER OF NON ZERO ELTS
!        real*8 , dimension(ai_nel), intent(in)  :: apr_a
!        integer, dimension(ai_nR+1), intent(in)  :: api_ia
!        integer, dimension(ai_nel), intent(in)  :: api_ja
!
!        CALL printlog("pyfem_set_assl_position : Start", ai_dtllevel = 1)
!
!        CALL assl_set_a ( mo_ass % opo_bbox_sp(ai_space) % oo_asslx, apr_a, api_ia, api_ja)
!
!!        IF (ai_term==0) THEN
!!            CALL set_a ( mo_ass % opo_bbox_sp(ai_space) % opo_asslx (ai_x, ai_d), apr_a, api_ia, api_ja)
!!        END IF
!!        IF (ai_term==1) THEN
!!            CALL set_r ( mo_ass % opo_bbox_sp(ai_space) % opo_asslx (ai_x, ai_d), apr_a, api_ia, api_ja)
!!        END IF
!!        IF (ai_term==2) THEN
!!            CALL set_b ( mo_ass % opo_bbox_sp(ai_space) % opo_asslx (ai_x, ai_d), apr_a, api_ia, api_ja)
!!        END IF
!
!        CALL printlog("pyfem_set_assl_position : End", ai_dtllevel = 1)
!
!    end subroutine pyfem_set_assl_position
!    !---------------------------------------------------------------
!    subroutine pyfem_set_tensor_gbasis ( ai_space, apr_gbasis, ai_nderiv, ai_nnz)
!        implicit none
!        integer :: ai_space
!        INTEGER :: ai_nderiv
!        INTEGER :: ai_nnz
!        REAL(8), DIMENSION(0:ai_nderiv, 1:ai_nnz) :: apr_gbasis
!        ! LOCAL
!
!        CALL printlog("pyfem_set_tensor_gbasis : Start", ai_dtllevel = 1)
!
!        mo_ass % opo_bbox_sp(ai_space) % opr_gBasis (0:ai_nderiv, 1:ai_nnz) = apr_gBasis(0:ai_nderiv, 1:ai_nnz)
!
!        CALL printlog("pyfem_set_tensor_gbasis : End", ai_dtllevel = 1)
!
!    end subroutine pyfem_set_tensor_gbasis
!    !---------------------------------------------------------------
!    subroutine pyfem_set_tensor_gbasis_t ( ai_space, apr_gbasis, ai_nderiv, ai_nnz)
!        implicit none
!        integer :: ai_space
!        INTEGER :: ai_nderiv
!        INTEGER :: ai_nnz
!        REAL(8), DIMENSION(0:ai_nderiv, 1:ai_nnz) :: apr_gbasis
!        ! LOCAL
!
!        CALL printlog("pyfem_set_tensor_gbasis_t : Start", ai_dtllevel = 1)
!
!        mo_ass % opo_bbox_sp(ai_space) % opr_gBasis_t (0:ai_nderiv, 1:ai_nnz) = apr_gBasis(0:ai_nderiv, 1:ai_nnz)
!
!        CALL printlog("pyfem_set_tensor_gbasis_t : End", ai_dtllevel = 1)
!
!    end subroutine pyfem_set_tensor_gbasis_t
!    !---------------------------------------------------------------
!    subroutine pyfem_set_tensor_lbasis ( ai_space, apr_lbasis, ai_nderiv_code, ai_kg, ai_nel)
!        implicit none
!        integer :: ai_space
!        INTEGER :: ai_nderiv_code
!        INTEGER :: ai_kg
!        INTEGER :: ai_nel
!        REAL(8), DIMENSION(0:ai_nderiv_code, 1:ai_kg, 1:ai_nel) :: apr_lbasis
!        ! LOCAL
!
!        CALL printlog("pyfem_set_tensor_lbasis : Start", ai_dtllevel = 1)
!
!        mo_ass % opo_bbox_sp(ai_space) % opr_lBasis (0:ai_nderiv_code, 1:ai_kg, 1:ai_nel) =  &
!        apr_lbasis(0:ai_nderiv_code, 1:ai_kg, 1:ai_nel)
!
!        CALL printlog("pyfem_set_tensor_lbasis : End", ai_dtllevel = 1)
!
!    end subroutine pyfem_set_tensor_lbasis
!    !---------------------------------------------------------------
!    subroutine pyfem_set_tensor_lbasis_t ( ai_space, apr_lbasis, ai_nderiv_code, ai_kg, ai_nel)
!        implicit none
!        integer :: ai_space
!        INTEGER :: ai_nderiv_code
!        INTEGER :: ai_kg
!        INTEGER :: ai_nel
!        REAL(8), DIMENSION(0:ai_nderiv_code, 1:ai_kg, 1:ai_nel) :: apr_lbasis
!        ! LOCAL
!
!        CALL printlog("pyfem_set_tensor_lbasis_t : Start", ai_dtllevel = 1)
!
!        mo_ass % opo_bbox_sp(ai_space) % opr_lBasis_t (0:ai_nderiv_code, 1:ai_kg, 1:ai_nel) =  &
!        apr_lbasis(0:ai_nderiv_code, 1:ai_kg, 1:ai_nel)
!
!        CALL printlog("pyfem_set_tensor_lbasis_t : End", ai_dtllevel = 1)
!
!    end subroutine pyfem_set_tensor_lbasis_t
    !---------------------------------------------------------------
    subroutine pyfem_set_X_matrix ( ai_space, ai_patch, apr_X, ai_Rd, ai_nel, ai_nnz)
        implicit none
        integer :: ai_space
        integer :: ai_patch
        INTEGER :: ai_Rd
        INTEGER :: ai_nel
        INTEGER :: ai_nnz
        REAL(8), DIMENSION(1:ai_Rd, 1:ai_nel, 1:ai_nnz) :: apr_X
        ! LOCAL

        CALL printlog("pyfem_set_X_matrix : Start", ai_dtllevel = 1)

        mo_fem % opo_spaces (ai_space) % oo_mapping % opo_geo (ai_patch+1) % opr_X (1:ai_Rd, 1:ai_nel, 1:ai_nnz) = &
            apr_X (1:ai_Rd, 1:ai_nel, 1:ai_nnz)

        CALL printlog("pyfem_set_X_matrix : End", ai_dtllevel = 1)

    end subroutine pyfem_set_X_matrix
!----------------------------------------------------------------------------------------------
    subroutine createcsrmatrix ( ai_ref, ai_nC, api_ia, api_ja, ai_nR, ai_nel )
        implicit none
        integer, intent(in)  :: ai_ref  !THE MATRIX FROM WICH WE WANT TO ADD CONTRIBUTION
        integer, intent(in)  :: ai_nC   !NUMBER OF COLUMNS
        integer, intent(in)  :: ai_nR   !NUMBER OF ROWS
        integer*8, intent(in)  :: ai_nel  !NUMBER OF NON ZERO ELTS
        integer, dimension(ai_nR+1), intent(in)  :: api_ia
        integer*8, dimension(ai_nel), intent(in)  :: api_ja
        ! LOCAL VARIABLES

        CALL printlog("createcsrmatrix : Begin", ai_dtllevel = 1)

        CALL createcsrmatrix_fem ( mo_fem, ai_ref, api_ia, api_ja, ai_nR, ai_nC, ai_nel )

        CALL printlog("createcsrmatrix : End", ai_dtllevel = 1)

    end subroutine createcsrmatrix
!----------------------------------------------------------------------------------------------
    subroutine setcsrmatrix ( ai_ref, apr_a, ai_nel )
        implicit none
        integer, intent(in)  :: ai_ref  !THE MATRIX FROM WICH WE WANT TO ADD CONTRIBUTION
        integer*8, intent(in)  :: ai_nel  !NUMBER OF NON ZERO ELTS
        real*8 , dimension(ai_nel), intent(in)  :: apr_a
        ! LOCAL VARIABLES

!#ifdef _TRACE
!        CALL printlog("setcsrmatrix : Begin", ai_dtllevel = 1)
!#endif

        CALL setcsrmatrix_fem ( mo_fem, ai_ref, apr_a, ai_nel )

!#ifdef _TRACE
!        CALL printlog("setcsrmatrix : End", ai_dtllevel = 1)
!#endif

    end subroutine setcsrmatrix
!----------------------------------------------------------------------------------------------
    subroutine getparamcsrmatrix ( ai_ref, ai_nR, ai_nC, ai_nel )
        implicit none
        integer, intent(in)  :: ai_ref  !THE MATRIX FROM WICH WE WANT TO ADD CONTRIBUTION
        integer, intent(out)  :: ai_nR   !NUMBER OF ROWS
        integer, intent(out)  :: ai_nC   !NUMBER OF COLUMNS
        integer*8, intent(out)  :: ai_nel  !NUMBER OF NON ZERO ELTS
        ! LOCAL
        INTEGER :: IERROR

        CALL SPM_GetnR (ai_ref, ai_nR , IERROR)
        CALL SPM_GetnC (ai_ref, ai_nC , IERROR)
        CALL SPM_Getnnz(ai_ref, ai_nel, IERROR)

    end subroutine getparamcsrmatrix
!----------------------------------------------------------------------------------------------
    subroutine getarrayparamcsrmatrix ( ai_ref, apr_a, api_ia, api_ja, ai_nR, ai_nel )
        implicit none
        integer, intent(in)  :: ai_ref  !THE MATRIX FROM WICH WE WANT TO ADD CONTRIBUTION
        integer, intent(in)  :: ai_nR   !NUMBER OF ROWS
        integer*8, intent(in)  :: ai_nel  !NUMBER OF NON ZERO ELTS
        real*8 , dimension(ai_nel), intent(out)  :: apr_a
        integer, dimension(ai_nR+1), intent(out)  :: api_ia
        integer*8, dimension(ai_nel), intent(out)  :: api_ja
        ! LOCAL VARIABLES

!#ifdef _TRACE
!        CALL printlog("getarrayparamcsrmatrix : Begin", ai_dtllevel = 1)
!#endif

        CALL getarrayparamcsrmatrix_fem ( mo_fem, ai_ref, apr_a, api_ia, api_ja, ai_nR, ai_nel )                

!#ifdef _TRACE
!        CALL printlog("getarrayparamcsrmatrix : End", ai_dtllevel = 1)
!#endif        
    end subroutine getarrayparamcsrmatrix    
!----------------------------------------------------------------------------------------------
    subroutine pyfem_getsizeijv( ai_ref, ai_size )
        implicit none
        integer, intent(in)  :: ai_ref  !THE MATRIX FROM WICH WE WANT TO ADD CONTRIBUTION
        integer*8, intent(out)  :: ai_size
        ! LOCAL
        INTEGER :: IERR

        CALL SPM_GetSIZEIJV(ai_ref, ai_size, IERR)

    end subroutine pyfem_getsizeijv
!----------------------------------------------------------------------------------------------
    subroutine pyfem_getdataijv( ai_ref, api_i, api_j, apr_v, ai_size )
        implicit none
        integer, intent(in)  :: ai_ref  !THE MATRIX FROM WICH WE WANT TO ADD CONTRIBUTION
        integer*8, intent(in)  :: ai_size  !NUMBER OF NON ZERO ELTS
        integer*8, dimension(ai_size), intent(out)      :: api_i
        integer*8, dimension(ai_size), intent(out)      :: api_j
        real*8 , dimension(ai_size), intent(out)        :: apr_v
        ! LOCAL VARIABLES
        integer  :: IERR
!#ifdef _TRACE
!        CALL printlog("pyfem_getdataijv: Begin", ai_dtllevel = 1)
!#endif

        CALL SPM_GetDATAIJV(ai_ref, ai_size, api_i, api_j, apr_v, IERR)

!#ifdef _TRACE
!        CALL printlog("pyfem_getdataijv: End", ai_dtllevel = 1)
!#endif        
    end subroutine pyfem_getdataijv
end module pyfem

