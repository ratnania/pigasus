!
! File:   fem_def.F90
! Author: root
!
! Created on January 10, 2012, 6:03 PM
!

module fem_def
    USE USEd_precision
    USE tracelog_module
    USE grids_def
    USE grids
    USE geometries_def
    USE bbox_def
    USE connectivities_def
    USE SPM 
    implicit none
!#ifdef _MURGE
!! will be replaced during compilation by replaceCOEF.sh
!INCLUDE "mpif.h"
!INCLUDE "murge.inc"
!#endif

    TYPE, PUBLIC :: SOLVER 
        !> maximum iteration of the iterative solver
        INTEGER  :: oi_maxiter
        !> number of iterations after convergence of the iterative solver
        INTEGER  :: oi_niter
        !> absolute error
        REAL(8)  :: or_atol
        !> relative error
        REAL(8)  :: or_rtol
        !> precision of the iterative solver
        REAL(8)  :: or_eps
        !> errors after each iteration of the iterative solver
        REAL(8), DIMENSION(:), POINTER :: opr_err 
        REAL(8) :: or_resnorm
        !> True if opr_err is allocated
        LOGICAL  :: ol_allocatedErr
        LOGICAL  :: ol_allocated
        !> coeff of the current field
        REAL(8), DIMENSION (:), POINTER :: opr_U ! solution
        REAL(8), DIMENSION (:), POINTER :: opr_B ! rhs
    END TYPE SOLVER 

    TYPE, PUBLIC :: GRAPH
        INTEGER :: oi_nR
        INTEGER :: oi_nC
        INTEGER*8 :: oi_nnz
!        INTEGER, DIMENSION (:), POINTER :: opi_ind_col ! COO
!        INTEGER, DIMENSION (:), POINTER :: opi_ind_row ! COO

        INTEGER, DIMENSION(:), POINTER :: opi_ia   !CSC/CSR format
        INTEGER*8, DIMENSION(:), POINTER :: opi_ja !CSC/CSR format
    END TYPE GRAPH 

    TYPE, PUBLIC :: OPER
        !>
        INTEGER :: oi_type
        !>
        INTEGER, DIMENSION (:), POINTER :: opi_matrices_toassembly
        !>
        INTEGER, DIMENSION (:), POINTER :: opi_matrices_toassembly_tmp
        !>
        REAL(wp), DIMENSION (:), POINTER :: opr_scale
    END TYPE OPER

    TYPE, PUBLIC :: FIELD
        INTEGER :: oi_size
        !> coeff of the current field
        real(wp), DIMENSION (:), POINTER :: opr_c
        real(wp), DIMENSION (:,:,:,:), POINTER :: opr_val
        !> for each index we can get the number of access
        INTEGER, DIMENSION (:), POINTER :: opi_counter 
       
    END TYPE FIELD

    TYPE, PUBLIC :: NORM
        real(wp), DIMENSION (:,:), POINTER :: opr_values
    END TYPE NORM
    
    TYPE, PUBLIC :: SPACE
        !> Quad/triangles/...
        INTEGER :: oi_TYPE
        !> tensorlevel
        INTEGER :: oi_tensorlevel
        !> maximum polynomial degree
        INTEGER :: oi_maxp
        INTEGER :: oi_maxnen
        !> the corresponding mapping if we USE the IGA approach
        TYPE(GEOMETRIES) :: oo_mapping
        !> the corresponding connectivity arrays
        TYPE(CONNECTIVITY) :: oo_con
    END TYPE SPACE

    TYPE, PUBLIC :: MAPPING
        !> maximum polynomial degree
        INTEGER :: oi_maxp
        INTEGER :: oi_maxnen
        !> the corresponding mapping if we USE the IGA approach
        TYPE(GEOMETRIES) :: oo_mapping
    END TYPE MAPPING

    TYPE, PUBLIC :: METRIC
        !> an array of points of the form Points[i, 0:Rd,0:Rd], i is the i^th point in the patch id and element elt
        !>  Points[id, elt,i,0,:] is the physical point
        !>  Points[id, elt,i,1,:] is the x-derivative
        !>  Points[id, elt,i,2,:] is the y-derivative
        !>  Points[id, elt,i,3,:] is the z-derivative
        real(wp), DIMENSION(:,:,:,:,:), POINTER :: opr_points
    END TYPE METRIC
   
    TYPE, PUBLIC :: FEM        
    ! ...
        !>
        LOGICAL :: ol_initialized = .FALSE.    
        !> DOMAIN'S DIMENSION
        !> FOR EACH GRIDS WE CAN HAVE A SPECIFIC DIMENSION; ex for THE FAST APPROACH
        INTEGER, DIMENSION(:), POINTER :: opi_dim
        !> generally Rd=dim, unless we treate a 2D-curve or 3D-curve/surface domain
        INTEGER, DIMENSION(:), POINTER :: opi_Rd
        !>
        INTEGER :: oi_nSolvers                
        !>
        INTEGER :: oi_nGraphs        
        !>
        INTEGER :: oi_nMatrices
        !>
        INTEGER :: oi_nOperators        
        !>
        INTEGER :: oi_nFields
        !>
        INTEGER :: oi_nNorms       
        !> number of involved spaces
        INTEGER :: oi_nspaces
        !> number of involved mappings
        INTEGER :: oi_nmappings
        !> number of involved metrics
        INTEGER :: oi_nmetrics
        !> number of involved grids
        INTEGER :: oi_nGrids
        !> the maximum number of patchs over all grids
        INTEGER :: oi_maxnpatchs
        !> the maximum number of elements over all grids
        INTEGER :: oi_maxnelts
        !> the maximum of params for all matrices and spaces
        !> \todo must be done for each space, to reduce memory
        INTEGER :: oi_maxnparams_operators
        INTEGER :: oi_maxnparams_fields
        INTEGER :: oi_maxnparams_norms
        !> The maximum derivative order for fields operators
        INTEGER :: oi_maxder_field_op
        !> 
        INTEGER :: oi_maxaddto
        !> 
        INTEGER :: oi_maxnen        
    ! ...

    ! ...
        !> we can have different GRIDS (ex. fot the fast approach)
        TYPE(GRIDS_DATA), DIMENSION(:), POINTER :: opo_grids
        !> INFORMATION OF THE i^th GRIDS
        INTEGER, DIMENSION (:,:), POINTER :: opi_InfoGrids
    ! ...

    ! ...
        !> INFORMATION OF THE i^th MAPPING
        INTEGER, DIMENSION (:,:), POINTER :: opi_InfoMapping
    ! ...

    ! ...
        !>
        TYPE(MAPPING), DIMENSION(:), POINTER :: opo_mappings
    ! ...

    ! ...
        !>
        TYPE(METRIC), DIMENSION(:), POINTER :: opo_metrics
    ! ...

    ! ...
        !>
        TYPE(SPACE), DIMENSION(:), POINTER :: opo_spaces
        !> INFORMATION OF THE i^th SPACE
        INTEGER, DIMENSION (:,:), POINTER :: opi_InfoSpace
    ! ...

    ! ...
        !> OPERATOR : THIS IS THE UNKNOWN DIFFERENTIAL OPERATOR
        TYPE(OPER), DIMENSION (:), POINTER :: opo_Op
        INTEGER, DIMENSION (:,:), POINTER :: opi_InfoOperator
    ! ...

    ! ...
        !> FIELDS : THIS IS THE UNKNOWN OF OUR PDE
        TYPE(FIELD), DIMENSION (:), POINTER :: opo_F
        INTEGER, DIMENSION (:,:), POINTER :: opi_InfoField
    ! ...
        
    ! ...
        INTEGER, DIMENSION (:,:), POINTER :: opi_InfoMatrix
    ! ...

    ! ...
        !> CONNECTIVITY GRAPH 
        INTEGER, DIMENSION (:,:), POINTER :: opi_InfoGraph
        TYPE(GRAPH), DIMENSION (:), POINTER :: opo_graph
    ! ...

    ! ...
        !> for each grids and patch we get an info
        INTEGER, DIMENSION (:,:,:), POINTER :: opi_InfoPatch
    ! ...

    ! ...
        !> NORMS
        INTEGER, DIMENSION (:,:), POINTER :: opi_InfoNorm
        TYPE(NORM), DIMENSION (:), POINTER :: opo_N
    ! ...

    ! ...
        INTEGER, DIMENSION (:,:), POINTER :: opi_InfoSolver
        TYPE(SOLVER), DIMENSION (:), POINTER :: opo_Solver
    ! ...
    END TYPE FEM

!    !> mpl_AssemblyMatrix(i) = .TRUE. IF WE WANT TO ASSEMBLY THE i^th MATRIX
!    logical, DIMENSION (:), allocatable :: mpl_AssemblyMatrix
!    !> mpl_AssemblyMatrix WILL BE USED ONLY IF ml_AssemblyMatrix IS TRUE
!    logical :: ml_AssemblyMatrix


! ***************************************************
!
! ***************************************************
    !> EACH MATRIX MAY HAVE DIFFERENT PARAMS, DEPENDING ON THE TYPE OF THE MATRIX
    !> Mass                 ->      1
    !> Stiffness            ->      4
    !> Advection            ->      2
    !> Second-Deriv         ->      3**2 for 2D and 6**2 in 3D
    INTEGER, PARAMETER :: NPARAM_MATRIX = 9 
! ***************************************************

! ***************************************************
!   INFOSOLVER
! ***************************************************
    INTEGER, PARAMETER :: NPARAM_INFOSOLVER        = 4 
    !> IF THE SOLVER IS ASSOCIATED TO A DEFINED MATRIX
    !> WE SET HERE ITS ID 
    INTEGER, PARAMETER :: INFOSOLVER_MATRIX        = 1
    !> IF THE SOLVER IS ASSOCIATED TO A DEFINED OPERATOR (RESIDUAL) 
    !> WE SET HERE ITS ID 
    INTEGER, PARAMETER :: INFOSOLVER_RESIDUAL      = 2
    !> WHICH SOLVER TO USE 
    INTEGER, PARAMETER :: INFOSOLVER_SOLVER        = 3
    !> 
    INTEGER, PARAMETER :: INFOSOLVER_MAXITER       = 4
!    !> 
!    INTEGER, PARAMETER :: INFOSOLVER_        = 
! ***************************************************

! ***************************************************
!       SOLVERS AND LIN-LIBRARIES IDs
! ***************************************************
! ... solvers
    INTEGER, PARAMETER :: CG       = 1
    INTEGER, PARAMETER :: BCG      = 2
    INTEGER, PARAMETER :: DBCG     = 3
    INTEGER, PARAMETER :: CNRG     = 4
    INTEGER, PARAMETER :: BCGSTAB  = 5
    INTEGER, PARAMETER :: TFQMR    = 6
    INTEGER, PARAMETER :: FOM      = 7
    INTEGER, PARAMETER :: GMRES    = 8
    INTEGER, PARAMETER :: FGMRES   = 9
    INTEGER, PARAMETER :: DQGMRES  = 10

! ... libraries
    INTEGER, PARAMETER :: BASIC_SLV    = 0
    INTEGER, PARAMETER :: PASTIX_SLV   = 1
    INTEGER, PARAMETER :: CUSP_SLV     = 2
! ***************************************************

! ***************************************************
!   INFOGRAPH
! ***************************************************
    INTEGER, PARAMETER :: NPARAM_INFOGRAPH         = 2 

    INTEGER, PARAMETER :: INFOGRAPH_SPACE_1        = 1
    INTEGER, PARAMETER :: INFOGRAPH_SPACE_2        = 2
! ***************************************************

! ***************************************************
!   INFOMATRIX
! ***************************************************
    !> THIS IS THE SIZE OF THE ARRAY INFO_MATRIX
    !> EACH MATRIX MUST HAVE A TYPE, A COUPLE OF FIELDS (i,j), ... AND OTHER RELEVANT PARAMS
    !> FOR THE MOMENT THE CORRESPONDANT VALUE FOR NPARAM_INFOMATRIX IS RESERVED
    INTEGER, PARAMETER :: NPARAM_INFOMATRIX         = 4 

    INTEGER, PARAMETER :: INFOMATRIX_TYPE           = 1
    !> the corresponding graph
    INTEGER, PARAMETER :: INFOMATRIX_GRAPH          = 2 
    !> 
    INTEGER, PARAMETER :: INFOMATRIX_TOASSEMBLY     = 3
    !> If true, we use IJV assembly
    !> the assembly process will generate tuple (i,j,v)
    !> the couples (i,j) are unsorted and redandont
    INTEGER, PARAMETER :: INFOMATRIX_IJVASSEMBLY    = 4 
! ***************************************************

! ***************************************************
!
! ***************************************************
    !>	TYPE MATRICES
    INTEGER, PARAMETER :: GENERIC_MATRIX            = 0  
    INTEGER, PARAMETER :: BLOCK_MATRIX              = 1
! ***************************************************

! ***************************************************
!
! ***************************************************
    !>	DIMENSION OF COVERED OPERATORS
    !>  FOR THE MOMENT WE ONLY HANDLE 2-OPERATORS : (MASS,STIFFNESS,ADVECTION,ROTATIONAL,R_MATRIX)
    !>  AND 3-OPERATORS (POISSON_BRACKET)
    INTEGER, PARAMETER :: MAX_NDIM_OPERATOR         = 3
! ***************************************************

! ***************************************************
!   INFOOPERATOR
! ***************************************************
    !> THIS IS THE SIZE OF THE ARRAY INFO_MATRIX
    !> EACH MATRIX MUST HAVE A TYPE, A COUPLE OF FIELDS (i,j), ... AND OTHER RELEVANT PARAMS
    !> FOR THE MOMENT THE CORRESPONDANT VALUE FOR NPARAM_INFOMATRIX IS RESERVED
    INTEGER, PARAMETER :: NPARAM_INFOOPERATOR         = 9

    INTEGER, PARAMETER :: INFOOPERATOR_TYPE           = 1
    INTEGER, PARAMETER :: INFOOPERATOR_SPACE_1        = 2
    INTEGER, PARAMETER :: INFOOPERATOR_SPACE_2        = 3
    !> LOCAL ID FOR THE GRIDS/SPACE
    INTEGER, PARAMETER :: INFOOPERATOR_LOCID          = 4 
    !> THE NUMBER OF PARAM
    INTEGER, PARAMETER :: INFOOPERATOR_NPARAM         = 5
    !> IF WE WANT TO EVALUATE BASIS IN THE PARAM DOMAIN (i.e. GRADIANT, DERIVATIVES, ...)
    INTEGER, PARAMETER :: INFOOPERATOR_PARAMEVAL      = 6
    !> if we want to transpose the matrix during the computations
    INTEGER, PARAMETER :: INFOOPERATOR_TRANSPOSE      = 7
    !> the higher order of derivative
    INTEGER, PARAMETER :: INFOOPERATOR_NDERIV         = 8 
    !> 
    INTEGER, PARAMETER :: INFOOPERATOR_TOASSEMBLY     = 9 
! ***************************************************

! ***************************************************
!
! ***************************************************
    !>	OPERATORS TYPES
    INTEGER, PARAMETER :: MASS                      = 0 
    INTEGER, PARAMETER :: STIFFNESS                 = 1
    INTEGER, PARAMETER :: ADVECTION                 = 2
    INTEGER, PARAMETER :: SECOND_DERIV              = 3
!    INTEGER, PARAMETER :: ROTATIONAL                = 6
!    INTEGER, PARAMETER :: ROTATIONAL_SCALAR         = 7
! ***************************************************

! ***************************************************
!   FIELDS OPERATORS
! ***************************************************
    INTEGER, PARAMETER :: IDENTITY                  = 1
    INTEGER, PARAMETER :: GRAD                      = 2
    INTEGER, PARAMETER :: CURL                      = 3
    INTEGER, PARAMETER :: IDENTITY_VECT             = 4
    INTEGER, PARAMETER :: SECOND_DERIV_FIELD        = 5
    INTEGER, PARAMETER :: GRAD_S                    = 6
    INTEGER, PARAMETER :: SECOND_DERIV_S_FIELD      = 8
    INTEGER, PARAMETER :: HESSIAN_FIELD             = 9
! ***************************************************

! ***************************************************
    !>	TYPE PROJECTORS
    !>	FOR THE MOMENT WE TREATE :
    !>		L2 Projection
    INTEGER, PARAMETER :: PROJECTION_L2 = 1
    INTEGER, PARAMETER :: FIELD_OPERATOR = 2
!    INTEGER, PARAMETER :: PROJECTION_L2_AXI = 2
!    INTEGER, PARAMETER :: ROJECTION_H1 = 3
!    INTEGER, PARAMETER :: PROJECTION_DIV = 4
! ***************************************************

! ***************************************************
!   INFOFIELD
! ***************************************************
    INTEGER, PARAMETER :: NPARAM_INFOFIELD = 10
    !> the DIMENSION of the field : scalar, vectorial: 2D/3D...
    INTEGER, PARAMETER :: INFOFIELD_NDOF  = 1
    !> the DIMENSION of the related discrete space
    INTEGER, PARAMETER :: INFOFIELD_SIZE  = 2
    !> in [4] we store the id of the corresponding space
    INTEGER, PARAMETER :: INFOFIELD_SPACE  = 3
    !> LOCAL ID FOR THE GRIDS/SPACE
    INTEGER, PARAMETER :: INFOFIELD_LOCID  = 4
    !> = 1 IF WE WANT TO ASSEMBLY THE CURRENT FIELD
    INTEGER, PARAMETER :: INFOFIELD_TOASSEMBLY  = 5
    !> PROJECTORS TYPE
    INTEGER, PARAMETER :: INFOFIELD_TYPE  = 6
    !> OPERATOR FIELD : IDENTITY, GRAD, CURL, ...
    INTEGER, PARAMETER :: INFOFIELD_OPERATOR  = 7
    !> OPERANDE FOR THE OPERATOR
    INTEGER, PARAMETER :: INFOFIELD_OPERANDE  = 8
    !> IF WE WANT TO EVALUATE BASIS IN THE PARAM DOMAIN (i.e. GRADIANT, DERIVATIVES, ...)
    INTEGER, PARAMETER :: INFOFIELD_PARAMEVAL  = 9
    !> NPARAM FOR FUNC EVAL
    INTEGER, PARAMETER :: INFOFIELD_NPARAM  = 10
! ***************************************************

! ***************************************************
!   INFOSPACE
! ***************************************************
    INTEGER, PARAMETER :: NPARAM_INFOSPACE = 13
    !> IF WE USE AN EXTERIOR MAPPING WITH THIS SPACE
    INTEGER, PARAMETER :: INFOSPACE_EXTMAPPING = 1
    !> ID OF THE EXTERIOR MAPPING IF INFOSPACE_EXTMAPPING == 1
    INTEGER, PARAMETER :: INFOSPACE_MAPPING = 2
    !> IF WE USE A DIFFERENT SPACE :  TRIANGLES, SIMPLEX, ...
    INTEGER, PARAMETER :: INFOSPACE_TYPE = 3
    !> IF WE USE A TENSOR PRODUCT TO EVALUATE BASIS FUNCTIONS : B-splines/NURBS
    !> = 1 : YES
    !> = 0 : NO
    INTEGER, PARAMETER :: INFOSPACE_TENSOR = 4
    !> THE ID OF THE CORRESPONDING GRIDS
    INTEGER, PARAMETER :: INFOSPACE_GRIDS = 5
    !> THIS ENABLE TO STORE DATA ONLY WHEN USING TENSOR DEFINITION
    INTEGER, PARAMETER :: INFOSPACE_STOREDDATA = 6
    !> if it is a scalar/vectorial field space
    !> this is in fact, the maximum of dof among all spaces related to the actual grids
    INTEGER, PARAMETER :: INFOSPACE_NDOF = 8
    !> the size (DIMENSION) of the space
    INTEGER, PARAMETER :: INFOSPACE_SIZE = 9
    !> the max of nen over all patchs
    INTEGER, PARAMETER :: INFOSPACE_MAXNEN = 10
    !> IF the space is composed
    INTEGER, PARAMETER :: INFOSPACE_COMPOSED = 11
    !> THE number of spaces USEd to compose the current one
    INTEGER, PARAMETER :: INFOSPACE_COMPOSED_NSP = 12
    !> SET TO 1 IF WE WANT TO ASSEMBLE THE CURRENT SPACE
    INTEGER, PARAMETER :: INFOSPACE_TOASSEMBLY = 13
!    INTEGER, PARAMETER :: INFOSPACE_ =
! ***************************************************

! ***************************************************
!   INFOMAPPING
! ***************************************************
    INTEGER, PARAMETER :: NPARAM_INFOMAPPING = 3
    !>
    INTEGER, PARAMETER :: INFOMAPPING_TENSOR = 1
    !> THIS ENABLE TO STORE DATA ONLY WHEN USING TENSOR DEFINITION
    INTEGER, PARAMETER :: INFOMAPPING_STOREDDATA = 2
    !> THIS THE ID OF THE CORRESPONDING SPACE
    INTEGER, PARAMETER :: INFOMAPPING_SPACE = 3
    !>
!    INTEGER, PARAMETER :: INFOMAPPING_ =
! ***************************************************

! ***************************************************
!   INFOPATCH
! ***************************************************
    INTEGER, PARAMETER :: NPARAM_INFOPATCH = 6
    !> NUMBER OF ELEMENTS IN THE CURRENT PATCH
    INTEGER, PARAMETER :: INFOPATCH_NEL  = 1
    !> MAXIMUM NUMBER OF POINTS IN THE GRID OF THE CURRENT PATCH
    INTEGER, PARAMETER :: INFOPATCH_MAXNPTS = 2
    !> = 1 <=> is(TENSOR) == TRUE
    INTEGER, PARAMETER :: INFOPATCH_TENSOR = 3
    !> MAXIMUM NUMBER OF POINTS IN THE GRID OF THE CURRENT PATCH AMONG ALL DIRECTIONS
    INTEGER, PARAMETER :: INFOPATCH_DIRMAXNPTS = 4
    !> SET TO 1 IF WE WANT TO ASSEMBLE THE CURRENT PATCH
    INTEGER, PARAMETER :: INFOPATCH_TOASSEMBLY = 5
    !>
!    INTEGER, PARAMETER :: INFOPATCH_ =
! ***************************************************

! ***************************************************
!   INFOGRIDS
! ***************************************************
    INTEGER, PARAMETER :: NPARAM_INFOGRIDS = 9
    !> NUMBER OF PATCHS IN THE CURRENT GRIDS
    INTEGER, PARAMETER :: INFOGRIDS_NPATCHS  = 1
    !> DIMENSION OF THE DOMAIN : 1D, 2D, 3D, ...
    INTEGER, PARAMETER :: INFOGRIDS_DIM = 2
    INTEGER, PARAMETER :: INFOGRIDS_RD = 3
    !> DIMENSION OF FIELDS : SCALAR/VECTOR...
    INTEGER, PARAMETER :: INFOGRIDS_DOF = 4
    !> THE NUMBER OF FIELDS OVER THE CURRENT SPACE
    INTEGER, PARAMETER :: INFOGRIDS_NFIELDS = 5
    !> THE NUMBER OF MATRICES OVER THE CURRENT SPACE
    INTEGER, PARAMETER :: INFOGRIDS_NMATRICES = 6
    !> USE OF THE METRIC IF INFOGRIDS_METRIC == 1
    INTEGER, PARAMETER :: INFOGRIDS_USEMETRIC = 7
    !> ID OF THE METRIC IF INFOGRIDS_USEMETRIC == 1
    INTEGER, PARAMETER :: INFOGRIDS_METRIC_ID = 8
    !> THE NUMBER OF NORMS OVER THE CURRENT SPACE
    INTEGER, PARAMETER :: INFOGRIDS_NNORMS = 9
    !>
!    INTEGER, PARAMETER :: INFOGRIDS_ =
! ***************************************************

! ***************************************************
    !>	TYPE NORMS
    !>	FOR THE MOMENT WE TREATE :
    !>		L2 NORM
    INTEGER, PARAMETER :: NORM_L2 = 1
    INTEGER, PARAMETER :: NORM_H1 = 2
! ***************************************************

! ***************************************************
!   INFONORM
! ***************************************************
    INTEGER, PARAMETER :: NPARAM_INFONORM = 6
    !> NORM TYPE : L2, H1, ...
    INTEGER, PARAMETER :: INFONORM_TYPE  = 1
    !> the related field
    INTEGER, PARAMETER :: INFONORM_FIELD  = 2
    !> if we want to compute the current norm
    INTEGER, PARAMETER :: INFONORM_TOASSEMBLY  = 3
    !> LOCAL ID FOR THE GRIDS/SPACE
    INTEGER, PARAMETER :: INFONORM_LOCID  = 4
! ***************************************************

END module fem_def


