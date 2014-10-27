!     
! File:   assembly_def.F90
! Author: root
!
! Created on January 2, 2012, 2:14 PM
!

module assembly_def
    use used_precision
    use bbox_def
!    use geometries_def
!    use connectivities_def
!    use grids_def
!    use fem_def
    implicit none

    TYPE, PUBLIC :: METRIC_INFO
        !> needed to compute the jacobian
        REAL(wp), DIMENSION(:,:,:), POINTER :: opr_points
        !> the jacobian for each point of the grid
        REAL(wp), DIMENSION (:), POINTER :: opr_jacobians
        !> the inverse of the jacobian matrix
        REAL(wp), DIMENSION (:,:,:), POINTER :: opr_invJacobian
        !> this is the non normalized inverse jacobian:
        !> we do not divide by the determinant
        REAL(wp), DIMENSION (:,:,:), POINTER :: opr_invJacobian_nn
    END TYPE METRIC_INFO

    TYPE, PUBLIC :: PHYSICAL_BASIS 
        REAL(wp), DIMENSION(:,:)  , POINTER :: opr_B
        REAL(wp), DIMENSION(:,:,:), POINTER :: opr_gradB
        REAL(wp), DIMENSION(:,:,:), POINTER :: opr_curlB
        REAL(wp), DIMENSION(:,:,:), POINTER :: opr_HessianB
    END TYPE PHYSICAL_BASIS

    TYPE, PUBLIC :: ASSEMBLY
        !>
        LOGICAL :: ol_initialized = .FALSE.
        !>
        INTEGER :: oi_maxnpts
        !>
        INTEGER :: oi_space_maxnen
        INTEGER :: oi_mapping_maxnen
        !>
        INTEGER :: oi_maxndof
        !> TRUE IF WE NEED TO COMPUTE POSITIONS AND BASIS FUNCTIONS AT THE GRID POINTS
        !> IT IS FALSE FOR POISSON BRACKET, TRUE OTHERWISE
        !> DEFAULT VALUE : FALSE
        LOGICAL :: ol_assembly_points
        LOGICAL :: ol_assembly_basis

        !> the highest order of derivatives among all used operators
        !> \todo for each operator, we must add a subroutine that gives the corresponding nderiv
        INTEGER :: oi_matrix_nderiv
        !> the highest order of derivatives among all mappings
        !> generally it is equal to 1, but we may need, in the PIC code, to compute the 2nd derivatives
        INTEGER :: oi_mapping_nderiv

        !>
        REAL(wp), DIMENSION (:,:,:), POINTER :: opr_Matrix_elt
        REAL(wp), DIMENSION (:), POINTER :: opr_Matrix_eltLine
        !>
        REAL(wp), DIMENSION (:,:), POINTER :: opr_Projection_elt
        !>
        REAL(wp), DIMENSION (:,:,:,:), POINTER :: opr_Fieldh_elt
        !>
        REAL(wp), DIMENSION(:), POINTER :: opr_Norm_elt

    ! ...
        !> for each mapping, we must have an approriate black-box
        TYPE(BBOX), DIMENSION(:), POINTER :: opo_bbox_mp
        !> for each space, we must have an approriate black-box
        TYPE(BBOX), DIMENSION(:), POINTER :: opo_bbox_sp
    ! ...

    ! ...
        !> points and jacobians will be stored here
        TYPE(METRIC_INFO), DIMENSION(:), POINTER :: opo_info_sp
        !> points and jacobians will be stored here
        TYPE(METRIC_INFO), DIMENSION(:), POINTER :: opo_info_mp
        !> points and jacobians will be stored here
        TYPE(METRIC_INFO), DIMENSION(:), POINTER :: opo_info_gr
    ! ...

    ! ...
        !> points and jacobians will be stored here
        TYPE(PHYSICAL_BASIS), DIMENSION(:), POINTER :: opo_pBasis
    ! ...

        !> Patchs to assembly
        INTEGER, DIMENSION(:), POINTER :: opi_patchs_toassembly
        !> Elements to assembly
        INTEGER, DIMENSION(:), POINTER :: opi_elts_toassembly
        !> Operators to assembly
        INTEGER, DIMENSION(:), POINTER :: opi_operators_toassembly
        !> Matrices to assembly
        INTEGER, DIMENSION(:), POINTER :: opi_matrices_toassembly        
        !> Fields to assembly
        INTEGER, DIMENSION(:), POINTER :: opi_fields_toassembly
        !> Norms to assembly
        INTEGER, DIMENSION(:), POINTER :: opi_norms_toassembly
        !> Spaces to assembly
        INTEGER, DIMENSION(:), POINTER :: opi_spaces_toassembly

        !> tmp : is used to store temp the ids because we need to assemble a temp terms
        !> Patchs to assembly
        INTEGER, DIMENSION(:), POINTER :: opi_patchs_toassembly_tmp
        !> Elements to assembly
        INTEGER, DIMENSION(:), POINTER :: opi_elts_toassembly_tmp
        !> Operators to assembly
        INTEGER, DIMENSION(:), POINTER :: opi_operators_toassembly_tmp
        !> Matrices to assembly
        INTEGER, DIMENSION(:), POINTER :: opi_matrices_toassembly_tmp        
        !> Fields to assembly
        INTEGER, DIMENSION(:), POINTER :: opi_fields_toassembly_tmp
        !> Norms to assembly
        INTEGER, DIMENSION(:), POINTER :: opi_norms_toassembly_tmp
        !> Spaces to assembly
        INTEGER, DIMENSION(:), POINTER :: opi_spaces_toassembly_tmp

    END TYPE ASSEMBLY

end module assembly_def


