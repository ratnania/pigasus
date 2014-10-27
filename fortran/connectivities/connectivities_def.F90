!     
! File:   connectivities_def.F90
! Author: root
!
! Created on December 20, 2011, 12:10 PM
!

module connectivities_def
    implicit none

    private

    TYPE, PUBLIC :: PATCH_CONNECTIVITY
        !> local ids for each patch
        INTEGER, DIMENSION(:), POINTER :: opi_ID_loc
    END TYPE PATCH_CONNECTIVITY

    TYPE, PUBLIC :: CONNECTIVITY
        INTEGER :: oi_npatch    
        INTEGER :: oi_nelt
        INTEGER :: oi_N
        !> this is the local connectivity array
        !> opi_IEN (id_p, b, elt)
        !> id_p : the id of the current patch
        !> b : the local num of the basis
        !> elt : the num of the element in the current patch
        INTEGER, DIMENSION(:,:,:), POINTER :: opi_IEN
        !> global ID 
        INTEGER, DIMENSION(:), POINTER :: opi_ID
        !> local ids for each patch
        TYPE(PATCH_CONNECTIVITY), DIMENSION(:), POINTER :: opi_info
        !> LM = ID ( IEN )
        INTEGER, DIMENSION(:,:,:), POINTER :: opi_LM
        !> nen : for each patch we give the number of element basis
        INTEGER, DIMENSION(:), POINTER :: opi_nen
        !> for each patch we give the correspondance
        !> between the real elt id and the loop-id
        !> this is to handle boundary grids, but can also be used
        !> in the case of local refinement, or to parallelize the assembling process
        !> when a local mesh has a lot of GL points (compared to others elts)
        INTEGER, DIMENSION(:,:), POINTER :: opi_real_elts
    END TYPE CONNECTIVITY

    TYPE, PUBLIC :: CONNECTIVITIES
        INTEGER :: oi_n
        TYPE(CONNECTIVITY), DIMENSION(:), POINTER :: opo_con
    END TYPE CONNECTIVITIES

end module connectivities_def
