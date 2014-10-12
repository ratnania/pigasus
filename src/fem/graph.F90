module qsort_c_module
use used_precision
implicit none
public :: QsortC
private :: Partition

contains

recursive subroutine QsortC(A)
  real(wp), intent(in out), dimension(:) :: A
  integer :: iq

  if(size(A) > 1) then
     call Partition(A, iq)
     call QsortC(A(:iq-1))
     call QsortC(A(iq:))
  endif
end subroutine QsortC

subroutine Partition(A, marker)
  real(wp), intent(in out), dimension(:) :: A
  integer, intent(out) :: marker
  integer :: i, j
  real(wp) :: temp
  real(wp) :: x      ! pivot point
  x = A(1)
  i= 0
  j= size(A) + 1

  do
     j = j-1
     do
        if (A(j) <= x) exit
        j = j-1
     end do
     i = i+1
     do
        if (A(i) >= x) exit
        i = i+1
     end do
     if (i < j) then
        ! exchange A(i) and A(j)
        temp = A(i)
        A(i) = A(j)
        A(j) = temp
     elseif (i == j) then
        marker = i+1
        return
     else
        marker = i
        return
     endif
  end do

end subroutine Partition

end module qsort_c_module





!**************************************************
!
!                   GRAPH MODULE
!
!**************************************************
module graph_module
    use fem_def
    use qsort_c_module
    implicit none

    integer, parameter, private :: mi_dtllevel_base = 1

contains

    !---------------------------------------------------------------
    subroutine SET_INFOGRAPH(self, ai_id, ai_param, ai_val)
        implicit none
        TYPE(FEM) :: self
        INTEGER :: ai_id
        INTEGER :: ai_param
        INTEGER :: ai_val

        self % opi_InfoGraph(ai_id, ai_param) = ai_val
    end subroutine SET_INFOGRAPH
    !----------------------------------------------------------------------------------------------
    !> \todo not finished yet
    subroutine allocate_simple_graph(self, ai_id)
        implicit none
        TYPE(FEM) :: self
        !> GRAPH REFERENCE IN THE DICTIONNARY
        integer :: ai_id
        ! LOCAL
        integer :: li_sp
        integer :: li_spprime
        integer :: li_size
        integer :: li_sizeprime
        integer :: li_grids
        integer :: li_gridsprime
        integer :: li_npatchs
        integer :: li_composed
        integer :: li_composed_sp1
        integer :: li_patch
        type(CONNECTIVITY), pointer :: lp_con
        type(CONNECTIVITY), pointer :: lp_conprime
        integer, dimension(:), pointer :: lpi_nel
        integer, dimension(:), pointer :: lpi_nrealel
        INTEGER, PARAMETER :: li_root = -1
        INTEGER :: ierr

#ifdef _TRACE
        CALL printlog("allocate_simple_graph : Start", ai_dtllevel = mi_dtllevel_base + 1)
#endif

        li_sp           = self % opi_InfoGraph (ai_id, INFOGRAPH_SPACE_1)
        li_spprime      = self % opi_InfoGraph (ai_id, INFOGRAPH_SPACE_2)

        li_size         = self % opi_InfoSpace(li_sp, INFOSPACE_SIZE)
        li_sizeprime    = self % opi_InfoSpace(li_spprime, INFOSPACE_SIZE)

        lp_con          => self % opo_spaces(li_sp) % oo_con
        lp_conprime     => self % opo_spaces(li_spprime) % oo_con

        li_grids    = self % opi_InfoSpace(li_sp, INFOSPACE_GRIDS)
        li_npatchs  = self % opi_infoGrids(li_grids, INFOGRIDS_NPATCHS)

        ALLOCATE(lpi_nel(0:li_npatchs - 1))
        ALLOCATE(lpi_nrealel(0:li_npatchs - 1))

        lpi_nel     = self % opi_InfoPatch(li_grids, 0:li_npatchs - 1, INFOPATCH_NEL)

        do li_patch = 0, li_npatchs-1
            lpi_nrealel (li_patch) = MAXVAL(lp_con % opi_real_elts (li_patch, :))
        end do


        call create_CSRGraph(self % opo_Graph(ai_id) &
        , li_size, li_sizeprime &
        , li_npatchs, lpi_nrealel &
        , lp_con, lp_conprime)

        DEALLOCATE(lpi_nel)
        DEALLOCATE(lpi_nrealel)

#ifdef _TRACE
        CALL printlog("allocate_simple_graph : End", ai_dtllevel = mi_dtllevel_base + 1)
#endif

    end subroutine allocate_simple_graph
    !----------------------------------------------------------------------------------------------
    !> \todo not finished yet
    subroutine allocate_composed_graph(self, ai_id)
        implicit none
        TYPE(FEM) :: self
        !> GRAPH REFERENCE IN THE DICTIONNARY
        integer :: ai_id
        ! LOCAL
        integer :: li_dim1
        integer :: li_dim2
        integer, dimension(:,:), pointer :: lpi_id
#ifdef _TRACE
        CALL printlog("allocate_composed_graph : Start", ai_dtllevel = mi_dtllevel_base + 1)
#endif

        print *, 'allocate_composed_graph: Not yet implemented'

#ifdef _TRACE
        CALL printlog("allocate_composed_graph : End", ai_dtllevel = mi_dtllevel_base + 1)
#endif

    end subroutine allocate_composed_graph
    !---------------------------------------------------------------
    subroutine allocate_graphs(self)
        implicit none
        TYPE(FEM) :: self
        ! LOCAL
        INTEGER :: li_id
        INTEGER :: ierr 

#ifdef _TRACE
        CALL printlog("allocate_graphs : Start", ai_dtllevel = mi_dtllevel_base + 1)
#endif
        ALLOCATE ( self % opo_Graph(0:self % oi_ngraphs - 1))

        DO li_id = 0, self % oi_nGraphs - 1
            CALL allocate_simple_graph(self, li_id)
        END DO


#ifdef _TRACE
        CALL printlog("allocate_graphs : End", ai_dtllevel = mi_dtllevel_base + 1)
#endif

    end subroutine allocate_graphs
    !---------------------------------------------------------------
    subroutine deallocate_graphs(self)
        implicit none
        TYPE(FEM) :: self
        ! LOCAL
        INTEGER :: li_id
        INTEGER :: ierr

#ifdef _TRACE
        CALL printlog("deallocate_graphs : Start", ai_dtllevel = mi_dtllevel_base + 1)
#endif

        DO li_id = 0, self % oi_nGraphs - 1
            DEALLOCATE(self % opo_Graph(li_id) % opi_ja) 
            DEALLOCATE(self % opo_Graph(li_id) % opi_ia) 
        END DO

        DEALLOCATE(self % opo_Graph)

#ifdef _TRACE
        CALL printlog("deallocate_graphs : End", ai_dtllevel = mi_dtllevel_base + 1)
#endif

    end subroutine deallocate_graphs
    !---------------------------------------------------------------------------------
    subroutine create_CSRGraph(self, ai_nR, ai_nC &
        , ai_npatch, api_nel    &
        , ao_con_C, ao_con_R    &
        , ai_COEF)
        implicit none
        !> param[inout] self : CSR MATRIX STRUCTURE
        type(GRAPH) :: self
        !> param[in] ai_nC : NUMBER OF COLUMNS
        integer :: ai_nC
        !> param[in] ai_nR : NUMBER OF ROWS
        integer :: ai_nR
        !> param[in] ai_npatch : NUMBER OF PATCHS
        integer :: ai_npatch
        !> param[in] api_nel : NUMBER OF NON ZERO ELEMENTS IN EACH PATCH
        integer, dimension(:) :: api_nel
        !>
        TYPE(CONNECTIVITY) :: ao_con_C
        TYPE(CONNECTIVITY) :: ao_con_R
        !>
        integer, optional :: ai_COEF
        !local var
        integer :: li_err, li_flag
        integer*8 :: li_nnz
        integer, dimension(:,:), pointer :: lpi_columns
        integer, dimension(:), pointer :: lpi_occ
        integer :: li_COEF
        integer :: li_maxnen_C
        integer*8 :: li_virtual_nnz

#ifdef _TRACE
        CALL printlog("create_CSRGraph: Start", ai_dtllevel = mi_dtllevel_base + 1)
#endif
        !> \todo il y a un bug ici, il faudra reduire li_COEF et voir
        li_COEF = 20
        if (present(ai_COEF)) then
            li_COEF = ai_COEF
        end if

        li_maxnen_C = MAXVAL(ao_con_C % opi_nen(:))

        allocate(lpi_columns(ai_nR, 0:li_COEF * li_maxnen_C))
        allocate(lpi_occ(ai_nR + 1))

        lpi_columns(:,:) = 0
        lpi_occ(:) = 0
    ! ...
    ! COUNTING NON ZERO ELEMENTS
    ! ...
        CALL count_non_zero_elts(ai_nR, ai_nC   &
        , ai_npatch, api_nel    &
        , ao_con_C, ao_con_R    &
        , lpi_columns, lpi_occ, li_virtual_nnz, li_nnz)
    ! ...
        self % oi_nR = ai_nR
        self % oi_nC = ai_nC
        self % oi_nnz = li_nnz

        allocate(self % opi_ia(self % oi_nR + 1))
        allocate(self % opi_ja(self % oi_nnz))

        call init_Graph(self, ai_nR, ai_nC   &
        , ai_npatch, api_nel    &
        , ao_con_C, ao_con_R    &
        , lpi_columns, lpi_occ  &
        , self % opi_ia, self % opi_ja)

        deallocate(lpi_columns)
        deallocate(lpi_occ)

#ifdef _TRACE
        CALL printlog("create_CSRGraph: End", ai_dtllevel = mi_dtllevel_base + 1)
#endif
    end subroutine create_CSRGraph
    !---------------------------------------------------------------------------------
    subroutine init_Graph(self, ai_nR, ai_nC   &
        , ai_npatch, api_nel    &
        , ao_con_C, ao_con_R    &
        , api_columns, api_occ  &
        , api_ia, api_ja)
        ! _C FOR ROWS
        ! _R FOR COLUMNS
        implicit none
        type(GRAPH) :: self
        integer :: ai_nC
        integer :: ai_nR
        integer :: ai_npatch
        integer, dimension(:) :: api_nel
        TYPE(CONNECTIVITY) :: ao_con_C
        TYPE(CONNECTIVITY) :: ao_con_R
        integer, dimension(:,:), pointer :: api_columns
        integer, dimension(:), pointer :: api_occ
        integer, dimension(:), pointer :: api_ia
        integer*8, dimension(:), pointer :: api_ja
        !local var
        integer :: li_nel
        integer :: li_id
        integer :: li_e
        integer :: li_b_C
        integer :: li_A_C
        integer :: li_b_R
        integer :: li_A_R
        integer :: li_index
        integer :: li_i
        integer :: li_size
        integer :: li_nen_C
        integer :: li_err
        integer :: li_flag        
        real(wp), dimension(:), pointer :: lpr_tmp
#ifdef _TRACE
        CALL printlog("init_Graph: Start", ai_dtllevel = mi_dtllevel_base + 1)
#endif
        ! INITIALIZING ia
        api_ia(1) = 1

        do li_i = 1, self % oi_nR

            api_ia(li_i + 1) = api_ia(1) + SUM(api_occ(1: li_i))

        end do

        ! INITIALIZING ja
        DO li_id = 1, ai_npatch

            li_nel = api_nel (li_id)

            ! WE FIRST COMPUTE, FOR EACH ROW, THE NUMBER OF COLUMNS THAT WILL BE USED
            do li_e = 1, li_nel

                li_nen_C = ao_con_C % opi_nen(li_id)
                do li_b_C = 1, li_nen_C

                    li_A_C = ao_con_C % opi_LM(li_id, li_b_C, li_e)

                    if (li_A_C == 0) then
                        cycle
                    end if

                    if (api_columns(li_A_C, 0) == 0) then
                        cycle
                    end if

                    li_size = api_columns(li_A_C, 0)

                    allocate ( lpr_tmp(li_size), stat = li_err)
                    if (li_err .ne. 0) li_flag = 10

                    lpr_tmp(1: li_size) = real( api_columns(li_A_C, 1: li_size))

                    call QsortC(lpr_tmp)

                    do li_i = 1, li_size

                        api_ja(api_ia(li_A_C) + li_i - 1) = int ( lpr_tmp(li_i))

                    end do

                    api_columns(li_A_C, 0) = 0
                    deallocate ( lpr_tmp)

                end do

            end do

        end do
#ifdef _TRACE
        CALL printlog("init_Graph: END", ai_dtllevel = mi_dtllevel_base + 1)
#endif

    end subroutine init_Graph    
    !---------------------------------------------------------------------------------
    subroutine count_non_zero_elts(ai_nR, ai_nC   &
        , ai_npatch, api_nel    &
        , ao_con_C, ao_con_R    &
        , api_columns, api_occ, ai_virtual_nnz, ai_nnz)
        ! _C FOR ROWS
        ! _R FOR COLUMNS
        implicit none
        integer :: ai_nC
        integer :: ai_nR
        integer :: ai_npatch
        integer, dimension(:) :: api_nel
        TYPE(CONNECTIVITY) :: ao_con_C
        TYPE(CONNECTIVITY) :: ao_con_R
        integer, dimension(:,:), pointer :: api_columns
        integer, dimension(:), pointer :: api_occ
        integer*8 :: ai_virtual_nnz
        integer*8 :: ai_nnz
        !local var
        integer :: li_e
        integer :: li_nel
        integer :: li_id
        integer :: li_nen_C
        integer :: li_nen_R
        integer :: li_b_C
        integer :: li_A_C
        integer :: li_b_R
        integer :: li_A_R
        integer :: li_i
        integer :: li_err
        integer :: li_flag
        integer*8 :: li_result
        logical :: ll_done
        integer, dimension(2) :: lpi_size
        real(wp), dimension(:), pointer :: lpr_tmp
        integer, dimension(:,:), pointer :: lpi_columns
#ifdef _TRACE
        CALL printlog("count_non_zero_elts: BEGIN", ai_dtllevel = mi_dtllevel_base + 1)
#endif        

        ai_virtual_nnz = 0

        DO li_id = 1, ai_npatch

            li_nel = api_nel (li_id)
            
            ! WE FIRST COMPUTE, FOR EACH ROW, THE NUMBER OF COLUMNS THAT WILL BE USED
            do li_e = 1, li_nel
!                print *,"li_id, li_e=",li_id, li_e
                li_nen_C = ao_con_C % opi_nen(li_id)
                do li_b_C = 1, li_nen_C

                    li_A_C = ao_con_C % opi_LM(li_id, li_b_C, li_e)
                    if (li_A_C == 0) then
                        cycle
                    end if

                    li_nen_R = ao_con_R % opi_nen(li_id)
                    do li_b_R = 1, li_nen_R

                        li_A_R = ao_con_R % opi_LM(li_id, li_b_R, li_e)
                        if (li_A_R == 0) then
                            cycle
                        end if

                        ai_virtual_nnz = ai_virtual_nnz + 1

                        ll_done = .false.
                        ! WE CHECK IF IT IS THE FIRST OCCURANCE OF THE COUPLE (li_A_C, li_A_R)
                        do li_i = 1, api_columns(li_A_C, 0)

                            if (api_columns(li_A_C, li_i) /= li_A_R) then
                                cycle
                            end if

                            ll_done = .true.
                            exit

                        end do

                        if (.not.ll_done) then

                            api_occ(li_A_C) = api_occ(li_A_C) + 1

                            ! li_A_C IS THE ROW NUM, li_A_R THE COLUMN NUM
                            ! INITIALIZATION OF THE SPARSE MATRIX
                            api_columns(li_A_C, 0) = api_columns(li_A_C, 0) + 1
                            api_columns(li_A_C, api_columns(li_A_C, 0)) = li_A_R

                            ! resizing the array
                            lpi_size(1) = SIZE(api_columns, 1)
                            lpi_size(2) = SIZE(api_columns, 2)
                            if (lpi_size(2) < api_columns(li_A_C, 0)) then
                                ALLOCATE(lpi_columns(lpi_size(1), lpi_size(2)))
                                lpi_columns = api_columns

                                DEALLOCATE(api_columns)

                                ALLOCATE(api_columns(lpi_size(1), 2 * lpi_size(2)))
                                api_columns(1:lpi_size(1), 1:lpi_size(2)) = lpi_columns(1:lpi_size(1), 1:lpi_size(2))

                                DEALLOCATE(lpi_columns)
                            end if


                        end if

                    end do

                end do

            end do

        END DO
        
        ! COUNT NON ZERO ELEMENTS
        li_result = SUM(api_occ(1: ai_nR))

        ai_nnz = li_result

#ifdef _TRACE
        CALL printlog("count_non_zero_elts: END", ai_dtllevel = mi_dtllevel_base + 1)
#endif        
    end subroutine count_non_zero_elts
    
end module graph_module
!**************************************************
