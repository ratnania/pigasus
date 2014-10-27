!**************************************************
!
!                   SPACES MODULE
!
!**************************************************
module spaces_module
    use fem_def
    use geometries_def
    use geometries_module
    use geometry_tools
    use connectivities_def
    use connectivity_module
    implicit none

    integer, parameter, private :: mi_dtllevel_base = 1

contains
    !---------------------------------------------------------------
    subroutine get_maxnel_space(self, ai_space, ai_maxnel)
        implicit none
        TYPE(FEM) :: self
        INTEGER :: ai_space
        INTEGER :: ai_maxnel
        ! LOCAL
        INTEGER :: li_patch
        INTEGER :: li_grids
        INTEGER :: li_npatchs
        type(CONNECTIVITY), pointer :: lp_con

        li_grids   = self % opi_InfoSpace(ai_space, INFOSPACE_GRIDS)
        li_npatchs = self % opi_InfoGrids(li_grids, INFOGRIDS_NPATCHS)
        lp_con => self % opo_spaces(ai_space) % oo_con

        ai_maxnel = 0
        do li_patch = 0, li_npatchs-1
            if ( ai_maxnel < MAXVAL(lp_con % opi_real_elts (li_patch, :)) ) then
                ai_maxnel = MAXVAL(lp_con % opi_real_elts (li_patch, :))
            end if
        end do
    end subroutine get_maxnel_space
    !---------------------------------------------------------------
    subroutine allocate_spaces(self)
        implicit none
        TYPE(FEM) :: self
        ! LOCAL
        INTEGER :: li_dim
        INTEGER :: li_nel
        INTEGER :: li_npts
        INTEGER :: li_dirnpts
        INTEGER :: li_tensorlevel
        INTEGER :: li_npatchs
        INTEGER :: li_patch
        LOGICAL :: ll_tensor
        INTEGER :: li_id
        INTEGER :: li_grids
        INTEGER :: li_composed
        INTEGER :: li_space
#ifdef _TRACE
        CALL printlog("allocate_spaces : Start", ai_dtllevel = mi_dtllevel_base + 1)
#endif

        DO li_id = 0, self % oi_nspaces - 1
            li_grids = self % opi_InfoSpace(li_id, INFOSPACE_GRIDS)
            li_npatchs = self % opi_InfoGrids(li_grids, INFOGRIDS_NPATCHS)
            CALL create_geometries(self % opo_spaces(li_id) % oo_mapping, li_npatchs)
        END DO

        !> \todo must allocate here the connectivity structure and not in pyfem

        DO li_id = 0, self % oi_nmappings - 1
            li_space   = self % opi_InfoMapping (li_id, INFOMAPPING_SPACE )
!            print *, 'li_space =' , li_space
            li_grids = self % opi_InfoSpace(li_space, INFOSPACE_GRIDS)
            li_npatchs = self % opi_InfoGrids(li_grids, INFOGRIDS_NPATCHS)

!            print *, 'li_npatchs=', li_npatchs
            CALL create_geometries(self % opo_mappings(li_id) % oo_mapping, li_npatchs)
        END DO
#ifdef _TRACE
        CALL printlog("allocate_spaces : End", ai_dtllevel = mi_dtllevel_base + 1)
#endif

    end subroutine allocate_spaces
    !---------------------------------------------------------------
    subroutine add_info_spaces(self)
        implicit none
        TYPE(FEM) :: self
        ! LOCAL
        INTEGER :: li_dim
        INTEGER :: li_nel
        INTEGER :: li_npts
        INTEGER :: li_dirnpts
        INTEGER :: li_npatchs
        LOGICAL :: ll_tensor
        INTEGER :: li_space
        INTEGER :: li_grids
        INTEGER :: li_maxnen
        INTEGER :: li_maxnel
        INTEGER :: li_tensorlevel
        INTEGER :: li_patch
        INTEGER :: li_elt
        INTEGER :: li_b
        INTEGER :: li_composed
        INTEGER :: li_map
        INTEGER :: li_nen
        INTEGER :: li_nrealel
        integer, dimension(:), pointer :: lpi_Pp1
        integer, dimension(:), pointer :: lpi_bi
        integer, dimension(:), pointer :: lpi_ni
        TYPE(GEOMETRY), POINTER :: lo_geo

#ifdef _TRACE
        CALL printlog("add_info_spaces : Start", ai_dtllevel = mi_dtllevel_base + 1)
#endif

        DO li_space = 0, self % oi_nspaces - 1

            !> must be here, in the case of composed spaces with a mapping/pulback
            if ( self % opi_InfoSpace ( li_space , INFOSPACE_EXTMAPPING ) == 1 ) then
                li_grids    = self % opi_InfoSpace(li_space, INFOSPACE_GRIDS)
                li_map      = self % opi_InfoSpace(li_space, INFOSPACE_MAPPING)
                li_npatchs = self % opi_InfoGrids(li_grids, INFOGRIDS_NPATCHS)

                CALL get_maxnel_space(self, li_space, li_maxnel)

                CALL get_maxnen_geos(self % opo_mappings(li_map) % oo_mapping, li_maxnen)
                li_dim      = self % opi_dim (li_grids)

!                print *, "li_maxnel=", li_maxnel
!                print *, "li_maxnen=", li_maxnen
!                print *, "li_grids=", li_grids
!                print *, "li_npatchs=", li_npatchs
!                print *, "li_dim=", li_dim
                li_tensorlevel   = self % opi_InfoMapping(li_map, INFOMAPPING_TENSOR)

                IF (li_tensorlevel==1) THEN
                    ALLOCATE(lpi_Pp1(li_dim))
                    ALLOCATE(lpi_bi(li_dim))
                    ALLOCATE(lpi_ni(li_dim))

                    DO li_patch = 0, li_npatchs - 1
                        li_nel = self % opi_InfoPatch(li_grids, li_patch, INFOPATCH_NEL)
                        lo_geo => self % opo_mappings(li_map) % oo_mapping % opo_geo (li_patch + 1)
                        lpi_Pp1 = lo_geo % opi_P + 1
                        ALLOCATE ( lo_geo % opi_bi (li_maxnel,li_maxnen,li_dim) )
                        ALLOCATE ( lo_geo % opi_ni (li_maxnel,li_dim) )

                        DO li_elt = 1, li_nel
                            CALL get_indices(lo_geo % opi_N, li_elt-1, lpi_ni)
    !                        print *, "lpi_ni=", lpi_ni
                            lo_geo % opi_ni (li_elt,1:li_dim) = lpi_ni(1:li_dim)
                            CALL get_nen_geometry(lo_geo, li_elt, li_nen)
                            DO li_b = 1, li_nen
                                CALL get_indices(lpi_Pp1, li_b-1, lpi_bi)
    !                            print *, "lpi_bi=", lpi_bi
                                lo_geo % opi_bi (li_elt,li_b,1:li_dim) = lpi_bi(1:li_dim)
                            END DO

                        END DO
                    END DO
                    DEALLOCATE(lpi_Pp1)
                    DEALLOCATE(lpi_bi)
                    DEALLOCATE(lpi_ni)
                END IF

            end if

            ! if it is a composed space, we do not need to defin the geometry
            li_composed = self % opi_InfoSpace(li_space, INFOSPACE_COMPOSED)
            IF (li_composed == 1) THEN
                CYCLE
            END IF

            li_grids    = self % opi_InfoSpace(li_space, INFOSPACE_GRIDS)
            li_npatchs  = self % opi_InfoGrids(li_grids, INFOGRIDS_NPATCHS)

            CALL get_maxnel_space(self, li_space, li_maxnel)

            CALL get_maxnen_geos(self % opo_spaces(li_space) % oo_mapping, li_maxnen)
            li_dim      = self % opi_dim (li_grids)
!            li_maxnen   = self % opo_spaces(li_space) % oi_maxnen
            self % opo_spaces(li_space) % oi_maxnen = li_maxnen
!            print *, "li_maxnel=", li_maxnel
!            print *, "li_maxnen=", li_maxnen
!            print *, "li_grids=", li_grids
!            print *, "li_npatchs=", li_npatchs
!            print *, "li_dim=", li_dim
            li_tensorlevel   = self % opi_InfoPatch(li_grids, 0, INFOPATCH_TENSOR)
            self % opo_spaces(li_space) % oi_tensorlevel = li_tensorlevel

            IF (li_tensorlevel==1) THEN
                ALLOCATE(lpi_Pp1(li_dim))
                ALLOCATE(lpi_bi(li_dim))
                ALLOCATE(lpi_ni(li_dim))

                DO li_patch = 0, li_npatchs - 1
                    li_nel = self % opi_InfoPatch(li_grids, li_patch, INFOPATCH_NEL)
                    li_nrealel = 0
                    DO li_elt = 1, li_nel 
                    IF (self % opo_spaces(li_space) % oo_con % opi_real_elts (li_patch, li_elt) > 0) THEN
                            li_nrealel = li_nrealel + 1
                    END IF
                    END DO
!                    li_nrealel = MAXVAL(self % opo_spaces(li_space) % oo_con % opi_real_elts (li_patch, :))

                    lo_geo => self % opo_spaces(li_space) % oo_mapping % opo_geo (li_patch + 1)
                    lpi_Pp1 = lo_geo % opi_P + 1
                    ALLOCATE ( lo_geo % opi_bi (li_maxnel,li_maxnen,li_dim) )
                    ALLOCATE ( lo_geo % opi_ni (li_maxnel,li_dim) )
!                    print *, 'li_nrealel=', li_nrealel
!                    print *, 'SIZE(lo_geo % opi_ni,1) =', SIZE(lo_geo % opi_ni,1)
!                    print *, 'SIZE(lo_geo % opi_ni,2) =', SIZE(lo_geo % opi_ni,2)
                    DO li_elt = 1, li_nrealel
                        CALL get_indices(lo_geo % opi_N, li_elt-1, lpi_ni)
!                        print *, "lpi_ni=", lpi_ni
                        lo_geo % opi_ni (li_elt,1:li_dim) = lpi_ni(1:li_dim)
                        DO li_b = 1, self % opo_spaces (li_space) % oo_con % opi_nen (li_patch + 1)
                            CALL get_indices(lpi_Pp1, li_b-1, lpi_bi)
!                            print *, "lpi_bi=", lpi_bi
                            lo_geo % opi_bi (li_elt,li_b,1:li_dim) = lpi_bi(1:li_dim)
                        END DO

                    END DO
                END DO
                DEALLOCATE(lpi_Pp1)
                DEALLOCATE(lpi_bi)
                DEALLOCATE(lpi_ni)
            END IF

        END DO
#ifdef _TRACE
        CALL printlog("add_info_spaces : End", ai_dtllevel = mi_dtllevel_base + 1)
#endif

    end subroutine add_info_spaces
    !---------------------------------------------------------------
    subroutine deallocate_spaces(self)
        implicit none
        TYPE(FEM) :: self
        ! LOCAL
        INTEGER :: li_dim
        INTEGER :: li_nel
        INTEGER :: li_npts
        INTEGER :: li_dirnpts
        INTEGER :: li_npatchs
        INTEGER :: li_patch
        LOGICAL :: ll_tensor
        INTEGER :: li_space
        INTEGER :: li_grids
        INTEGER :: li_maxnen
        INTEGER :: li_maxnel
        INTEGER :: li_tensorlevel
#ifdef _TRACE
        CALL printlog("deallocate_spaces : Start", ai_dtllevel = mi_dtllevel_base + 1)
#endif

        DO li_space = 0, self % oi_nspaces - 1

            li_grids    = self % opi_InfoSpace(li_space, INFOSPACE_GRIDS)
            li_npatchs  = self % opi_InfoGrids(li_grids, INFOGRIDS_NPATCHS)

            CALL get_maxnel_space(self, li_space, li_maxnel)

            li_maxnen   = self % opo_spaces(li_space) % oi_maxnen
            li_dim      = self % opi_dim (li_grids)
            li_tensorlevel   = self % opi_InfoPatch(li_grids, 0, INFOPATCH_TENSOR)

            IF (li_tensorlevel==1) THEN
                DO li_patch = 0, li_npatchs - 1
                    DEALLOCATE ( self % opo_spaces (li_space) % oo_mapping % opo_geo (li_patch + 1) % opi_bi )
                    DEALLOCATE ( self % opo_spaces (li_space) % oo_mapping % opo_geo (li_patch + 1) % opi_ni )
                END DO
            END IF

        END DO
#ifdef _TRACE
        CALL printlog("deallocate_spaces : End", ai_dtllevel = mi_dtllevel_base + 1)
#endif

    end subroutine deallocate_spaces
    !-----------------------------------------------------------
    subroutine set_space_weights_fem ( self, ai_space, ai_patch, ai_direction, apr_values, ai_dim)
        implicit none
        TYPE(FEM) :: self
        integer, intent(in)  :: ai_space
        integer, intent(in)  :: ai_patch
        integer, intent(in)  :: ai_direction
        integer, intent(in)  :: ai_dim
        real(wp) , dimension(ai_dim) :: apr_values
        ! LOCAL VARIABLES
        integer :: li_space
        integer :: li_grids
        integer :: li_locid
#ifdef _TRACE
        CALL printlog("set_space_weights_fem : Start", ai_dtllevel = 1)
#endif

        CALL set_geometry_weights( self % opo_spaces(ai_space) % oo_mapping % opo_geo (ai_patch+1) &
        , ai_direction, apr_values, ai_dim)
#ifdef _TRACE
        CALL printlog("set_space_weights_fem : End", ai_dtllevel = 1)
#endif

    end subroutine set_space_weights_fem
    !-----------------------------------------------------------
end module spaces_module
!**************************************************
