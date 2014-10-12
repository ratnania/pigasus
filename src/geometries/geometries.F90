!     
! File:   geometries.F90
! Author: ratnani
!
! Created on November 13, 2011, 5:50 PM
!

module geometries_module
    use tracelog_module
    use geometries_def
    use geometry_module
    implicit none

    integer, parameter, private  :: mi_dtllevel_base = 0
contains
    !---------------------------------------------------------------
    subroutine create_geometries(self, ai_ngeo)
        implicit none
        type(GEOMETRIES), intent(inout) :: self
        integer :: ai_ngeo

        CALL printlog("create_geometries : Start", ai_dtllevel = mi_dtllevel_base + 1)

        self % oi_ngeo = ai_ngeo
        ALLOCATE(self % opo_geo(self % oi_ngeo))
!        print*,'self % oi_ngeo=',self % oi_ngeo

        CALL printlog("create_geometries : End", ai_dtllevel = mi_dtllevel_base + 1)
    end subroutine create_geometries
    !---------------------------------------------------------------
    subroutine free_geometries(self)
        implicit none
        type(GEOMETRIES), intent(inout) :: self
        ! LOCAL
        integer :: li_id

        CALL printlog("free_geometries : Start", ai_dtllevel = mi_dtllevel_base + 1)

        do li_id = 1, self % oi_ngeo
            CALL free_geometry(self % opo_geo(li_id))
        end do

        DEALLOCATE(self % opo_geo)

        CALL printlog("free_geometries : End", ai_dtllevel = mi_dtllevel_base + 1)
    end subroutine free_geometries
    !---------------------------------------------------------------
    subroutine add_geo(self, ao_geo, ai_id)
        implicit none
        type(GEOMETRIES), intent(inout) :: self
        type(GEOMETRY), intent(in) :: ao_geo
        integer :: ai_id
        ! LOCAL
        integer :: li_i

        CALL printlog("add_geo : Start", ai_dtllevel = mi_dtllevel_base + 1)

        CALL create_geometry(self % opo_geo(ai_id) &
        , ao_geo % oi_dim &
        , ao_geo % oi_Rd &
        , ao_geo % oi_rational &
        , ao_geo % opi_N &
        , ao_geo % opi_P &
        , ao_geo % opr_u &
        , ao_geo % opr_P &
        , ao_geo % opr_W)

        CALL printlog("add_geo : End", ai_dtllevel = mi_dtllevel_base + 1)

    end subroutine add_geo
    !---------------------------------------------------------------
    subroutine update_geo(self, ao_geo, ai_id)
        implicit none
        type(GEOMETRIES), intent(inout) :: self
        type(GEOMETRY), intent(in) :: ao_geo
        integer :: ai_id
        ! LOCAL
        integer :: li_i

        CALL printlog("update_geo : Start", ai_dtllevel = mi_dtllevel_base + 1)

        CALL free_geometry(self % opo_geo(ai_id))

        CALL create_geometry(self % opo_geo(ai_id) &
        , ao_geo % oi_dim &
        , ao_geo % oi_Rd &
        , ao_geo % oi_rational &
        , ao_geo % opi_N &
        , ao_geo % opi_P &
        , ao_geo % opr_u &
        , ao_geo % opr_P &
        , ao_geo % opr_W)

        CALL printlog("update_geo : End", ai_dtllevel = mi_dtllevel_base + 1)

    end subroutine update_geo
    !---------------------------------------------------------------
    subroutine get_maxnen_geos(self, ai_maxnen)
        implicit none
        type(GEOMETRIES), intent(in) :: self
        INTEGER, intent(inout) :: ai_maxnen
        ! LOCAL
        integer :: li_i
        type(GEOMETRY_INFO) :: lo_geoinfo

        CALL printlog("get_maxnen_geos : Start", ai_dtllevel = mi_dtllevel_base + 3)

        CALL get_info(self % opo_geo(1), lo_geoinfo)
        ai_maxnen = lo_geoinfo % oi_nen

        DO li_i = 2, self % oi_ngeo
            CALL get_info(self % opo_geo(li_i), lo_geoinfo)
            IF (ai_maxnen < lo_geoinfo % oi_nen) THEN
                ai_maxnen = lo_geoinfo % oi_nen
            END IF
        END DO

        CALL printlog("get_maxnen_geos : End", ai_dtllevel = mi_dtllevel_base + 3)

    end subroutine get_maxnen_geos

end module geometries_module
