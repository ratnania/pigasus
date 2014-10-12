##########################################################################
# geometries directory
##########################################################################
set(HEADER_FILES ${CMAKE_BINARY_DIR}/tools/used_precision.mod ${CMAKE_BINARY_DIR}/tracelog/tracelog_module.mod)
include_directories(${CMAKE_BINARY_DIR}/tools ${CMAKE_BINARY_DIR}/tracelog)

set(geo_files igalib.F90 geometries_def.F90 geometry_tools.F90 geometry.F90  geometries.F90)
set(geo_mods bspline.mod  bsp.mod  crv.mod  iga.mod  srf.mod  vol.mod)
#add_library(geometries STATIC ${geo_files})
add_library(geometries SHARED ${geo_files})
target_link_libraries(geometries tools tracelog)
install(TARGETS geometries DESTINATION lib)
install(FILES ${geo_mods} DESTINATION include)
##########################################################################
