##########################################################################
# coonectivities directory
##########################################################################
set(HEADER_FILES ${CMAKE_BINARY_DIR}/tools/used_precision.mod 
	${CMAKE_BINARY_DIR}/tracelog/tracelog_module.mod 
	${CMAKE_BINARY_DIR}/blackboxes/bbox_def.mod 
	${CMAKE_BINARY_DIR}/blackboxes/bbox_modules.mod 
	${CMAKE_BINARY_DIR}/geometries/bspline.mod
	${CMAKE_BINARY_DIR}/geometries/bsp.mod
	${CMAKE_BINARY_DIR}/geometries/crv.mod
	${CMAKE_BINARY_DIR}/geometries/vol.mod
	${CMAKE_BINARY_DIR}/geometries/iga.mod
	${CMAKE_BINARY_DIR}/geometries/srf.mod
	${CMAKE_BINARY_DIR}/geometries/geometries_def.mod
	${CMAKE_BINARY_DIR}/geometries/geometrie_module.mod
	${CMAKE_BINARY_DIR}/geometries/geometries_module.mod
	${CMAKE_BINARY_DIR}/geometries/geometry_tools.mod)
include_directories(${CMAKE_BINARY_DIR}/tools 
	${CMAKE_BINARY_DIR}/tracelog
	${CMAKE_BINARY_DIR}/blackboxes
	${CMAKE_BINARY_DIR}/geometries)

set(grids_files grids_def.F90 elements.F90 elements_tensor.F90 grids.F90)
set(grids_mods grids_def.mod)
#add_library(grids STATIC ${grids_files})
add_library(grids SHARED ${grids_files})
target_link_libraries(grids tools tracelog geometries)
install(TARGETS grids DESTINATION lib)
install(FILES ${grids_mods} DESTINATION include)
##########################################################################
