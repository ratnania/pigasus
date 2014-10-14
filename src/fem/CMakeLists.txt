##########################################################################
# fem directory
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
	${CMAKE_BINARY_DIR}/geometries/geometry_tools.mod
	${CMAKE_BINARY_DIR}/connectivities/connectivities_def.mod
	${CMAKE_BINARY_DIR}/connectivities/connectivities.mod
	${CMAKE_BINARY_DIR}/connectivities/connectivity_module.mod
	${CMAKE_BINARY_DIR}/grids/grids_def.mod
	${CMAKE_BINARY_DIR}/grids/elements_module.mod
	${CMAKE_BINARY_DIR}/grids/elements_tensor_module.mod
	${CMAKE_BINARY_DIR}/grids/grids.mod
	${CMAKE_BINARY_DIR}/spm/spm_def.mod
	${CMAKE_BINARY_DIR}/spm/spm.mod)
include_directories(${CMAKE_BINARY_DIR}/tools 
	${CMAKE_BINARY_DIR}/tracelog
	${CMAKE_BINARY_DIR}/blackboxes
	${CMAKE_BINARY_DIR}/geometries
	${CMAKE_BINARY_DIR}/connectivities
	${CMAKE_BINARY_DIR}/grids
	${CMAKE_BINARY_DIR}/spm)

set(fem_files fem_def.F90 metric.F90 grids.F90 space.F90 field.F90 graph.F90 operator.F90 matrix.F90 norm.F90 solver.F90 fem.F90)
set(fem_mods fem_def.mod fem_module.mod field_module.mod graph_module.mod grids_module.mod matrix_module.mod metric_module.mod norm_module.mod operator_module.mod)
add_library(fem SHARED ${fem_files})
target_link_libraries(fem tools tracelog geometries connectivities grids spm)
install(TARGETS fem DESTINATION lib)
install(FILES ${fem_mods} DESTINATION include)
##########################################################################
