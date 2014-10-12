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
	${CMAKE_BINARY_DIR}/spm/spm.mod
	${CMAKE_BINARY_DIR}/fem/fem_def.mod)	
include_directories(${CMAKE_BINARY_DIR}/tools 
	${CMAKE_BINARY_DIR}/tracelog
	${CMAKE_BINARY_DIR}/blackboxes
	${CMAKE_BINARY_DIR}/geometries
	${CMAKE_BINARY_DIR}/connectivities
	${CMAKE_BINARY_DIR}/grids
	${CMAKE_BINARY_DIR}/spm
	${CMAKE_BINARY_DIR}/fem)

set(assembly_files assembly_def.F90 
	mass.F90 stiffness.F90 advection.F90 second_deriv.F90 operators.F90 
	l2_projector.F90 l2_projector_grad.F90 l2_projector_curl.F90 projectors.F90 
	identity_field.F90 grad_field.F90 curl_field.F90 second_deriv_field.F90 grad_s_field.F90 second_deriv_s_field.F90 hessian_field.F90 fields_operators.F90
	l2_norm.F90 H1_norm.F90 norms.F90
	assembly_tools.F90 assembly.F90 fields_tools.F90)
set(assembly_mods advection_module.mod l2_norm_module.mod assembly_def.mod l2_projector_curl_module.mod assembly_module.mod l2_projector_grad_module.mod assembly_tools.mod l2_projector_module.mod curl_field_module.mod mass_module.mod fields_operators_module.mod norms_module.mod fields_tools_module.mod operators_module.mod grad_field_module.mod projectors_module.mod grad_s_field_module.mod second_deriv_field_module.mod h1_norm_module.mod second_deriv_module.mod hessian_field_module.mod second_deriv_s_field_module.mod identity_field_module.mod stiffness_module.mod)
#add_library(assembly STATIC ${assembly_files})
add_library(assembly SHARED ${assembly_files})
target_link_libraries(assembly tools tracelog geometries connectivities grids blackboxes spm fem)
install(TARGETS assembly DESTINATION lib)
install(FILES ${assembly_mods} DESTINATION include)
##########################################################################
