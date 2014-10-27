##########################################################################
# fem directory
##########################################################################
set(HEADER_FILES 
	${PIGASUS_DIRECTORY_TOOLS}/used_precision.mod 
	${PIGASUS_DIRECTORY_TRACELOG}/tracelog_module.mod 
	${PIGASUS_DIRECTORY_BLACKBOXES}/bbox_def.mod 
	${PIGASUS_DIRECTORY_BLACKBOXES}/bbox_modules.mod 
	${PIGASUS_DIRECTORY_GEOMETRIES}/bspline.mod
	${PIGASUS_DIRECTORY_GEOMETRIES}/bsp.mod
	${PIGASUS_DIRECTORY_GEOMETRIES}/crv.mod
	${PIGASUS_DIRECTORY_GEOMETRIES}/vol.mod
	${PIGASUS_DIRECTORY_GEOMETRIES}/iga.mod
	${PIGASUS_DIRECTORY_GEOMETRIES}/srf.mod
	${PIGASUS_DIRECTORY_GEOMETRIES}/geometries_def.mod
	${PIGASUS_DIRECTORY_GEOMETRIES}/geometrie_module.mod
	${PIGASUS_DIRECTORY_GEOMETRIES}/geometries_module.mod
	${PIGASUS_DIRECTORY_GEOMETRIES}/geometry_tools.mod
	${PIGASUS_DIRECTORY_CONNECTIVITIES}/connectivities_def.mod
	${PIGASUS_DIRECTORY_CONNECTIVITIES}/connectivities.mod
	${PIGASUS_DIRECTORY_CONNECTIVITIES}/connectivity_module.mod
	${PIGASUS_DIRECTORY_GRIDS}/grids_def.mod
	${PIGASUS_DIRECTORY_GRIDS}/elements_module.mod
	${PIGASUS_DIRECTORY_GRIDS}/elements_tensor_module.mod
	${PIGASUS_DIRECTORY_GRIDS}/grids.mod
	${PIGASUS_DIRECTORY_SPM}/spm_def.mod
	${PIGASUS_DIRECTORY_SPM}/spm.mod)	
include_directories(
	${PIGASUS_DIRECTORY_TOOLS}
	${PIGASUS_DIRECTORY_TRACELOG}
	${PIGASUS_DIRECTORY_BLACKBOXES}
	${PIGASUS_DIRECTORY_GEOMETRIES}
	${PIGASUS_DIRECTORY_CONNECTIVITIES}
	${PIGASUS_DIRECTORY_GRIDS}
	${PIGASUS_DIRECTORY_SPM})

set(fem_files 
	fem_def.F90 
	metric.F90 
	grids.F90 
	space.F90 
	field.F90 
	graph.F90 
	operator.F90 
	matrix.F90 
	norm.F90 
	solver.F90 
	fem.F90)
set(fem_mods 
	fem_def.mod 
	fem_module.mod 
	field_module.mod 
	graph_module.mod 
	grids_module.mod 
	matrix_module.mod 
	metric_module.mod 
	norm_module.mod 
	operator_module.mod)

add_library(fem SHARED ${fem_files})
target_link_libraries(fem tools tracelog geometries connectivities grids spm)
install(TARGETS fem DESTINATION lib)
##########################################################################
