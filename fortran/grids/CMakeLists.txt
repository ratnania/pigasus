##########################################################################
# grids directory
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
	${PIGASUS_DIRECTORY_GEOMETRIES}/geometry_tools.mod)
include_directories(
	${PIGASUS_DIRECTORY_TOOLS}
	${PIGASUS_DIRECTORY_TRACELOG}
	${PIGASUS_DIRECTORY_BLACKBOXES}
	${PIGASUS_DIRECTORY_GEOMETRIES})

set(grids_files 
	grids_def.F90 
	elements.F90 
	elements_tensor.F90 
	grids.F90)
set(grids_mods 
	grids_def.mod)

add_library(grids SHARED ${grids_files})
target_link_libraries(grids tools tracelog geometries)
install(TARGETS grids DESTINATION lib)
##########################################################################
