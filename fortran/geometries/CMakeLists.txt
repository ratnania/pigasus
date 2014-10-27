##########################################################################
# geometries directory
##########################################################################
set(HEADER_FILES 
	${PIGASUS_DIRECTORY_TOOLS}/used_precision.mod 
	${PIGASUS_DIRECTORY_TRACELOG}/tracelog_module.mod)
include_directories(
	${PIGASUS_DIRECTORY_TOOLS} 
	${PIGASUS_DIRECTORY_TRACELOG})

set(geo_files 
	igalib.F90 
	geometries_def.F90 
	geometry_tools.F90 
	geometry.F90  
	geometries.F90)
set(geo_mods 
	bspline.mod  
	bsp.mod  
	crv.mod  
	iga.mod  
	srf.mod  
	vol.mod)
#
add_library(geometries SHARED ${geo_files})
target_link_libraries(geometries tools tracelog)
install(TARGETS geometries DESTINATION lib)
##########################################################################
