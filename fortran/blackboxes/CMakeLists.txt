##########################################################################
# BlackBox directory
##########################################################################
set(HEADER_FILES 
	${PIGASUS_DIRECTORY_TOOLS}/used_precision.mod 
	${PIGASUS_DIRECTORY_TRACELOG}/tracelog_module.mod)
include_directories(
	${PIGASUS_DIRECTORY_TOOLS} 
	${PIGASUS_DIRECTORY_TRACELOG})

set(bbox_files 
	bbox_def.F90 
	bbox.F90)
set(bbox_mods 
	bbox_def.mod)

add_library(blackboxes SHARED ${bbox_files})
target_link_libraries(blackboxes tools tracelog)
install(TARGETS blackboxes DESTINATION lib)
##########################################################################
