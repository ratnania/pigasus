##########################################################################
# coonectivities directory
##########################################################################
set(HEADER_FILES 
	${PIGASUS_DIRECTORY_TOOLS}/used_precision.mod 
	${PIGASUS_DIRECTORY_TRACELOG}/tracelog_module.mod)
include_directories(
	${PIGASUS_DIRECTORY_TOOLS}
	${PIGASUS_DIRECTORY_TRACELOG})

set(conn_files 
	connectivities_def.F90 
	connectivity.F90 
	connectivities.F90)
set(conn_mods 
	connectivities_def.mod)

add_library(connectivities SHARED ${conn_files})
target_link_libraries(connectivities tools tracelog)
install(TARGETS connectivities DESTINATION lib)
##########################################################################
