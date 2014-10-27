##########################################################################
# tracelog directory
##########################################################################
set(HEADER_FILES 
	${PIGASUS_DIRECTORY_TOOLS}/used_precision.mod)
include_directories(
	${PIGASUS_DIRECTORY_TOOLS})

set(tracelog_files 
	tracelog.F90)
set(tracelog_mods 
	tracelog_module.mod)

add_library(tracelog SHARED ${tracelog_files})
target_link_libraries(tracelog tools)
install(TARGETS tracelog DESTINATION lib)
##########################################################################
