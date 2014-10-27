##########################################################################
# tracelog directory
##########################################################################
set(tools_files 
	precision.F90)
set(tools_mods 
	used_precision.mod)

add_library(tools SHARED ${tools_files})
install(TARGETS tools DESTINATION lib)
##########################################################################
