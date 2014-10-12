##########################################################################
# coonectivities directory
##########################################################################
set(HEADER_FILES ${CMAKE_BINARY_DIR}/tools/used_precision.mod ${CMAKE_BINARY_DIR}/tracelog/tracelog_module.mod)
include_directories(${CMAKE_BINARY_DIR}/tools ${CMAKE_BINARY_DIR}/tracelog)

set(conn_files connectivities_def.F90 connectivity.F90 connectivities.F90)
set(conn_mods connectivities_def.mod)
#add_library(connectivities STATIC ${conn_files})
add_library(connectivities SHARED ${conn_files})
target_link_libraries(connectivities tools tracelog)
install(TARGETS connectivities DESTINATION lib)
install(FILES ${conn_mods} DESTINATION include)
##########################################################################
