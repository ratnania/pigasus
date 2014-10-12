##########################################################################
# tracelog directory
##########################################################################
set(HEADER_FILES ${CMAKE_BINARY_DIR}/tools/used_precision.mod)
include_directories(${CMAKE_BINARY_DIR}/tools)

set(tracelog_files tracelog.F90)
set(tracelog_mods tracelog_module.mod)
#add_library(tracelog STATIC ${tracelog_files})
add_library(tracelog SHARED ${tracelog_files})
target_link_libraries(tracelog tools)
install(TARGETS tracelog DESTINATION lib)
install(FILES ${tracelog_mods} DESTINATION include)
##########################################################################
