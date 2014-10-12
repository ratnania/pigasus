##########################################################################
# geometries directory
##########################################################################
set(HEADER_FILES ${CMAKE_BINARY_DIR}/tools/used_precision.mod ${CMAKE_BINARY_DIR}/tracelog/tracelog_module.mod)
include_directories(${CMAKE_BINARY_DIR}/tools ${CMAKE_BINARY_DIR}/tracelog)

set(bbox_files bbox_def.F90 bbox.F90)
set(bbox_mods bbox_def.mod)
#add_library(bbox SHARED ${bbox_files})
#target_link_libraries(bbox tools tracelog)
#install(TARGETS bbox DESTINATION lib)
add_library(blackboxes SHARED ${bbox_files})
target_link_libraries(blackboxes tools tracelog)
install(TARGETS blackboxes DESTINATION lib)
install(FILES ${bbox_mods} DESTINATION include)
##########################################################################
