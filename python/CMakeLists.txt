##########################################################################
# python package directory
##########################################################################
find_program(PYTHON "python")

if (PYTHON)
    set(SETUP_PY_IN "${CMAKE_CURRENT_SOURCE_DIR}/setup.py.in")
    set(SETUP_PY    "${CMAKE_CURRENT_BINARY_DIR}/setup.py")
    set(DEPS        "${CMAKE_CURRENT_SOURCE_DIR}")
    set(OUTPUT      "${CMAKE_CURRENT_BINARY_DIR}/timestamp")

    set(SETUP_README_IN "${CMAKE_SOURCE_DIR}/README.md")
    set(SETUP_README    "${CMAKE_CURRENT_BINARY_DIR}/README.md")

    file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/fem DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
    file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/fit DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
    file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/gallery DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
    file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/interpolate DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
    file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/multigrid DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
    file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/pigasus DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
    file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/plugin DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
    file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/utils DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
    file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/core DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

    configure_file(${SETUP_PY_IN} ${SETUP_PY})
    configure_file(${SETUP_README_IN} ${SETUP_README})

    add_custom_command(OUTPUT ${OUTPUT}
                       COMMAND ${PYTHON} ${SETUP_PY} build
                       COMMAND ${CMAKE_COMMAND} -E touch ${OUTPUT}
                       DEPENDS ${DEPS})

    add_custom_target(target ALL DEPENDS ${OUTPUT})

    install(CODE "execute_process(COMMAND ${PYTHON} ${SETUP_PY} install -f --prefix=${CMAKE_INSTALL_PREFIX} WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})")

#    install( DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/fem
#       DESTINATION ${CMAKE_CURRENT_BINARY_DIR}
#       COMPONENT myinstall
#       PATTERN .svn EXCLUDE
#    )
#    
#    add_custom_target(myinstall
#        COMMAND ${CMAKE_COMMAND}
#        -D COMPONENT=myinstall
#        -P cmake_install.cmake
#    )

endif()
##########################################################################
