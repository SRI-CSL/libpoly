add_subdirectory(polyxx)

if(LIBPOLY_BUILD_PYTHON_API)
    add_custom_target(check
        COMMAND
            PYTHONPATH=${CMAKE_BINARY_DIR}/python ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check.py
        DEPENDS
            polypy
        WORKING_DIRECTORY
            ${CMAKE_CURRENT_SOURCE_DIR}
    )
endif()
