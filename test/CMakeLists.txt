if (LIBPOLY_BUILD_PYTHON3_API) 
add_custom_target(check
    COMMAND 
        PYTHONPATH=${CMAKE_BINARY_DIR}/python ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check.py
    DEPENDS 
        polypy
    WORKING_DIRECTORY 
        ${CMAKE_CURRENT_SOURCE_DIR}
)
else()
add_custom_target(check
    COMMAND 
        PYTHONPATH=${CMAKE_BINARY_DIR}/python ${Python2_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check.py
    DEPENDS 
        polypy
    WORKING_DIRECTORY 
        ${CMAKE_CURRENT_SOURCE_DIR}
)
endif()

