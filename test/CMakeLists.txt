add_custom_target(check
    COMMAND 
        PYTHONPATH=${CMAKE_BINARY_DIR}/python ${CMAKE_CURRENT_SOURCE_DIR}/check.py
    DEPENDS 
        polypy
    WORKING_DIRECTORY 
        ${CMAKE_CURRENT_SOURCE_DIR}
)
   