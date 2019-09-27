set(polypy_SOURCES
    polypyAlgebraicNumber.c
    polypyAssignment.c
    polypyInteger.c
    polypyPolynomial.c
    polypy.c
    polypyUPolynomial.c
    utils.c
    polypyVariable.c
    polypyVariableOrder.c
    polypyValue.c
    polypyInterval.c
    polypyFeasibilitySet.c
)

find_package(PythonInterp)
execute_process(COMMAND
  "${PYTHON_EXECUTABLE}" "-c" "from sysconfig import get_paths; print(get_paths()['include'])"
  OUTPUT_VARIABLE PYTHON_INCLUDE_PATH
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
message(STATUS "Using python headers at: ${PYTHON_INCLUDE_PATH}")

include_directories(${GMP_INCLUDE_DIR})
include_directories(${libpoly_SOURCE_DIR}/include)
include_directories(${PYTHON_INCLUDE_PATH})

add_library(polypy MODULE ${polypy_SOURCES})
set_target_properties(polypy PROPERTIES PREFIX "")
if (APPLE)
  set_target_properties(polypy PROPERTIES LINK_FLAGS "-undefined dynamic_lookup")
endif()
target_link_libraries(polypy poly)

string (REPLACE ";" "', '" polypy_SOURCES_STR "${polypy_SOURCES}")
configure_file(setup.py.in ${CMAKE_CURRENT_SOURCE_DIR}/setup.py)
