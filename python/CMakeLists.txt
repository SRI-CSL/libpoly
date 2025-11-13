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

string (REPLACE ";" "', '" polypy_SOURCES_STR "${polypy_SOURCES}")
configure_file(setup.py.in ${CMAKE_CURRENT_SOURCE_DIR}/setup.py)

Python3_add_library(polypy MODULE ${polypy_SOURCES})

include_directories(${Python3_INCLUDE_DIRS})
include_directories(${GMP_INCLUDE_DIR})
include_directories(${libpoly_SOURCE_DIR}/include)

set_target_properties(polypy PROPERTIES PREFIX "")
if (APPLE)
  set_target_properties(polypy PROPERTIES LINK_FLAGS "-undefined dynamic_lookup")
endif()
target_link_libraries(polypy PRIVATE poly)
