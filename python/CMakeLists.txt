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

include_directories(${GMP_INCLUDE_DIR})
include_directories(${libpoly_SOURCE_DIR}/include)

add_library(polypy MODULE ${polypy_SOURCES})
target_include_directories(polypy PRIVATE ${Python3_INCLUDE_DIRS})
target_link_libraries(polypy PRIVATE poly ${Python3_LIBRARIES})

set_target_properties(polypy PROPERTIES PREFIX "")
#set_target_properties(polypy PROPERTIES SUFFIX "${Python3_EXTENSION_SUFFIX}")

if (APPLE)
  set_target_properties(polypy PROPERTIES LINK_FLAGS "-undefined dynamic_lookup")
endif()

string (REPLACE ";" "', '" polypy_SOURCES_STR "${polypy_SOURCES}")
configure_file(setup.py.in ${CMAKE_CURRENT_SOURCE_DIR}/setup.py)
