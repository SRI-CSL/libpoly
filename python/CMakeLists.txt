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

add_library(polypy MODULE ${polypy_SOURCES})

include_directories(${GMP_INCLUDE_DIR})
include_directories(${libpoly_SOURCE_DIR}/include)
set_target_properties(polypy PROPERTIES PREFIX "")
target_link_libraries(polypy poly)

if (LIBPOLY_BUILD_PYTHON3_API) 
  find_package (Python3 COMPONENTS Interpreter Development)
  include_directories(${Python3_INCLUDE_DIRS})
  target_link_libraries(polypy ${Python3_LIBRARIES})  
else()
  find_package (Python2 COMPONENTS Interpreter Development)
  include_directories(${Python2_INCLUDE_DIRS})
  target_link_libraries(polypy ${Python2_LIBRARIES})
endif()

string (REPLACE ";" "', '" polypy_SOURCES_STR "${polypy_SOURCES}")
configure_file(setup.py.in ${CMAKE_CURRENT_SOURCE_DIR}/setup.py)
