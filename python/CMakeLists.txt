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
find_package(PythonLibs)

include_directories(${GMP_INCLUDE_DIR})
include_directories(${libpoly_SOURCE_DIR}/include)
include_directories(${PYTHON_INCLUDE_PATH})

add_library(polypy MODULE ${polypy_SOURCES})
set_target_properties(polypy PROPERTIES PREFIX "")
target_link_libraries(polypy poly)
