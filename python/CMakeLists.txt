set(polypy_SOURCES 
    AlgebraicNumber.c
    Assignment.c
    Integer.c
    Polynomial.c
    polypy.c
    UPolynomial.c
    utils.c
    Variable.c
    VariableOrder.c
)

find_package(PythonInterp)
find_package(PythonLibs)

include_directories(${libpoly_SOURCE_DIR}/include)
include_directories(${PYTHON_INCLUDE_PATH})

add_library(polypy MODULE ${polypy_SOURCES})
set_target_properties(polypy PROPERTIES PREFIX "")
target_link_libraries(polypy poly)
