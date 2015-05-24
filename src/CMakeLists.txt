set(poly_SOURCES 
    utils/debug_trace.c
    utils/assignment.c
    utils/statistics.c
    utils/output.c
    number/integer.c
    number/rational.c
    number/dyadic_rational.c
    number/algebraic_number.c
    number/value.c
    interval/interval.c
    interval/arithmetic.c
    variable/variable_db.c
    variable/variable_list.c
    variable/variable_order.c
    upolynomial/umonomial.c
    upolynomial/upolynomial.c
    upolynomial/output.c
    upolynomial/upolynomial_dense.c
    upolynomial/bounds.c
    upolynomial/gcd.c
    upolynomial/factors.c
    upolynomial/factorization.c
    upolynomial/root_finding.c
    polynomial/monomial.c
    polynomial/coefficient.c
    polynomial/output.c
    polynomial/gcd.c
    polynomial/psc.c
    polynomial/factorization.c
    polynomial/polynomial.c
    polynomial/polynomial_context.c
    polynomial/feasibility_set.c
    poly.c
)

# Add warnings as errors for good practice
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-unused-function -Wall -Werror -Wextra -std=gnu99")

# Find the GMP number library
find_path(GMP_INCLUDE_DIR gmp.h)
find_library(GMP_LIBRARY gmp)
if (GMP_INCLUDE_DIR AND GMP_LIBRARY)
    message(STATUS "GMP headers: " ${GMP_INCLUDE_DIR})
    message(STATUS "GMP library: " ${GMP_LIBRARY})
else()
    MESSAGE(FATAL_ERROR "Could not the GMP number library (sudo apt-get install libgmp-dev)")
endif()

include_directories(${GMP_INCLUDE_DIR})
include_directories(${libpoly_SOURCE_DIR}/include)
include_directories(${libpoly_SOURCE_DIR}/src)

add_library(poly SHARED ${poly_SOURCES})
target_link_libraries(poly ${GMP_LIBRARY})

install(TARGETS poly LIBRARY DESTINATION lib)
