set(poly_SOURCES 
    polynomial/polynomial.c
    polynomial/polynomial_context.c
    polynomial/coefficient.c
    polynomial/factorization.c
    polynomial/monomial.c
    upolynomial/bounds.c
    upolynomial/factorization.c
    upolynomial/upolynomial.c
    upolynomial/gcd.c
    upolynomial/factors.c
    upolynomial/upolynomial_dense.c
    upolynomial/root_finding.c
    upolynomial/upolynomial.c
    upolynomial/umonomial.c
    number/rational.c
    number/algebraic_number.c
    number/dyadic_rational.c
    number/integer.c
    variable/variable.c
    variable/variable_order.c
    interval/interval.c
    interval/arithmetic.c
    utils/debug_trace.c
    utils/assignment.c
    utils/statistics.c
)

# Add warnings as errors for good practice
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Werror -Wextra -std=gnu99")

# Find the GMP number library
find_path(GMP_INCLUDE_DIR gmp.h)
find_library(GMP_LIBRARY gmp)
if (GMP_INCLUDE_DIR AND GMP_LIBRARY)
    message(STATUS "GMP headers: "${GMP_INCLUDE_DIR})
    message(STATUS "GMP library: "${GMP_LIBRARY})
else()
    MESSAGE(FATAL_ERROR "Could not the GMP number library (sudo apt-get install libgmp-dev)")
endif()

include_directories(${GMP_INCLUDE_DIR})
include_directories(${libpoly_SOURCE_DIR}/include)
include_directories(${libpoly_SOURCE_DIR}/src)

add_library(poly SHARED ${poly_SOURCES})
target_link_libraries(poly ${GMP_LIBRARY})
