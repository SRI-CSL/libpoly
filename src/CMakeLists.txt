set(poly_SOURCES 
  utils/debug_trace.c
  utils/assignment.c
  utils/statistics.c
  utils/output.c
  utils/sign_condition.c
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
  polynomial/polynomial_hash_set.c
  polynomial/polynomial_vector.c
  poly.c
)

if (NOT HAVE_OPEN_MEMSTREAM) 
  set(poly_SOURCES "utils/open_memstream.c ${poly_SOURCES}")
endif()

set(polyxx_SOURCES
  polyxx/context.cpp
  polyxx/variable.cpp
)

# Add warnings as errors for good practice
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Werror -Wextra -std=gnu99")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Werror -Wextra -std=c++11")

include_directories(${GMP_INCLUDE_DIR})
include_directories(${libpoly_SOURCE_DIR}/include)
include_directories(${libpoly_SOURCE_DIR}/src)

add_library(poly SHARED ${poly_SOURCES})

set_target_properties(poly PROPERTIES
  VERSION ${LIBPOLY_VERSION}
  SOVERSION ${LIBPOLY_VERSION_MAJOR}
)

target_link_libraries(poly ${GMP_LIBRARY})

add_library(polyxx SHARED ${polyxx_SOURCES})
set_target_properties(polyxx PROPERTIES 
  VERSION ${LIBPOLY_VERSION} 
  SOVERSION ${LIBPOLY_VERSION_MAJOR}
) 

target_link_libraries(polyxx poly ${GMP_LIBRARY})

#
# For windows and friends: the DLL goes into <install>/bin
# and the lib goes into <install>/lib
# cmake complains if the RUNTIME DESTINATION bin is missing
#
if (CYGWIN OR WIN32 OR MINGW)
 install(TARGETS poly LIBRARY DESTINATION lib RUNTIME DESTINATION bin)
 install(TARGETS polyxx LIBRARY DESTINATION lib RUNTIME DESTINATION bin)
else()
 install(TARGETS poly LIBRARY DESTINATION lib)
 install(TARGETS polyxx LIBRARY DESTINATION lib)
endif()

#
# BD: also build the archive libpoly.a
#
if (LIBPOLY_BUILD_STATIC)
  add_library(static_poly STATIC ${poly_SOURCES})
  set_target_properties(static_poly PROPERTIES OUTPUT_NAME poly)
  install(TARGETS static_poly ARCHIVE DESTINATION lib)

  add_library(static_polyxx STATIC ${polyxx_SOURCES})
  set_target_properties(static_polyxx PROPERTIES OUTPUT_NAME polyxx)
  install(TARGETS static_polyxx ARCHIVE DESTINATION lib)
  target_link_libraries(static_polyxx static_poly)
endif()

#
# If libpoly.a is not position independent, we may need
# another version that's build to be position-indenpedent code.
# This allows the library to be embedded into something like libyices.so
#
# Linux, FreeBSD, and old versions of Darwin (and probably other systems) don't use fPIC by default.
# We just check for Linux and FreeBSD here.
# 
if(LIBPOLY_BUILD_STATIC_PIC AND ${CMAKE_SYSTEM_NAME} MATCHES "Linux" OR ${CMAKE_SYSTEM_NAME} MATCHES "FreeBSD")
  add_library(static_pic_poly STATIC ${poly_SOURCES})
  set_target_properties(static_pic_poly PROPERTIES OUTPUT_NAME picpoly POSITION_INDEPENDENT_CODE true)

  add_library(static_pic_polyxx STATIC ${polyxx_SOURCES})
  set_target_properties(static_pic_polyxx PROPERTIES OUTPUT_NAME picpolyxx POSITION_INDEPENDENT_CODE true)
  target_link_libraries(static_pic_polyxx static_pic_poly)
endif()
