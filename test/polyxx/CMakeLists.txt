set(polyxx_tests
    test_algebraic_number
    test_assignment
    test_context
    test_dyadic_interval
    test_dyadic_rational
    test_integer
    test_interval
    test_interval_assignment
    test_polynomial
    test_rational
    test_rational_interval
    test_upolynomial
    test_value
    test_variable
)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Werror -Wextra -std=c++11")

include_directories(${GMP_INCLUDE_DIR})

foreach(file ${polyxx_tests})
    add_executable(${file} ${file}.cpp)
    target_include_directories(${file} PUBLIC ${CMAKE_SOURCE_DIR}/include)
    add_test(NAME ${file} COMMAND ${file})
    target_link_libraries(${file} polyxx poly)
endforeach()
