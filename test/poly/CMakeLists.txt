set(poly_tests
    test_feasible_int_set
)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Werror -Wextra -std=c++11")

include_directories(${GMP_INCLUDE_DIR})

foreach(file ${poly_tests})
    add_executable(${file} ${file}.cpp)
    target_include_directories(${file} PUBLIC ${CMAKE_SOURCE_DIR}/include)
    add_test(NAME ${file} COMMAND ${file})
    target_link_libraries(${file} poly polyxx)
endforeach()
