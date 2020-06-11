set(polyxx_tests
    test_variable
)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Werror -Wextra -std=c++11")

foreach(file ${polyxx_tests})
    add_executable(${file} ${file}.cpp)
    target_include_directories(${file} PUBLIC ${CMAKE_SOURCE_DIR}/include)
    add_test(NAME ${file} COMMAND ${file})
    target_link_libraries(${file} polyxx poly)
endforeach()
