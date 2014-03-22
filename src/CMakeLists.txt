add_subdirectory(interval)
add_subdirectory(number)
add_subdirectory(polynomial)
add_subdirectory(upolynomial)
add_subdirectory(utils)
add_subdirectory(variable)

add_library(poly 
    $<TARGET_OBJECTS:interval>
    $<TARGET_OBJECTS:number>
    $<TARGET_OBJECTS:polynomial>
    $<TARGET_OBJECTS:upolynomial>
    $<TARGET_OBJECTS:utils>
    $<TARGET_OBJECTS:variable>
)