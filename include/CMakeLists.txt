# Configure the version
configure_file(version.h.in ${CMAKE_CURRENT_SOURCE_DIR}/version.h)

# Find all the base headers and add them to install target
file(GLOB headers *.h)
install(FILES ${headers} DESTINATION include/poly)

# Find all the cxx headers and add them to install target
file(GLOB headers_cxx polyxx/*.h)
install(FILES ${headers_cxx} DESTINATION include/poly/polyxx)