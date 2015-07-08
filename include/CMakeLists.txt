# Configure the version
configure_file(version.h.in ${CMAKE_CURRENT_SOURCE_DIR}/version.h)

# Find all the headers and add them to install terget
file(GLOB headers *.h)
install(FILES ${headers} DESTINATION include/poly)