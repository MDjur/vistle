vistle_find_package(NETCDF)
if (NOT NETCDF_FOUND)
    return()
endif()

set(HEADERS
    ReadTsunami.h
)

set(SOURCES
    ReadTsunami.cpp
)

add_module(ReadTsunami ${SOURCES} ${HEADERS})

target_link_libraries(ReadTsunami ${NETCDF_LIBRARIES})
target_include_directories(ReadTsunami SYSTEM PRIVATE ${NETCDF_INCLUDE_DIRS})

