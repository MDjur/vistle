vistle_find_package(PROJ)
if (NOT PROJ_FOUND)
  message("Skipping because of missing PROJ")
  return()
endif()

include_directories(SYSTEM ${PROJ_INCLUDE_DIR})
add_definitions(-DACCEPT_USE_OF_DEPRECATED_PROJ_API_H)

set(HEADERS
  MapDrape.h
)

set(SOURCES
  MapDrape.cpp
)

add_module(MapDrape ${SOURCES} ${HEADERS})
target_link_libraries(MapDrape ${PROJ_LIBRARIES})

