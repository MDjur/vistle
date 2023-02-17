if(${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated")
elseif(${CMAKE_CXX_COMPILER_ID} MATCHES "Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated -Wno-inconsistent-missing-override")
endif()

set(HEADERS GhostCellGenerator.h ghostcell.h)
set(SOURCES GhostCellGenerator.cpp ghostcell.cpp)

add_module(GhostCellGenerator "Generate Ghost Cells for UnstructuredGrids" ${SOURCES} ${HEADERS})
