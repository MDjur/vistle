if (WIN32)
   set(SPAWN_VISTLE spawn_vistle.bat)
else()
   set(SPAWN_VISTLE spawn_vistle.sh)
endif()

add_custom_target(vistle_spawner ALL
   COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_CURRENT_SOURCE_DIR}/${SPAWN_VISTLE} ${EXECUTABLE_OUTPUT_PATH}/${SPAWN_VISTLE}
   DEPENDS ${SPAWN_VISTLE}
   SOURCES ${SPAWN_VISTLE}
)

if (UNIX)
add_custom_target(vistle_tui ALL
   COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_CURRENT_SOURCE_DIR}/vistle_tui.sh ${EXECUTABLE_OUTPUT_PATH}/vistle_tui
   DEPENDS vistle_tui.sh
   SOURCES vistle_tui.sh
)
endif()

install(PROGRAMS ${SPAWN_VISTLE} DESTINATION bin)
