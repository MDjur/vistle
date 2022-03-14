if(WIN32)
    set(SPAWN_VISTLE spawn_vistle.bat)
else()
    set(SPAWN_VISTLE spawn_vistle.sh)
endif()

add_custom_target(
    vistle_spawner ALL
    COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_CURRENT_SOURCE_DIR}/${SPAWN_VISTLE} ${EXECUTABLE_OUTPUT_PATH}/${SPAWN_VISTLE}
    DEPENDS ${SPAWN_VISTLE}
    SOURCES ${SPAWN_VISTLE})

# if(UNIX)
#     add_custom_target(
#         vistle_tui ALL
#         COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_CURRENT_SOURCE_DIR}/vistle_tui.sh ${EXECUTABLE_OUTPUT_PATH}/vistle_tui
#         DEPENDS vistle_tui.sh
#         SOURCES vistle_tui.sh)
# endif()

if(UNIX)
    add_custom_target(
        vistle_opener ALL
        COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_CURRENT_SOURCE_DIR}/vistle_opener.sh ${EXECUTABLE_OUTPUT_PATH}/vistle_opener
        DEPENDS vistle_opener.sh
        SOURCES vistle_opener.sh)
endif()

install(PROGRAMS ${SPAWN_VISTLE} DESTINATION bin)
if(UNIX)
    # install(
    #     PROGRAMS vistle_tui.sh
    #     RENAME vistle_tui
    #     DESTINATION bin)
    install(
        PROGRAMS vistle_opener.sh
        RENAME vistle_opener
        DESTINATION bin)
endif()

find_program(PYINSTALLER pyinstaller)
if(NOT PYINSTALLER)
    return()
endif()

set(VISTLE_PYINSTALLER_COMPILE
    pyinstaller
    --noconfirm
    --console
    --distpath
    ${EXECUTABLE_OUTPUT_PATH}
    --onefile
    --name)

set(VISTLE_TUI_SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/vistle_tui_cmd.py)
set(VISTLE_TUI_BIN ${CMAKE_CURRENT_BINARY_DIR}/vistle_tui)

add_custom_target(build_vistle_tui ALL DEPENDS ${VISTLE_TUI_BIN})

add_custom_command(
    OUTPUT ${VISTLE_TUI_BIN}
    COMMENT "Create vistle_tui"
    COMMAND ${VISTLE_PYINSTALLER_COMPILE} vistle_tui ${VISTLE_TUI_SOURCE}
    COMMAND ${CMAKE_COMMAND} -E touch ${VISTLE_TUI_BIN}
    DEPENDS ${VISTLE_TUI_SOURCE}
    VERBATIM)
