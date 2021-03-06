set(control_SOURCES
   hub.cpp
   uimanager.cpp
   uiclient.cpp
   fileinfocrawler.cpp
   scanmodules.cpp
)

set(control_HEADERS
   export.h
   hub.h
   uimanager.h
   uiclient.h
   fileinfocrawler.h
   scanmodules.h
)

if(Python_FOUND)
   set(control_SOURCES ${control_SOURCES} ../userinterface/pythoninterface.cpp ../userinterface/pythonmodule.cpp pythoninterpreter.cpp)
   set(control_HEADERS ${control_HEADERS} ../userinterface/pythoninterface.h ../userinterface/pythonmodule.h pythoninterpreter.h)
endif()

vistle_add_library(vistle_control ${VISTLE_LIB_TYPE} ${control_SOURCES} ${control_HEADERS})
vistle_export_library(vistle_control ${VISTLE_LIB_TYPE} ${control_SOURCES} ${control_HEADERS})
target_link_libraries(vistle_control
    PRIVATE Boost::system
    PRIVATE Boost::program_options
    PRIVATE vistle_util
    PRIVATE vistle_core
    PRIVATE vistle_net
    PRIVATE vistle_libsim_connect
    PRIVATE Threads::Threads
)

if(Python_FOUND)
   target_compile_definitions(vistle_control
    PRIVATE -DHAVE_PYTHON
    PRIVATE -DVISTLE_CONTROL
    PRIVATE -Dvistle_userinterface_EXPORTS
    PRIVATE -DEMBED_PYTHON)
   target_link_libraries(vistle_control
    PRIVATE Python::Python)
endif()


if(UNIX)
    target_link_libraries(vistle_control PRIVATE ${CMAKE_DL_LIBS})
endif()

add_custom_target(vistle_python ALL
   COMMAND ${CMAKE_COMMAND} -E copy_if_different
   ${CMAKE_CURRENT_SOURCE_DIR}/vistle.py
   ${PROJECT_BINARY_DIR}/share/vistle/vistle.py
   DEPENDS vistle.py
   SOURCES vistle.py
)

install(FILES vistle.py DESTINATION share/vistle)
