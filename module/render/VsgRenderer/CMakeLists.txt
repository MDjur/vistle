vistle_find_package(vsg)
vistle_find_package(vsgXchange)

set(utils_source utils/TimestepSwitch.cpp)
set(utils_header utils/TimestepSwitch.h utils/visitors.h utils/operations.h)

set(SOURCE VSGRenderer.cpp VSGRenderObject.cpp VSGTimestepHandler.cpp ${utils_source})

set(HEADERS VSGRenderer.h VSGRenderObject.h VSGTimestepHandler.h ${utils_header})

add_module(VsgRenderer "VulkanSceneGraph remote renderer" ${SOURCE} ${HEADERS})

target_link_libraries(VsgRenderer vsg::vsg vistle_module vistle_renderer vistle_rhr)

if(vsgXchange_FOUND)
    target_compile_definitions(VsgRenderer PRIVATE vsgXchange_FOUND)
    target_link_libraries(VsgRenderer vsgXchange::vsgXchange)
endif()
