vistle_find_package(vsg)
vistle_find_package(vsgXchange)

add_module(
    VsgRenderer
    "VulkanSceneGraph remote renderer"
    VSGRenderer.cpp
    VSGRenderer.h
    VSGRenderObject.h
    VSGRenderObject.cpp
    VSGTimestepHandler.cpp
    VSGTimestepHandler.h
    TimestepSwitch.h)

target_link_libraries(VsgRenderer vsg::vsg vistle_module vistle_renderer vistle_rhr)

if(vsgXchange_FOUND)
    target_compile_definitions(VsgRenderer PRIVATE vsgXchange_FOUND)
    target_link_libraries(VsgRenderer vsgXchange::vsgXchange)
endif()
