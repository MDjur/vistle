#ifndef VISITORS_H
#define VISITORS_H

#include <vsg/all.h>

/* @brief Set the polygon mode of the rasterization state to line mode (wireframe).
*/
struct SetPipelineStates: public vsg::Inherit<vsg::Visitor, SetPipelineStates> {
    void apply(vsg::Object &object) override { object.traverse(*this); }
    void apply(vsg::RasterizationState &rs) override { rs.polygonMode = VK_POLYGON_MODE_LINE; }
};

#endif
