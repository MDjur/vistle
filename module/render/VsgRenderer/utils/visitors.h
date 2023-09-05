#ifndef VISITORS_H
#define VISITORS_H

#include <vsg/all.h>
#include "TimestepSwitch.h"

/* @brief Set the polygon mode of the rasterization state to line mode (wireframe).
*/
struct SetPipelineStates: public vsg::Inherit<vsg::Visitor, SetPipelineStates> {
    void apply(vsg::Object &object) override { object.traverse(*this); }
    void apply(vsg::RasterizationState &rs) override { rs.polygonMode = VK_POLYGON_MODE_LINE; }
};

class AnimateTimestepSwitch: public vsg::Inherit<vsg::Visitor, AnimateTimestepSwitch> {
public:
    AnimateTimestepSwitch(vsg::ref_ptr<TimestepSwitch> in_timestepSwitch, vsg::clock::time_point in_start)
    : m_timestepSwitch(in_timestepSwitch), start(in_start)
    {}

    vsg::observer_ptr<TimestepSwitch> m_timestepSwitch;
    vsg::clock::time_point start;
    double time = 0.0;

    void apply(vsg::FrameEvent &frame) override
    {
        vsg::ref_ptr<TimestepSwitch> ts_ptr = m_timestepSwitch;
        time = std::chrono::duration<double, std::chrono::seconds::period>(frame.time - start).count();
        if (time > ts_ptr->stepWith()) {
            ts_ptr->traverseTime();
            start = frame.time;
        }
        if (ts_ptr)
            ts_ptr->accept(*this);
    }
};

#endif
