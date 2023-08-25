#include "VSGTimestepHandler.h"

VSGTimestepHandler::VSGTimestepHandler()
{
    m_root = vsg::MatrixTransform::create();
    m_fixed = vsg::Group::create();
    m_animated = TimestepSwitch::create(0);

    m_root->addChild(m_fixed);
    m_root->addChild(m_animated);
}

void VSGTimestepHandler::addVSGObject(vsg::ref_ptr<vsg::Node> geo, const int step)
{
    // timestep -1 is static
    if (step < 0)
        m_fixed->addChild(geo);
    else
        m_animated->addChild(step, geo);
}

bool VSGTimestepHandler::setTimestep(const int timestep)
{
    if (!m_animated.valid() && m_animated->children.size() == 0)
        return false;

    m_animated->setSingleChildOn(timestep);
    return true;
}
