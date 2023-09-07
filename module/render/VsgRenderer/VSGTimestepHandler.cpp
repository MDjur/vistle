#include "VSGTimestepHandler.h"
#include "utils/operations.h"

VSGTimestepHandler::VSGTimestepHandler(vsg::ref_ptr<vsg::Viewer> in_viewer, int in_stepWith, int in_firstTimestep)
: m_viewer(in_viewer)
{
    m_root = vsg::MatrixTransform::create();
    m_fixed = vsg::Group::create();
    m_animated = TimestepSwitch::create(in_stepWith, in_firstTimestep);

    m_root->addChild(m_fixed);
    m_root->addChild(m_animated);
}

void VSGTimestepHandler::addNode(vsg::ref_ptr<vsg::Node> geo, const int step)
{
    // timestep -1 is static
    if (step < 0)
        addThreadSafe(m_fixed, geo);
    else {
        vsg::ref_ptr<vsg::Group> timestepGroup;
        if (m_animated->getNumTimesteps() <= step) {
            timestepGroup = vsg::Group::create();
            // add timestep as meta data
            timestepGroup->setValue("timestep", step);
            timestepGroup->addChild(geo);
            addThreadSafe(m_animated, timestepGroup);
        } else {
            timestepGroup = m_animated->children[step].node->cast<vsg::Group>();
            addThreadSafe(timestepGroup, geo);
        }
    }
}

void VSGTimestepHandler::removeNode(vsg::ref_ptr<vsg::Node> geo, const int step)
{
    // timestep -1 is static
    if (step < 0)
        m_fixed->children.clear();
    else
        m_animated->removeChild(step);
}

template<typename VSGGroupNodeType>
void VSGTimestepHandler::addThreadSafe(vsg::ref_ptr<VSGGroupNodeType> attachTo, vsg::ref_ptr<vsg::Node> node)
{
    vsg::ref_ptr<vsg::Viewer> ref_viewer = m_viewer;
    auto result = ref_viewer->compileManager->compile(node);
    if (result)
        ref_viewer->addUpdateOperation(Merge<VSGGroupNodeType>::create(m_viewer, attachTo, node, result));
}
