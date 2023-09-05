#include "VSGTimestepHandler.h"
#include "utils/operations.h"
#include <vsg/app/Viewer.h>

VSGTimestepHandler::VSGTimestepHandler(vsg::ref_ptr<vsg::Viewer> in_viewer, int in_numTimesteps, int in_stepWith,
                                       int in_firstTimestep)
: m_viewer(in_viewer)
{
    m_root = vsg::MatrixTransform::create();
    m_fixed = vsg::Group::create();
    m_animated = TimestepSwitch::create(in_numTimesteps, in_stepWith, in_firstTimestep);

    m_root->addChild(m_fixed);
    m_root->addChild(m_animated);
}

void VSGTimestepHandler::addNode(vsg::ref_ptr<vsg::Node> geo, const int step)
{
    // add timestep as meta data
    geo->setValue("timestep", step);

    // timestep -1 is static
    if (step < 0) {
        //m_fixed->addChild(geo);
        addThreadSafe(m_fixed, geo);
    } else {
        //m_animated->addChild(geo);
        addThreadSafe(m_animated, geo);
    }
}

void VSGTimestepHandler::removeNode(vsg::ref_ptr<vsg::Node> geo, const int step)
{
    // timestep -1 is static
    if (step < 0)
        m_fixed->children.clear();
    else {
        //FIXME: the problem with this is that if you remove a timestep in the middle, the indices of the following timesteps are shifted
        m_animated->children.erase(m_animated->children.begin() + step);
        /* m_animated->children.shrink_to_fit(); */
        /* m_animated->children[step].node.reset(); */
        /* auto children = m_animated->children; */
        /* children.erase(std::remove(children.begin(), children.end(), geo), children.end()); */
    }
}

template<typename VSGGroupNodeType>
void VSGTimestepHandler::addThreadSafe(vsg::ref_ptr<VSGGroupNodeType> attachTo, vsg::ref_ptr<vsg::Node> node)
{
    vsg::ref_ptr<vsg::Viewer> ref_viewer = m_viewer;
    auto result = ref_viewer->compileManager->compile(node);
    if (result)
        ref_viewer->addUpdateOperation(Merge<VSGGroupNodeType>::create(m_viewer, attachTo, node, result));
}
