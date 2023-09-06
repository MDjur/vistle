
#ifndef VSGTIMESTEPHANDLER_H
#define VSGTIMESTEPHANDLER_H

#include "utils/TimestepSwitch.h"
#include <vistle/core/vectortypes.h>
#include <vsg/all.h>

/** 
 * Implementation of vsg::Object for managing timesteps.
 *
 * VSGTimestepHandler is a class to handle the different timesteps of the
 * scenegraph. It is used by the VsgRenderer and the VsgGeometryGenerator.
 *
 */
class VSGTimestepHandler: public vsg::Inherit<vsg::Object, VSGTimestepHandler> {
public:
    explicit VSGTimestepHandler(vsg::ref_ptr<vsg::Viewer> in_viewer, int in_numTimesteps, int in_stepWith,
                                int in_firstTimestep);

    /**
     * @brief Add a geometry to root, either to the static (step = -1) or animated part (step > -1).
     *
     * @param geo the geometry to add
     * @param step the timestep of geo
    */
    void addNode(vsg::ref_ptr<vsg::Node> geo, const int step);

    /**
     * @brief Remove node from root.
     *
     * @param geo the geometry to remove
     * @param step the timestep of geo
    */
    void removeNode(vsg::ref_ptr<vsg::Node> geo, const int step);
    vsg::ref_ptr<vsg::MatrixTransform> root() const { return m_root; }
    vsg::ref_ptr<TimestepSwitch> animated() const { return m_animated; }

private:
    template<typename VSGNodeType>
    void addThreadSafe(vsg::ref_ptr<VSGNodeType> attachTo, vsg::ref_ptr<vsg::Node> node);

    vsg::ref_ptr<vsg::MatrixTransform> m_root;
    vsg::ref_ptr<vsg::Group> m_fixed;
    vsg::ref_ptr<TimestepSwitch> m_animated;
    vsg::observer_ptr<vsg::Viewer> m_viewer;

protected:
    virtual ~VSGTimestepHandler() {}
};
#endif
