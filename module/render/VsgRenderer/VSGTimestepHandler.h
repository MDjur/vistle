//@brief VSGTimestepHandler is a class to handle the different timesteps of the
//scenegraph. It is used by the VsgRenderer and the VsgGeometryGenerator.

// Created by Marko Djuric on 24.08.23.
//

#ifndef VSGTIMESTEPHANDLER_H
#define VSGTIMESTEPHANDLER_H

#include <vistle/core/vectortypes.h>
#include <vsg/all.h>
#include "TimestepSwitch.h"

class VSGTimestepHandler: public vsg::Inherit<vsg::Object, VSGTimestepHandler> {
public:
    explicit VSGTimestepHandler();

    /*@brief Add a geometry to the scenegraph, either to the static or animated part.
     *
     * @param geo the geometry to add
     * @param step the timestep to add the geometry to
    */
    void addVSGObject(vsg::ref_ptr<vsg::Node> geo, const int step);
    vsg::ref_ptr<vsg::MatrixTransform> root() const { return m_root; }
    vsg::ref_ptr<TimestepSwitch> animated() const { return m_animated; }

    bool setTimestep(const int timestep);

private:
    vsg::ref_ptr<vsg::MatrixTransform> m_root;
    vsg::ref_ptr<vsg::Group> m_fixed;
    vsg::ref_ptr<TimestepSwitch> m_animated;

protected:
    virtual ~VSGTimestepHandler() {}
};
#endif
