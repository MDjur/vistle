#ifndef TIMESTEPSWITCH_H
#define TIMESTEPSWITCH_H

#include <vsg/all.h>

namespace vistle {
class ExtendedVSGBuilder: public vsg::Inherit<vsg::Builder, ExtendedVSGBuilder> {
public:
    vsg::ref_ptr<vsg::Node> createTriangles(const vsg::GeometryInfo &info = {}, const vsg::StateInfo &stateInfo = {});

protected:
    GeometryMap _triangles;
};
} // namespace vistle

#endif
