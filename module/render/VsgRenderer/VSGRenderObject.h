#ifndef VSGRENDEROBJECT_H
#define VSGRENDEROBJECT_H
#include <vsg/all.h>
#include <vistle/renderer/renderer.h>

/**
 * Implementation of vistle::RenderObject to convert vistle::Object to vsg::Geometry.
 */
class VsgRenderObject: public vistle::RenderObject {
public:
    VsgRenderObject(int senderId, const std::string &senderPort, vistle::Object::const_ptr container,
                    vistle::Object::const_ptr geometry, vistle::Object::const_ptr normals,
                    vistle::Object::const_ptr texture);

    auto geo() const { return m_geometry; }

private:
    vsg::ref_ptr<vsg::Geometry> m_geometry;
};

#endif
