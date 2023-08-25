#ifndef VSGRENDEROBJECT_H
#define VSGRENDEROBJECT_H
#include <vsg/all.h>
#include <vistle/renderer/renderer.h>


class VsgRenderObject: public vistle::RenderObject {
public:
    VsgRenderObject(int senderId, const std::string &senderPort, vistle::Object::const_ptr container,
                    vistle::Object::const_ptr geometry, vistle::Object::const_ptr normals,
                    vistle::Object::const_ptr texture, vsg::ref_ptr<vsg::Node> _node)
    : vistle::RenderObject(senderId, senderPort, container, geometry, normals, texture), node(_node)
    {}

    vsg::ref_ptr<vsg::Node> node;
};

#endif
