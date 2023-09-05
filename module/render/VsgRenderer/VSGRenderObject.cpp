#include "VSGRenderObject.h"
#include <sstream>
#include <vistle/core/triangles.h>
#include <vsg/core/Data.h>

using namespace vistle;

VsgRenderObject::VsgRenderObject(int senderId, const std::string &senderPort, vistle::Object::const_ptr container,
                                 vistle::Object::const_ptr geometry, vistle::Object::const_ptr normals,
                                 vistle::Object::const_ptr texture)
: vistle::RenderObject(senderId, senderPort, container, geometry, normals, texture)
{
    updateBounds();
    std::stringstream debug;

    m_geometry = vsg::Geometry::create();
    /* auto vertices = vsg::vec3Array::create(); */
    auto texcoords = vsg::vec2Array::create();
    /* auto colors = vsg::vec4Array::create(); */
    /* auto indices = vsg::ushortArray::create(); */
    vsg::ref_ptr<vsg::vec3Array> vertices;
    /* vsg::ref_ptr<vsg::vec2Array> texcoords; */
    vsg::ref_ptr<vsg::vec4Array> colors;
    vsg::ref_ptr<vsg::ushortArray> indices;

    // fill texture coordinates
    if (auto t = Texture1D::as(mapdata)) {
        auto coord_ptr = t->coords();
        for (size_t i = 0; i < t->getNumCoords(); ++i)
            texcoords->set(i, vsg::vec2(coord_ptr[i], 0.0f));

        debug << texcoords;
    }
    /* } else if (auto s = Vec<Scalar, 1>::as(this->mapdata)) { // check for scalar data 1D (Velocity in x) */
    /*     data->texCoords = &s->x()[0]; */

    /*     std::cerr << "texcoords from scalar field" << std::endl; */

    /* } else if (auto vec = Vec<Scalar, 3>::as(this->mapdata)) { // check for scalar data 3D (Velocity) */
    /*     tcoord.resize(vec->getSize()); */
    /*     data->texCoords = tcoord.data(); */
    /*     const Scalar *x = &vec->x()[0]; */
    /*     const Scalar *y = &vec->y()[0]; */
    /*     const Scalar *z = &vec->z()[0]; */
    /*     for (auto it = tcoord.begin(); it != tcoord.end(); ++it) { */
    /*         *it = sqrtf(*x * *x + *y * *y + *z * *z); */
    /*         ++x; */
    /*         ++y; */
    /*         ++z; */
    /*     } */
    /* } else if (auto iscal = Vec<Index>::as(this->mapdata)) { // check if mapped data is only integer type (partition) */
    /*     tcoord.resize(iscal->getSize()); */
    /*     data->texCoords = tcoord.data(); */
    /*     vistle::Scalar *d = tcoord.data(); */
    /*     for (const Index *i = &iscal->x()[0], *end = i + iscal->getSize(); i < end; ++i) { */
    /*         *d++ = *i; */
    /*     } */
    /* } */
    // convert vistle geometry to vsg geometry
    if (auto tri = Triangles::as(geometry)) {
        auto numElem = tri->getNumElements();
        auto numCoords = tri->getNumCoords();
        auto numCorners = tri->getNumCorners();

        vertices = vsg::vec3Array::create(numCoords);
        indices = vsg::ushortArray::create(numElem * 3);

        //unref after transfer to gpu
        vertices->properties.dataVariance = vsg::STATIC_DATA_UNREF_AFTER_TRANSFER;
        indices->properties.dataVariance = vsg::STATIC_DATA_UNREF_AFTER_TRANSFER;
        auto triX = tri->x();
        auto triY = tri->y();
        auto triZ = tri->z();

        for (auto i = 0; i < numCoords; ++i)
            vertices->set(i, vsg::vec3(triX[i], triY[i], triZ[i]));

        auto triConnectivityList = tri->cl();
        for (auto i = 0; i < numElem; ++i) {
            for (auto j = 0; j < 3; ++j) {
                uint32_t idx = j + i * 3;
                indices->at(idx) = triConnectivityList[idx];
            }
        }

        debug << "Tri: #: " << numElem << ", #corners: " << numCorners << ", #coord: " << numCoords << std::endl;
    }

    m_geometry->assignArrays(vsg::DataList{vertices, colors, texcoords});
    m_geometry->assignIndices(indices);
    m_geometry->commands.push_back(vsg::DrawIndexed::create(indices->size(), 1, 0, 0, 0));

    std::cerr << debug.str() << std::endl;
}
