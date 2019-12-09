// fewer files compile faster
#ifndef TEMPLATES_IN_HEADERS

#define VISTLE_IMPL
#include "archives.h"
#include "vec.cpp"
#include "coords.cpp"
#include "coordswradius.cpp"
#include "normals.cpp"
#include "points.cpp"
#include "spheres.cpp"
#include "tubes.cpp"
#include "indexed.cpp"
#include "lines.cpp"
#include "ngons.cpp"
#include "polygons.cpp"
#include "texture1d.cpp"
#include "empty.cpp"
#include "placeholder.cpp"
#include "celltree.cpp"
#include "vertexownerlist.cpp"
#include "unstr.cpp"
#include "structuredgridbase.cpp"
#include "uniformgrid.cpp"
#include "rectilineargrid.cpp"
#include "structuredgrid.cpp"
#include "findobjectreferenceoarchive.cpp"

#else

#include "object.h"
#include "empty.h"
#include "archives.h"
#include "vec.h"
#include "coords.h"
#include "coordswradius.h"
#include "normals.h"
#include "points.h"
#include "spheres.h"
#include "tubes.h"
#include "indexed.h"
#include "lines.h"
#include "triangles.h"
#include "quads.h"
#include "polygons.h"
#include "texture1d.h"
#include "empty.h"
#include "placeholder.h"
#include "celltree.h"
#include "celltree_impl.h"
#include "vertexownerlist.h"
#include "unstr.h"
#include "structuredgridbase.h"
#include "uniformgrid.h"
#include "rectilineargrid.h"
#include "structuredgrid.h"
#include "findobjectreferenceoarchive.h"

#endif

namespace vistle {

#define REGISTER_TYPE(ObjType, id) \
do { \
   ObjectTypeRegistry::registerType<ObjType>(id); \
} while (true)

#define REGISTER_VEC_TYPE(t) \
do { \
   ObjectTypeRegistry::registerType<Vec<t,1>>(Vec<t,1>::type()); \
   ObjectTypeRegistry::registerType<Vec<t,3>>(Vec<t,3>::type()); \
} while (true)

void registerTypes() {

   using namespace vistle;
   REGISTER_TYPE(Empty, Object::EMPTY);
   REGISTER_TYPE(PlaceHolder, Object::PLACEHOLDER);
   REGISTER_TYPE(Texture1D, Object::TEXTURE1D);
   REGISTER_TYPE(Points, Object::POINTS);
   REGISTER_TYPE(Spheres, Object::SPHERES);
   REGISTER_TYPE(Lines, Object::LINES);
   REGISTER_TYPE(Tubes, Object::TUBES);
   REGISTER_TYPE(Triangles, Object::TRIANGLES);
   REGISTER_TYPE(Quads, Object::QUADS);
   REGISTER_TYPE(Polygons, Object::POLYGONS);
   REGISTER_TYPE(UniformGrid, Object::UNIFORMGRID);
   REGISTER_TYPE(RectilinearGrid, Object::RECTILINEARGRID);
   REGISTER_TYPE(StructuredGrid, Object::STRUCTUREDGRID);
   REGISTER_TYPE(UnstructuredGrid, Object::UNSTRUCTUREDGRID);
   REGISTER_TYPE(VertexOwnerList, Object::VERTEXOWNERLIST);
   REGISTER_TYPE(Celltree1, Object::CELLTREE1);
   REGISTER_TYPE(Celltree2, Object::CELLTREE2);
   REGISTER_TYPE(Celltree3, Object::CELLTREE3);
   REGISTER_TYPE(Normals, Object::NORMALS);

   REGISTER_VEC_TYPE(char);
   REGISTER_VEC_TYPE(signed char);
   REGISTER_VEC_TYPE(unsigned char);
   REGISTER_VEC_TYPE(int32_t);
   REGISTER_VEC_TYPE(uint32_t);
   REGISTER_VEC_TYPE(int64_t);
   REGISTER_VEC_TYPE(uint64_t);
   REGISTER_VEC_TYPE(float);
   REGISTER_VEC_TYPE(double);
}

} // namespace vistle
