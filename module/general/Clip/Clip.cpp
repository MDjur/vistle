#include <sstream>
#include <iomanip>

#include <vistle/core/object.h>
#include <vistle/core/polygons.h>
#include <vistle/core/triangles.h>
#include <vistle/util/math.h>

#include "Clip.h"
#include "PlaneClip.h"
#include "../IsoSurface/IsoDataFunctor.h"

MODULE_MAIN(Clip)

using namespace vistle;

Clip::Clip(const std::string &name, int moduleID, mpi::communicator comm)
: Module(name, moduleID, comm), isocontrol(this)
{
    isocontrol.init();

    setDefaultCacheMode(ObjectCache::CacheDeleteLate);

    createInputPort("grid_in");
    createOutputPort("grid_out");
}

Clip::~Clip()
{}

Object::ptr Clip::clip(Object::const_ptr object) const
{
    auto coords = Coords::as(object);
    if (!coords)
        return Object::ptr();

    IsoDataFunctor decider =
        isocontrol.newFunc(object->getTransform(), &coords->x()[0], &coords->y()[0], &coords->z()[0]);

    switch (object->getType()) {
    case Object::TRIANGLES: {
        PlaneClip cutter(Triangles::as(object), decider);
        cutter.process();
        return cutter.result();
    }

    case Object::POLYGONS: {
        PlaneClip cutter(Polygons::as(object), decider);
        cutter.process();
        return cutter.result();
    }

    default:
        break;
    }
    return Object::ptr();
}

bool Clip::changeParameter(const Parameter *param)
{
    bool ok = isocontrol.changeParameter(param);
    return Module::changeParameter(param) && ok;
}

bool Clip::compute(std::shared_ptr<BlockTask> task) const
{
    Object::const_ptr oin = task->expect<Object>("grid_in");
    if (!oin)
        return false;

    Object::ptr object = clip(oin);
    if (object) {
        object->copyAttributes(oin);
        updateMeta(object);
        task->addObject("grid_out", object);
    }

    return true;
}
