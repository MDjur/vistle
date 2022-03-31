#ifndef CLIP_H
#define CLIP_H

#include <vistle/module/module.h>
#include <vistle/core/vector.h>
#include "../IsoSurface/IsoDataFunctor.h"

class Clip: public vistle::Module {
public:
    Clip(const std::string &name, int moduleID, mpi::communicator comm);
    ~Clip();

    vistle::Object::ptr clip(vistle::Object::const_ptr object) const;

private:
    bool compute(std::shared_ptr<vistle::BlockTask> task) const override;
    virtual bool changeParameter(const vistle::Parameter *param) override;
    IsoController isocontrol;
};

#endif
