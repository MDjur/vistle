//@brief TimestepSwitch is an extension of vsg::Switch to handle the different timesteps.
//It is used by the VSGTimestepHandler.
//
// Created by Marko Djuric on 25.08.23.
//

#ifndef TIMESTEPSWITCH_H
#define TIMESTEPSWITCH_H

#include <vsg/all.h>

class TimestepSwitch: public vsg::Inherit<vsg::Switch, TimestepSwitch> {
public:
    TimestepSwitch(int _timestep): timestep(_timestep) {}
    void addChild(int timestep, vsg::ref_ptr<vsg::Node> child)
    {
        children.insert(children.begin() + timestep, Child{vsg::boolToMask(false), child});
        setTimestep(timestep);
    }
    int getTimestep() const { return timestep; }

private:
    void setTimestep(int _timstep) { timestep = _timstep; }
    int timestep;

protected:
    virtual ~TimestepSwitch() {}
};
#endif
