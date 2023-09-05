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
    TimestepSwitch(int _timestep): m_currentTimestep(_timestep) {}
    void addChild(vsg::ref_ptr<vsg::Node> child)
    {
        int timestep;
        if (child->getValue("timestep", timestep)) {
            children.insert(children.begin() + timestep, Child{vsg::boolToMask(false), child});
            setTimestep(timestep);
        }
    }
    int getCurrentTimestep() const { return m_currentTimestep; }

private:
    void setTimestep(int _timstep) { m_currentTimestep = _timstep; }
    int m_currentTimestep;

protected:
    virtual ~TimestepSwitch() {}
};
#endif
