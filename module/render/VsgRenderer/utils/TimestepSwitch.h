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
    TimestepSwitch(int in_numTimesteps = 1, int in_stepWith = 1, int in_firstTimestep = -1)
    : m_numTimesteps(in_numTimesteps), m_stepWith(in_stepWith), m_currentTimestep(in_firstTimestep)
    {}

    auto getCurrentTimestep() const { return m_currentTimestep; }
    auto numTimesteps() { return m_numTimesteps; }
    auto numTimesteps() const { return m_numTimesteps; }
    auto stepWith() const { return m_stepWith; }

    void addChild(vsg::ref_ptr<vsg::Node> child)
    {
        int timestep;
        if (child->getValue("timestep", timestep)) {
            children.insert(children.begin() + timestep, Child{vsg::boolToMask(false), child});
            setTimestep(timestep);
        }
    }

    bool traverseTime()
    {
        auto numAllBlocks = children.size();
        if (numAllBlocks == 0)
            return true;
        m_currentTimestep += m_stepWith;
        if (m_currentTimestep > m_numTimesteps)
            // reset timestep
            if (!children[0].node->getValue("timestep", m_currentTimestep))
                return false;

        setAllChildren(false);

        // enable all blocks for current timestep
        auto itr = children.begin() + m_currentTimestep * numAllBlocks / m_numTimesteps;
        int blockTimestep;
        while (itr->node->getValue("timestep", blockTimestep) && blockTimestep == m_currentTimestep)
            (*itr++).mask = vsg::boolToMask(true);

        return true;
    }


private:
    void setTimestep(int _timstep) { m_currentTimestep = _timstep; }
    int m_numTimesteps;
    int m_stepWith;
    int m_currentTimestep;

protected:
    virtual ~TimestepSwitch() {}
};
#endif
