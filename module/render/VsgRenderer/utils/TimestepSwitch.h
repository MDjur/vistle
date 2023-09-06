#ifndef TIMESTEPSWITCH_H
#define TIMESTEPSWITCH_H

#include <vsg/all.h>

/**
 * Implementation of vsg::Switch to handle timestep animation.
 *
 * TimestepSwitch is an extension of vsg::Switch to handle timestep animation.
 * It is used by the VSGTimestepHandler.
 *
 */
class TimestepSwitch: public vsg::Inherit<vsg::Switch, TimestepSwitch> {
public:
    TimestepSwitch(int in_numTimesteps = 1, int in_stepWith = 1, int in_firstTimestep = -1)
    : m_numTimesteps(in_numTimesteps), m_stepWith(in_stepWith), m_currentTimestep(in_firstTimestep)
    {}

    auto getCurrentTimestep() const { return m_currentTimestep; }
    auto numTimesteps() { return m_numTimesteps; }
    auto numTimesteps() const { return m_numTimesteps; }
    auto stepWith() const { return m_stepWith; }

    /**
     * @brief Add a child to the switch which holds meta data "timestep". E.g. node->setValue("timestep", timestep).
     * The child will be inserted at the position timestep if timestep is available.
     *
     * @param child the child to add with meta data "timestep".
     */
    void addChild(vsg::ref_ptr<vsg::Node> child)
    {
        int timestep;
        if (child->getValue("timestep", timestep)) {
            children.insert(children.begin() + timestep, Child{vsg::boolToMask(false), child});
            setTimestep(timestep);
        }
    }

    /**
     * @brief Traverse the switch and enable all children for the current timestep.
     *
     * @return true if the current timestep is valid, false otherwise.
     */
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
    virtual ~TimestepSwitch() {
    } // hide from shared_ptr und explicit deletion. E.g. use "vsg::ref_ptr<TimestepSwitch> ts_ptr = TimestepSwitch::create();""
    // => "std::shared_ptr<TimestepSwitch> ts_ptr = std::make_shared<TimestepSwitch>(...)" or "delete ts_ptr" will no longer compile
};
#endif
