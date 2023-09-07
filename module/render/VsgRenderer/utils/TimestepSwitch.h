#ifndef TIMESTEPSWITCH_H
#define TIMESTEPSWITCH_H

#include <vsg/all.h>
#include <vsg/core/ref_ptr.h>

/**
 * Implementation of vsg::Switch to handle timestep animation.
 *
 * TimestepSwitch is an extension of vsg::Switch to handle timestep animation.
 *
 */
class TimestepSwitch: public vsg::Inherit<vsg::Switch, TimestepSwitch> {
public:
    TimestepSwitch(int in_stepWith = 1, int in_firstTimestep = 0)
    : m_stepWith(in_stepWith), m_currentTimestep(in_firstTimestep)
    {}

    auto getNumTimesteps() const { return children.size(); }
    auto getCurrentTimestep() const { return m_currentTimestep; }
    auto getStepWith() const { return m_stepWith; }

    /**
     * @brief Set the current timestep and activate the corresponding children.
     */
    void setCurrentTimestep(int timestep);

    /**
     * @brief Add a child to the switch which holds meta data "timestep". E.g. node->setValue("timestep", timestep).
     * The child will be inserted at the position timestep + current block if timestep is available.
     *
     * @param child the child to add with meta data "timestep".
     */
    void addChild(vsg::ref_ptr<vsg::Node> child);

    /**
     * @brief Remove a child.
     * The child will be removed with std::erase by checking timestep.
     *
     * @param child the child to add with meta data "timestep".
     */
    void removeChild(int step);

    /**
     * @brief Traverse the switch and enable all children for the current timestep.
     *
     * @return true if the current timestep is valid, false otherwise.
     */
    bool traverseTime();

private:
    /**
     * @brief Enable all childs for current timestep and disable all others.
     */
    void enableTimestepChilds();
    int m_stepWith;
    int m_currentTimestep;

protected:
    virtual ~TimestepSwitch() {
    } // hide from shared_ptr und explicit deletion. E.g. use "vsg::ref_ptr<TimestepSwitch> ts_ptr = TimestepSwitch::create();""
    // => "std::shared_ptr<TimestepSwitch> ts_ptr = std::make_shared<TimestepSwitch>(...)" or "delete ts_ptr" will no longer compile
};
#endif
