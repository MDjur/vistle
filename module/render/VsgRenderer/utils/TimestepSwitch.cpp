#include "TimestepSwitch.h"

void TimestepSwitch::addChild(vsg::ref_ptr<vsg::Node> child)
{
    int timestep;
    if (child->getValue("timestep", timestep))
        children.insert(children.begin() + timestep, Child{vsg::boolToMask(false), child});
}

void TimestepSwitch::removeChild(int step)
{
    if (auto group = children[step].node->cast<vsg::Group>())
        if (group)
            group->children.clear();

    children.erase(std::remove_if(children.begin(), children.end(),
                                  [step](const Child &c) {
                                      int timestep;
                                      return c.node->getValue("timestep", timestep) == step;
                                  }),
                   children.end());
}

void TimestepSwitch::setCurrentTimestep(int timestep)
{
    m_currentTimestep = timestep;
    enableTimestepChilds();
}

void TimestepSwitch::enableTimestepChilds()
{
    setAllChildren(false);

    // enable all blocks for current timestep
    children[m_currentTimestep].mask = vsg::boolToMask(true);
}

bool TimestepSwitch::traverseTime()
{
    auto numTimesteps = children.size();
    if (numTimesteps == 0)
        return true;
    m_currentTimestep += m_stepWith;
    if (m_currentTimestep > numTimesteps)
        // reset timestep
        if (!children[0].node->getValue("timestep", m_currentTimestep))
            return false;

    enableTimestepChilds();
    return true;
}
