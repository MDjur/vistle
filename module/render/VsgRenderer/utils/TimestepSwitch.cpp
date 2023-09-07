#include "TimestepSwitch.h"

void TimestepSwitch::addChild(vsg::ref_ptr<vsg::Node> child)
{
    static int currentBlock = 0;
    int timestep;
    if (child->getValue("timestep", timestep)) {
        children.insert(children.begin() + timestep + currentBlock, Child{vsg::boolToMask(false), child});
        // as long as timestep stays the same increment block counter
        if (timestep != m_currentTimestep) {
            ++m_numTimesteps;
            m_currentTimestep = timestep;
            if (m_blocksPerTimestep != currentBlock)
                m_blocksPerTimestep = currentBlock;
            currentBlock = -1;
        }
        ++currentBlock;
    }
}

void TimestepSwitch::removeChild(vsg::ref_ptr<vsg::Node> child)
{
    children.erase(std::remove_if(children.begin(), children.end(),
                                  [child](const Child &c) { return c.node.get() == child.get(); }),
                   children.end());
}

bool TimestepSwitch::traverseTime()
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
