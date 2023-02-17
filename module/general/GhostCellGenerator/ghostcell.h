#ifndef GHOSTCELLS_H
#define GHOSTCELLS_H

#include "vistle/core/index.h"
#include "vistle/core/object.h"
#include "vistle/core/scalar.h"
#include "vistle/core/unstr.h"

#include <boost/mpi/communicator.hpp>
#include <boost/serialization/access.hpp>
#include <boost/smart_ptr/make_shared_array.hpp>
#include <boost/smart_ptr/shared_ptr.hpp>
#include <vector>

namespace vistle {

class GhostCell {
    friend class boost::serialization::access;

public:
    GhostCell(vistle::Index el, int timestep, std::vector<vistle::Scalar> x, std::vector<vistle::Scalar> y,
              std::vector<vistle::Scalar> z, std::vector<vistle::Index> connectivityList);

    vistle::Index el_id() const { return m_el; };
    std::vector<vistle::Scalar> X() { return m_x; };
    std::vector<vistle::Scalar> Y() { return m_y; };
    std::vector<vistle::Scalar> Z() { return m_z; };
    std::vector<vistle::Index> connectivityList() { return m_connectivityList; };

    template<class Archive>
    void serialize(Archive &ar, const unsigned int version)
    {
        ar &m_el;
        ar &m_x;
        ar &m_y;
        ar &m_z;
        ar &m_connectivityList;
    }

private:
    int m_timestep;

    // index of cell of current ghostcell
    vistle::Index m_el;

    // vertices
    std::vector<vistle::Scalar> m_x;
    std::vector<vistle::Scalar> m_y;
    std::vector<vistle::Scalar> m_z;
    std::vector<vistle::Index> m_connectivityList;
};

class GhostZones {
    friend class boost::serialization::access;
    friend class GhostCell;

    GhostZones(): m_block(vistle::InvalidIndex) {}
    GhostZones(vistle::Index block, int rank, int layer, const std::vector<Index> &domainCells,
               UnstructuredGrid::const_ptr ugrid);

    int rank() { return m_rank; };
    int numLayer() { return m_numLayer; };
    std::vector<GhostCell> &ghostzone() { return m_ghostcells; };

    template<class Archive>
    void serialize(Archive &ar, const unsigned int version)
    {
        ar &m_block;
        ar &m_ghostcells;
    }

    void clear()
    {
        m_block = vistle::InvalidIndex;
        m_ghostcells.clear();
    }

private:
    vistle::Index m_block;
    std::vector<GhostCell> m_ghostcells;
    int m_rank;
    int m_numLayer;
};

}; // namespace vistle

#endif
