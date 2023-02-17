#include "ghostcell.h"
#include "vistle/core/index.h"
#include "vistle/core/unstr.h"
#include <algorithm>
#include <vector>

using namespace vistle;

GhostCell::GhostCell(vistle::Index el, int timestep, std::vector<vistle::Scalar> x, std::vector<vistle::Scalar> y,
                     std::vector<vistle::Scalar> z, std::vector<vistle::Index> connectivityList)
: m_timestep(timestep), m_el(el), m_x(x), m_y(y), m_z(z), m_connectivityList(connectivityList)

{}

GhostZones::GhostZones(vistle::Index block, int rank, int layer, const std::vector<Index> &domainCells,
                       UnstructuredGrid::const_ptr ugrid)
: m_block(block), m_rank(rank), m_numLayer(layer)
{
    const auto cellList = ugrid->el();
    const auto connectivityList = ugrid->cl();
    const auto x = ugrid->x();
    const auto y = ugrid->y();
    const auto z = ugrid->z();
    for (auto cellId: domainCells) {
        if (ugrid->tl()[cellId] == UnstructuredGrid::BAR)
            continue;
        auto connectivityStart = cellList[cellId];
        const auto connectivityEnd = cellList[cellId + 1];
        std::vector<vistle::Index> cellConnectivity;
        std::copy_n(connectivityList, connectivityEnd, cellConnectivity.begin());

        std::array<std::vector<vistle::Scalar>, 3> cellCoords;
        while (connectivityStart < connectivityEnd) {
            auto vertexId = connectivityList[connectivityStart++];
            cellCoords[0].push_back(x[vertexId]);
            cellCoords[1].push_back(y[vertexId]);
            cellCoords[2].push_back(z[vertexId]);
        }
        m_ghostcells.push_back(GhostCell(cellId, -1, cellCoords[0], cellCoords[1], cellCoords[2], cellConnectivity));
    }
}
