#ifndef GHOSTCELLGENERATOR_H
#define GHOSTCELLGENERATOR_H

#include <vector>

#include <vistle/core/index.h>
#include <vistle/core/unstr.h>
#include <vistle/core/object.h>
#include <vistle/core/polygons.h>
#include <vistle/core/quads.h>
#include <vistle/core/structuredgridbase.h>
#include <vistle/module/module.h>

/* #include "ghostcell.h" */
namespace vistle {
class GhostCellGenerator: public Module {
public:
    GhostCellGenerator(const std::string &name, int moduleID, mpi::communicator comm);
    ~GhostCellGenerator();

    typedef std::vector<Index> DataMapping;

private:
    struct BoundaryCell {
        BoundaryCell(Index vertexId, Index cellId, Index faceId): vertexId(vertexId), cellId(cellId), faceId(faceId){};
        Index vertexId; // vertex which was used to pick this cell
        Index cellId;
        Index faceId;
    };
    bool prepare() override;
    bool reduce(int timestep) override;
    bool compute() override;
    void linkNeighbors(std::shared_ptr<mpi::communicator> neighbors) const;
    void addCellsContainingVert(std::vector<BoundaryCell> &boundary, UnstructuredGrid::const_ptr ugrid, Index lastCell,
                                Index vertexId, Index faceId);
    void addCellsContainingFaceVert(std::vector<BoundaryCell> &boundary, UnstructuredGrid::const_ptr ugrid,
                                    Index lastCell, Index faceId);
    std::vector<Index> blockdomainCells(UnstructuredGrid::const_ptr ugrid);
    BoundaryCell findFirstBoundaryCell(UnstructuredGrid::const_ptr ugrid);

    IntParameter *m_celltree;
    IntParameter *m_constCellSize;
    std::vector<std::vector<Index>> m_blockDomain;
};
} // namespace vistle
#endif
