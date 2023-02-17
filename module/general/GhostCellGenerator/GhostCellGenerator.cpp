/**************************************************************************\
 **                                                                      **
 **                                                                      **
 ** Description: GhostCellGenerator for unstructured grids.              **
 **                                                                      **
 ** Based on paper: Parallel Multi-Layer Ghost Cell Generation for       **
 **                 Distributed UnstructuredGrid Grids.                  **
 ** DOI: 10.1109/LDAV.2017.8231854                                       **
 **                                                                      **
 **                                                                      **

 TODO:
 [x] 1. Extract external boundary of local partition
    => Domain Surface.
 [ ] 2. Share boundary with potential partition neighbors
    => share with every partition
    [ ] optimize later.
 [ ] 3. Calculate actual partition neighbors
 [ ] 4. create cell list to send to each partition
 [ ] 5. send cells to partition neighbors.
 [ ] 6. receive cells from partition.
 [ ] 7. integrate cells into local partition.

 **                                                                      **
 **                                                                      **
 **                                                                      **
 ** Author:    Marko Djuric                                              **
 **                                                                      **
 **                                                                      **
 **                                                                      **
 ** Date:  10.05.2021                                                    **
\**************************************************************************/

// std
#include <boost/mpi/environment.hpp>
#include <memory>
#include <sstream>
#include <iomanip>
#include <vector>

// vistle
#include <vistle/core/object.h>
#include <vistle/core/vec.h>
#include <vistle/core/unstr.h>
#include <vistle/core/serialize.h>
#include <vistle/core/polygons.h>
#include <vistle/core/structuredgrid.h>

// boost
#include <boost/mpi.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/collectives/all_reduce.hpp>
#include <boost/mpi/collectives/all_gather.hpp>
#include <boost/mpi/collectives/broadcast.hpp>
#include <boost/smart_ptr/shared_ptr.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/map.hpp>

//header
#include "GhostCellGenerator.h"
#include "ghostcell.h"
#include "vistle/core/index.h"
#include "vistle/core/messages.h"

using namespace vistle;
namespace mpi = boost::mpi;
namespace {

typedef GhostCellGenerator::DataMapping DataMapping;

Object::const_ptr isNeighbor(Object::const_ptr send, Object::const_ptr recv)
{
    return send;
}

unsigned numEdges(UnstructuredGrid::const_ptr ugrid, Index elem)
{
    auto t = ugrid->tl()[elem] & UnstructuredGrid::TYPE_MASK;
    switch (t) {
    case UnstructuredGrid::NONE:
    case UnstructuredGrid::POINT:
        return 0;
    case UnstructuredGrid::BAR:
        return 1;
    case UnstructuredGrid::TRIANGLE:
        return 3;
    case UnstructuredGrid::QUAD:
        return 4;
    case UnstructuredGrid::TETRAHEDRON:
        return 4;
    case UnstructuredGrid::PYRAMID:
        return 8;
    case UnstructuredGrid::PRISM:
        return 9;
    case UnstructuredGrid::HEXAHEDRON:
        return 12;
    }

    return -1;
}
}; // namespace

GhostCellGenerator::GhostCellGenerator(const std::string &name, int moduleID, mpi::communicator comm)
: Module(name, moduleID, comm)
{
    createInputPort("data_in", "input grid");
    createOutputPort("data_out", "grid with added ghost/halo cells");
    m_celltree = addIntParameter("create_celltree", "create celltree", 0, Parameter::Boolean);
    m_constCellSize = addIntParameter("const_cellSize", "const cellSize", 1, Parameter::Boolean);

    // policies
    setReducePolicy(message::ReducePolicy::Locally);
    setSchedulingPolicy(message::SchedulingPolicy::Gang);
}

GhostCellGenerator::~GhostCellGenerator()
{}

void GhostCellGenerator::linkNeighbors(std::shared_ptr<mpi::communicator> neighbors) const
{
    auto world = comm();
    *neighbors = world.split(1); //TODO: function to identify neighbors
}

bool GhostCellGenerator::prepare()
{
    m_blockDomain.clear();

    return true;
}

bool GhostCellGenerator::reduce(int timestep)
{
    std::vector<std::vector<std::vector<Index>>> neighborDomains;
    mpi::all_gather(comm(), m_blockDomain, neighborDomains);
    for (int neighborRank = 0; neighborRank < neighborDomains.size(); ++neighborRank) {
        if (neighborRank == rank())
            continue;
        auto blockdomains = neighborDomains[neighborRank];
        for (int block = 0; block < blockdomains.size(); ++block)
            for (auto &potentialGhost: blockdomains[block]) {
                sendInfo("Block %d send from rank %d Index potential ghost: %d", block, neighborRank, potentialGhost);
            }
    }
    return true;
}

GhostCellGenerator::BoundaryCell GhostCellGenerator::findFirstBoundaryCell(UnstructuredGrid::const_ptr ugrid)
{
    auto numElem = ugrid->getNumElements();
    auto cellList = ugrid->el();
    auto connectivityList = ugrid->cl();
    for (Index cellId = 0; cellId < numElem; ++cellId) {
        auto &nf = ugrid->getNeighborFinder();
        auto connectivityStart = cellList[cellId];
        auto connectivityEnd = cellList[cellId + 1];

        auto type = ugrid->tl()[cellId];
        if (type == UnstructuredGrid::BAR)
            break;

        // const auto &faces = UnstructuredGrid::FaceVertices[type];
        // const auto &sizes = UnstructuredGrid::FaceSizes[type];

        /* for (Index faceId = 0; faceId < sizes; ) */

        // find obvious face neighbor in cell
        while (connectivityStart < connectivityEnd) {
            auto vertexId1 = connectivityList[connectivityStart++];
            auto vertexId2 = connectivityList[connectivityStart++];
            auto vertexId3 = connectivityList[connectivityStart++];
            auto faceNeighbor = nf.getNeighborElement(cellId, vertexId1, vertexId2, vertexId3);
            if (faceNeighbor == InvalidIndex) {
                Index faceId = connectivityEnd - connectivityStart;
                return BoundaryCell(cellId, vertexId1, faceId);
            }
        }
    }
    return BoundaryCell(InvalidIndex, InvalidIndex, InvalidIndex);
}

void GhostCellGenerator::addCellsContainingFaceVert(std::vector<BoundaryCell> &boundary,
                                                    UnstructuredGrid::const_ptr ugrid, Index lastCell, Index faceId)
{
    auto type = ugrid->tl()[lastCell] & UnstructuredGrid::TYPE_MASK;
    auto connectivityList = ugrid->cl();
    auto cellList = ugrid->el();
    const auto &connectivityStart = cellList[lastCell];
    const auto &faces = UnstructuredGrid::FaceVertices[type];
    const auto &sizes = UnstructuredGrid::FaceSizes[type];
    const auto &face = faces[faceId];
    int N = sizes[faceId];
    for (int i = 0; i < N; ++i) {
        Index vertexFace = face[i];
        Index vertexId = connectivityList[connectivityStart + vertexFace];
        addCellsContainingVert(boundary, ugrid, lastCell, vertexId, faceId);
    }
}

void GhostCellGenerator::addCellsContainingVert(std::vector<BoundaryCell> &boundary, UnstructuredGrid::const_ptr ugrid,
                                                Index lastCell, Index vertexId, Index faceId)
{
    auto &nf = ugrid->getNeighborFinder();
    auto cells = nf.getContainingElements(vertexId);
    for (auto &cell: cells) {
        if (cell == lastCell)
            continue;
        boundary.push_back(BoundaryCell(vertexId, cell, faceId));
        lastCell = cell;
    }
}

DataMapping GhostCellGenerator::blockdomainCells(UnstructuredGrid::const_ptr ugrid)
{
    std::vector<Index> domain;
    auto numElem = ugrid->getNumElements();
    auto cellList = ugrid->el();
    for (Index elem = 0; elem < numElem; ++elem) {
        auto neighbor = ugrid->getNeighborElements(elem);
        auto numEdge = numEdges(ugrid, elem);
        const auto numFaces = ugrid->cellNumFaces(elem);
        const auto numCorners = cellList[elem + 1] - cellList[elem];
        auto minNumNeighbor = numEdge + numFaces + numCorners;

        if (neighbor.size() < minNumNeighbor)
            domain.push_back(elem);
    }
    return domain;
}

bool GhostCellGenerator::compute()
{
    /*
     * 1. gather all Block vec<obj_const_ptr> per timestep (compute)
     * 2. send vec to MPI_RANKS (reduce)
     * 3. go through each received vec and look for neighbors and safe rank and cellIndex (reduce)
     * 4. send neighbor info (reduce)
     * 5. attach ghostcells (reduce)
     * 6. pass to outputport
     */
    DataBase::const_ptr data;
    UnstructuredGrid::const_ptr unstr = accept<UnstructuredGrid>("data_in");
    if (!unstr) {
        data = expect<DataBase>("data_in");
        if (!data) {
            sendError("no grid and no data received");
            return true;
        }
        if (!data->grid()) {
            sendError("no grid attached to data");
            return true;
        }
        unstr = UnstructuredGrid::as(data->grid());
        if (!unstr) {
            sendError("no valid grid attached to data");
            return true;
        }
    }

    //init celltree
    if (m_celltree->getValue())
        if (!unstr->hasCelltree())
            unstr->getCelltree();

    //init vertexownerlist
    unstr->getNeighborElements(InvalidIndex);

    if (m_blockDomain.empty())
        m_blockDomain.resize(unstr->meta().numBlocks());

    const auto &block = unstr->meta().block();
    m_blockDomain[block] = blockdomainCells(unstr);

    /* Object::const_ptr surface_in = surface; */
    /* assert(surface_in); */

    /* const auto &block = surface_in->meta().block(); */
    /* const auto &numBlock = surface_in->meta().numBlocks(); */
    /* const auto &procRank = rank(); */

    /* int nblocks = surface_in->meta().numBlocks(); */
    /* int threads = hardware_concurrency(); */
    /* int in_memory = 1; */

    // TODO: DomainSurface implementation here => for now needs domainsurface linked before
    // gather surfaces from other ranks
    /* auto v_comm = comm(); */
    /* std::vector<Object::const_ptr> surfaces_recv(size()); */
    /* std::vector<std::string> test_vec(size()); */
    /* mpi::gather(v_comm, surface_in, surfaces_recv, rank()); */
    /* std::stringstream out; */
    /* out << "rank " << rank() << " thread: " */
    /*     << "\n"; */
    /* out.clear(); */
    /* mpi::gather(v_comm, out.str(), test_vec, rank()); */
    /* out << "rank: " << rank() << " size of vec:" << surfaces_recv.size() << "\n"; */
    /* sendInfo(out.str()); */

    /* std::shared_ptr<mpi::communicator> neighbors; */
    /* linkNeighbors(neighbors); */

#if 0
    bool haveElementData = false;
    if (data && data->guessMapping(surface_in) == DataBase::Element) {
        haveElementData = true;
    }
#endif

    //gather all domainsurfaces => at the moment all to all comm
    /* std::vector<Polygons::const_ptr> surfaceNeighbors; */
    /* b_mpi::all_gather(comm(), surface, surfaceNeighbors); */
    /* MPI_Barrier(comm()); */

    /* const auto num_elem = surface->getNumElements(); */
    /* const auto num_corners = surface->getNumCorners(); */
    /* const auto num_vertices = surface->getNumVertices(); */

    /* int proc{0}; */
    /* for (auto surface_neighbor: surfaceNeighbors) { */
    /*     /1*    _ _ _ _             _ _ _ _ */
    /*          /_/_/_/_/|         /_/_/_/_/               /| */
    /*         /_/_/_/_/ |        /_/_/_/_/               / | */
    /*        /_/_/_/_/  |       /_/_/_/_/               /  | */
    /*       /_/_/_/_/   |      /_/_/_/_/    _ _ _ _    /   | */
    /*       |_|_|_|_|   / =>               |_|_|_|_|   |   / */
    /*       |_|_|_|_|  /                   |_|_|_|_|   |  / */
    /*       |_|_|_|_| /                    |_|_|_|_|   | / */
    /*       |_|_|_|_|/                     |_|_|_|_|   |/   *1/ */

    /*     ++proc; */
    /*     std::vector<int> neighbors; */

    /*     //extract neighbors via vertices + send back to proc */
    /*     const auto n_num_elem = surface_neighbor->getNumElements(); */
    /*     const auto n_num_corners = surface_neighbor->getNumCorners(); */
    /*     const auto n_num_vertices = surface_neighbor->getNumVertices(); */

    /*     Polygons::ptr poly_build{ */
    /*         new Polygons(num_elem + n_num_elem, num_corners + n_num_corners, num_vertices + n_num_vertices)}; */

    /* const Index *el = &surface->el()[0]; */
    /* const Index *cl = &surface->cl()[0]; */
    /* Polygons::VertexOwnerList::const_ptr poly_vol = surface->getVertexOwnerList(); */

    /* Polygons::ptr m_grid_out(new Polygons(0, 0, 0)); */
    /* auto &pl = m_grid_out->el(); */
    /* auto &pcl = m_grid_out->cl(); */

    /* auto nf = m_grid_in->getNeighborFinder(); */
    /* for (Index i=0; i<num_elem; ++i) { */
    /*    const Index elStart = el[i], elEnd = el[i+1]; */
    /*    bool ghost = tl[i] & UnstructuredGrid::GHOST_BIT; */
    /*    if (!showgho && ghost) */
    /*        continue; */
    /*    Byte t = tl[i] & UnstructuredGrid::TYPE_MASK; */
    /*    if (t == UnstructuredGrid::POLYHEDRON) { */
    /*        if (showpol) { */
    /*            Index facestart = InvalidIndex; */
    /*            Index term = 0; */
    /*            for (Index j=elStart; j<elEnd; ++j) { */
    /*                if (facestart == InvalidIndex) { */
    /*                    facestart = j; */
    /*                    term = cl[j]; */
    /*                } else if (cl[j] == term) { */
    /*                    Index numVert = j - facestart; */
    /*                    if (numVert >= 3) { */
    /*                        auto face = &cl[facestart]; */
    /*                        Index neighbour = nf.getNeighborElement(i, face[0], face[1], face[2]); */
    /*                        if (neighbour == InvalidIndex) { */
    /*                            const Index *begin = &face[0], *end=&face[numVert]; */
    /*                            auto rbegin = std::reverse_iterator<const Index *>(end), rend = std::reverse_iterator<const Index *>(begin); */
    /*                            std::copy(rbegin, rend, std::back_inserter(pcl)); */
    /*                            if (haveElementData) */
    /*                                em.emplace_back(i); */
    /*                            pl.push_back(pcl.size()); */
    /*                        } */
    /*                    } */
    /*                    facestart = InvalidIndex; */
    /*                } */
    /*            } */
    /*        } */
    /*    } else { */
    /*        bool show = false; */
    /*        switch(t) { */
    /*        case UnstructuredGrid::PYRAMID: */
    /*            show = showpyr; */
    /*            break; */
    /*        case UnstructuredGrid::PRISM: */
    /*            show = showpri; */
    /*            break; */
    /*        case UnstructuredGrid::TETRAHEDRON: */
    /*            show = showtet; */
    /*            break; */
    /*        case UnstructuredGrid::HEXAHEDRON: */
    /*            show = showhex; */
    /*            break; */
    /*        case UnstructuredGrid::TRIANGLE: */
    /*            show = showtri; */
    /*            break; */
    /*        case UnstructuredGrid::QUAD: */
    /*            show = showqua; */
    /*            break; */
    /*        default: */
    /*            break; */
    /*        } */

    /*        if (show) { */
    /*          const auto numFaces = UnstructuredGrid::NumFaces[t]; */
    /*          const auto &faces = UnstructuredGrid::FaceVertices[t]; */
    /*          for (int f=0; f<numFaces; ++f) { */
    /*             const auto &face = faces[f]; */
    /*             Index neighbour = 0; */
    /*             if (UnstructuredGrid::Dimensionality[t] == 3) */
    /*                 neighbour = nf.getNeighborElement(i, cl[elStart + face[0]], cl[elStart + face[1]], cl[elStart + face[2]]); */
    /*             if (UnstructuredGrid::Dimensionality[t] == 2 || neighbour == InvalidIndex) { */
    /*                const auto facesize = UnstructuredGrid::FaceSizes[t][f]; */
    /*                for (unsigned j=0;j<facesize;++j) { */
    /*                   pcl.push_back(cl[elStart + face[j]]); */
    /*                } */
    /*                if (haveElementData) */
    /*                    em.emplace_back(i); */
    /*                pl.push_back(pcl.size()); */
    /*             } */
    /*          } */
    /*       } */
    /*    } */
    /* } */

    /* if (m_grid_out->getNumElements() == 0) { */
    /*    return Polygons::ptr(); */
    /* } */

    /* return m_grid_out; */
    /* } */
    /* Object::ptr surface; */
    /* broadcastObject(comm(), ugrid, ); */
    /* DataMapping vm; */
    /* DataMapping em; */
    /* if (ugrid) { */
    /*     auto poly = createSurface(ugrid, em, haveElementData); */
    /*     surface = poly; */
    /*     if (!poly) */
    /*         return true; */
    /*     renumberVertices(ugrid, poly, vm); */
    /* } else if (sgrid) { */
    /*     auto quad = createSurface(sgrid, em, haveElementData); */
    /*     surface = quad; */
    /*     if (!quad) */
    /*         return true; */
    /*     if (auto coords = Coords::as(grid_in)) { */
    /*         renumberVertices(coords, quad, vm); */
    /*     } else { */
    /*         createVertices(sgrid, quad, vm); */
    /*     } */
    /* } */

    /* surface->setMeta(grid_in->meta()); */
    /* surface->copyAttributes(grid_in); */

    /* if (!data) { */
    /*     updateMeta(surface); */
    /*     task->addObject("data_out", surface); */
    /*     return true; */
    /* } */

    /* if (!haveElementData && data->guessMapping(grid_in) != DataBase::Vertex) { */
    /*     sendError("data mapping not per vertex and not per element"); */
    /*     return true; */
    /* } */

    /* if (!haveElementData && vm.empty()) { */
    /*     DataBase::ptr dout = data->clone(); */
    /*     dout->setGrid(surface); */
    /*     updateMeta(dout); */
    /*     task->addObject("data_out", dout); */
    /*     return true; */
    /* } */

    /* DataBase::ptr data_obj_out; */
    /* const auto &dm = haveElementData ? em : vm; */
    /* if(auto data_in = Vec<Scalar, 3>::as(data)) { */
    /*     data_obj_out = remapData<Scalar,3>(data_in, dm); */
    /* } else if(auto data_in = Vec<Scalar,1>::as(data)) { */
    /*     data_obj_out = remapData<Scalar,1>(data_in, dm); */
    /* } else if(auto data_in = Vec<Index,3>::as(data)) { */
    /*     data_obj_out = remapData<Index,3>(data_in, dm); */
    /* } else if(auto data_in = Vec<Index,1>::as(data)) { */
    /*     data_obj_out = remapData<Index,1>(data_in, dm); */
    /* } else if(auto data_in = Vec<Byte,3>::as(data)) { */
    /*     data_obj_out = remapData<Byte,3>(data_in, dm); */
    /* } else if(auto data_in = Vec<Byte,1>::as(data)) { */
    /*     data_obj_out = remapData<Byte,1>(data_in, dm); */
    /* } else { */
    /*     std::cerr << "WARNING: No valid 1D or 3D element data on input Port" << std::endl; */
    /* } */

    /* if (data_obj_out) { */
    /*     data_obj_out->setGrid(surface); */
    /*     data_obj_out->setMeta(data->meta()); */
    /*     data_obj_out->copyAttributes(data); */
    /*     updateMeta(data_obj_out); */
    /*     task->addObject("data_out", data_obj_out); */
    /* } */

    return true;
}

MODULE_MAIN(GhostCellGenerator)
