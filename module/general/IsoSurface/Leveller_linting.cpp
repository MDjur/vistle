//
//This code is used for both IsoCut and IsoSurface!
//

#include <memory>
#include <sstream>
#include <iomanip>
#include <thrust/execution_policy.h>
#include <thrust/device_vector.h>
#include <thrust/transform.h>
#include <thrust/for_each.h>
#include <thrust/scan.h>
#include <thrust/sequence.h>
#include <thrust/copy.h>
#include <thrust/count.h>
#include <thrust/iterator/zip_iterator.h>
#include <thrust/iterator/constant_iterator.h>
#include <thrust/tuple.h>

#include "tables.h"
#include "HostData.h"
#include "DeviceData.h"
#include "Leveller.h"
#include "vistle/core/unstr.h"

using namespace vistle;

constexpr int MaxNumData = 6;

template<class Data>
struct ComputeOutput {
    Data &m_data;

    ComputeOutput(Data &data): m_data(data) { m_data.prepareForOutput(); }

    __host__ __device__ void operator()(Index ValidCellIndex)
    {
        const Index CellNr = m_data.SelectedCellVector()[ValidCellIndex];
        m_data.copyCellInToOutAtIdx(ValidCellIndex, CellNr);

        if (m_data.IsUnstructured()) {
            const Index Cellbegin = m_data.El()[CellNr];
            const Index Cellend = m_data.El()[CellNr + 1];
            const auto &unstrType = m_data.Tl()[CellNr] & ~UnstructuredGrid::CONVEX_BIT;

            switch (unstrType) {
            case UnstructuredGrid::HEXAHEDRON:
                m_data.computeGeoOutput(CellNr, unstrType, ValidCellIndex, 0, hexaTriTable, hexaEdgeTable);
                break;
            case UnstructuredGrid::TETRAHEDRON:
                m_data.computeGeoOutput(CellNr, unstrType, ValidCellIndex, 0, tetraTriTable, tetraEdgeTable);
                break;
            case UnstructuredGrid::PYRAMID:
                m_data.computeGeoOutput(CellNr, unstrType, ValidCellIndex, 0, pyrTriTable, pyrEdgeTable);
                break;
            case UnstructuredGrid::PRISM:
                m_data.computeGeoOutput(CellNr, unstrType, ValidCellIndex, 0, prismTriTable, prismEdgeTable);
                break;
            /* case UnstructuredGrid::VPOLYHEDRON: */
            /*     m_data.interpolUnstrPolyhedron(Cellbegin, Cellend, ValidCellIndex); */
            /*     break; */
            case UnstructuredGrid::CPOLYHEDRON:
                m_data.interpolUnstrPolyhedron(Cellbegin, Cellend, ValidCellIndex);
                break;
            }
        } else if (m_data.m_isTri) {
            const Index Cellbegin = CellNr * 3;
            m_data.computeGeoOutput(Cellbegin, 3, ValidCellIndex, 0, triLineTable, triEdgeTable);
        } else if (m_data.m_isQuad) {
            const Index Cellbegin = CellNr * 4;
            m_data.computeGeoOutput(Cellbegin, 4, ValidCellIndex, 0, quadLineTable, quadEdgeTable);
        } else if (m_data.m_isPoly) {
        } else {
            m_data.computeStructOutput(CellNr);
        }
    }
};

template<class Data>
struct SelectCells {
    typedef float argument_type;
    typedef float result_type;
    Data &m_data;
    SelectCells(Data &data): m_data(data) {}

    // for unstructured grids
    __host__ __device__ int operator()(const thrust::tuple<Index, Index, Byte> iCell) const
    {
        int havelower = 0;
        int havehigher = 0;
        Index Cell = iCell.get<0>();
        Index nextCell = iCell.get<1>();
        Byte cellType = iCell.get<2>();
        if (cellType & UnstructuredGrid::GHOST_BIT)
            return 0;
        if ((cellType & UnstructuredGrid::TYPE_MASK) == UnstructuredGrid::VPOLYHEDRON) {
            for (Index i = Cell; i < nextCell; i++) {
                Index nv = m_data.Cl()[i];
                for (Index k = 0; k < nv; k++) {
                    ++i;
                    float val = m_data.IsoFunc()(m_data.Cl()[i]);
                    if (val > m_data.Isovalue()) {
                        havelower = 1;
                        if (havehigher)
                            return 1;
                    } else {
                        havehigher = 1;
                        if (havelower)
                            return 1;
                    }
                }
            }
        } else {
            // also for CPOLYHEDRON
            for (Index i = Cell; i < nextCell; i++) {
                float val = m_data.IsoFunc()(m_data.Cl()[i]);
                if (val > m_data.Isovalue()) {
                    havelower = 1;
                    if (havehigher)
                        return 1;
                } else {
                    havehigher = 1;
                    if (havelower)
                        return 1;
                }
            }
        }
        return 0;
    }

    // for all types of structured grids
    __host__ __device__ int operator()(const Index Cell) const
    {
        auto cc = vistle::StructuredGridBase::cellCoordinates(Cell, m_data.Nvert());
        for (int c = 0; c < 3; ++c) {
            if (cc[c] < m_data.Nghost()[c][0])
                return 0;
            if (cc[c] + m_data.Nghost()[c][1] + 1 >= m_data.Nvert()[c])
                return 0;
        }

        auto verts = vistle::StructuredGridBase::cellVertices(Cell, m_data.Nvert());
        int havelower = 0;
        int havehigher = 0;
        for (int i = 0; i < 8; ++i) {
            float val = m_data.IsoFunc()(verts[i]);
            if (val > m_data.Isovalue()) {
                havelower = 1;
                if (havehigher)
                    return 1;
            } else {
                havehigher = 1;
                if (havelower)
                    return 1;
            }
        }
        return 0;
    }
};

template<class Data>
struct SelectCells2D {
    typedef float argument_type;
    typedef float result_type;
    Data &m_data;
    SelectCells2D(Data &data): m_data(data) {}

    // for polygons
    __host__ __device__ int operator()(const thrust::tuple<Index, Index, Byte> iCell) const
    {
        int havelower = 0;
        int havehigher = 0;
        Index Cell = iCell.get<0>();
        Index nextCell = iCell.get<1>();
        for (Index i = Cell; i < nextCell; i++) {
            Index nv = m_data.Cl()[i];
            for (Index k = 0; k < nv; k++) {
                ++i;
                float val = m_data.IsoFunc()(m_data.Cl()[i]);
                if (val > m_data.Isovalue()) {
                    havelower = 1;
                    if (havehigher)
                        return 1;
                } else {
                    havehigher = 1;
                    if (havelower)
                        return 1;
                }
            }
        }
        return 0;
    }

    // for triangles and quads
    __host__ __device__ int operator()(const Index Cell) const
    {
        int havelower = 0;
        int havehigher = 0;
        Index begin = Cell * m_data.NumVertPerCell(), end = begin + m_data.NumVertPerCell();
        if (m_data.Cl()) {
            for (Index i = begin; i < end; ++i) {
                float val = m_data.IsoFunc()(m_data.Cl()[i]);
                if (val > m_data.Isovalue()) {
                    havelower = 1;
                    if (havehigher)
                        return 1;
                } else {
                    havehigher = 1;
                    if (havelower)
                        return 1;
                }
            }
        } else {
            for (Index i = begin; i < end; ++i) {
                float val = m_data.IsoFunc()(i);
                if (val > m_data.Isovalue()) {
                    havelower = 1;
                    if (havehigher)
                        return 1;
                } else {
                    havehigher = 1;
                    if (havelower)
                        return 1;
                }
            }
        }
        return 0;
    }
};


template<class Data>
struct ComputeOutputSizes {
    ComputeOutputSizes(Data &data): m_data(data) {}

    Data &m_data;

    __host__ __device__ thrust::tuple<Index, Index> operator()(Index CellNr)
    {
        int tableIndex = 0;
        unsigned numVerts = 0;
        if (m_data.IsUnstructured()) {
            const auto &cl = m_data.Cl();

            Index begin = m_data.El()[CellNr], end = m_data.El()[CellNr + 1];
            Index nvert = end - begin;
            Byte CellType = m_data.Tl()[CellNr] & ~UnstructuredGrid::CONVEX_BIT;
            if (CellType != UnstructuredGrid::VPOLYHEDRON && CellType != UnstructuredGrid::CPOLYHEDRON) {
                for (Index idx = 0; idx < nvert; idx++) {
                    tableIndex += (((int)(m_data.IsoFunc()(m_data.Cl()[begin + idx]) > m_data.Isovalue())) << idx);
                }
            }
            switch (CellType) {
            case UnstructuredGrid::HEXAHEDRON:
                numVerts = hexaNumVertsTable[tableIndex];
                break;

            case UnstructuredGrid::TETRAHEDRON:
                numVerts = tetraNumVertsTable[tableIndex];
                break;

            case UnstructuredGrid::PYRAMID:
                numVerts = pyrNumVertsTable[tableIndex];
                break;

            case UnstructuredGrid::PRISM:
                numVerts = prismNumVertsTable[tableIndex];
                break;

            case UnstructuredGrid::VPOLYHEDRON: {
                Index vertcounter = 0;
                for (Index i = begin; i < end; i += cl[i] + 1) {
                    const Index N = cl[i];
                    Index prev = cl[i + N];
                    for (Index k = i + 1; k < i + N + 1; ++k) {
                        Index v = cl[k];

                        if (m_data.IsoFunc()(prev) <= m_data.Isovalue() && m_data.IsoFunc()(v) > m_data.Isovalue()) {
                            ++vertcounter;
                        } else if (m_data.IsoFunc()(prev) > m_data.Isovalue() &&
                                   m_data.IsoFunc()(v) <= m_data.Isovalue()) {
                            ++vertcounter;
                        }

                        prev = v;
                    }
                }
                numVerts = vertcounter + vertcounter / 2;
                break;
            }

            case UnstructuredGrid::CPOLYHEDRON: {
                Index vertcounter = 0;
                Index facestart = InvalidIndex;
                Index term = 0;
                for (Index i = begin; i < end; ++i) {
                    if (facestart == InvalidIndex) {
                        facestart = i;
                        term = cl[i];
                    } else if (term == cl[i]) {
                        const Index N = i - facestart;
                        Index prev = cl[facestart + N - 1];
                        for (Index k = facestart; k < facestart + N; ++k) {
                            Index v = cl[k];

                            if (m_data.IsoFunc()(prev) <= m_data.Isovalue() &&
                                m_data.IsoFunc()(v) <= m_data.Isovalue()) {
                                ++vertcounter;
                            } else if (m_data.IsoFunc()(prev) > m_data.Isovalue() &&
                                       m_data.IsoFunc()(v) <= m_data.Isovalue()) {
                                ++vertcounter;
                            }

                            prev = v;
                        }
                        facestart = InvalidIndex;
                    }
                }
                numVerts = vertcounter + vertcounter / 2;
                break;
            }
            }
        } else if (m_data.IsTri() || m_data.IsQuad()) {
            const auto &cl = m_data.Cl();
            const Index begin = CellNr * m_data.NumVertPerCell(), end = begin + m_data.NumVertPerCell();
            if (cl) {
                int idx = 0;
                for (Index i = begin; i < end; ++i) {
                    tableIndex += (((int)(m_data.IsoFunc()(cl[i]) > m_data.Isovalue())) << idx);
                    ++idx;
                }
            } else {
                int idx = 0;
                for (Index i = begin; i < end; ++i) {
                    tableIndex += (((int)(m_data.IsoFunc()(i) > m_data.Isovalue())) << idx);
                    ++idx;
                }
            }
            numVerts = m_data.IsQuad() ? quadNumVertsTable[tableIndex] : triNumVertsTable[tableIndex];
        } else if (m_data.IsPoly()) {
            const auto &cl = m_data.Cl();
            Index begin = m_data.El()[CellNr], end = m_data.El()[CellNr + 1];
            Index vertcounter = 0;
            Index prev = cl[end];
            for (Index i = begin; i < end; ++i) {
                Index v = cl[i];
                if (m_data.IsoFunc()(prev) <= m_data.Isovalue() && m_data.IsoFunc()(v) > m_data.Isovalue()) {
                    ++vertcounter;
                } else if (m_data.IsoFunc()(prev) > m_data.Isovalue() && m_data.IsoFunc()(v) <= m_data.Isovalue()) {
                    ++vertcounter;
                }

                prev = v;
            }
            numVerts = vertcounter + vertcounter / 2;
        } else {
            auto verts = vistle::StructuredGridBase::cellVertices(CellNr, m_data.Nvert());
            assert(verts.size() <= 8);
            for (unsigned idx = 0; idx < verts.size(); ++idx) {
                tableIndex += (((int)(m_data.IsoFunc()(verts[idx]) > m_data.Isovalue())) << idx);
            }
            numVerts = hexaNumVertsTable[tableIndex];
        }
        return thrust::make_tuple<Index, Index>(tableIndex, numVerts);
    }
};

Leveller::Leveller(const IsoController &isocontrol, Object::const_ptr grid, const Scalar isovalue, Index processortype)
: m_isocontrol(isocontrol)
, m_grid(grid)
, m_uni(UniformGrid::as(grid))
, m_lg(LayerGrid::as(grid))
, m_rect(RectilinearGrid::as(grid))
, m_str(StructuredGrid::as(grid))
, m_unstr(UnstructuredGrid::as(grid))
, m_strbase(StructuredGridBase::as(grid))
, m_poly(Polygons::as(grid))
, m_quad(Quads::as(grid))
, m_tri(Triangles::as(grid))
, m_coord(Coords::as(grid))
, m_isoValue(isovalue)
, m_processortype(processortype)
, gmin(std::numeric_limits<Scalar>::max())
, gmax(-std::numeric_limits<Scalar>::max())
, m_objectTransform(grid->getTransform())
{
    if (m_strbase || m_unstr) {
        m_triangles = Triangles::ptr(new Triangles(Object::Initialized));
        m_triangles->setMeta(grid->meta());
    } else if (m_poly || m_tri || m_quad) {
        m_lines = Lines::ptr(new Lines(Object::Initialized));
        m_lines->setMeta(grid->meta());
    }
}

template<class Data, class pol>
Index Leveller::calculateSurface(Data &data)
{
    Index nelem = 0;
    if (m_strbase) {
        nelem = m_strbase->getNumElements();
    } else if (m_unstr) {
        nelem = m_unstr->getNumElements();
    } else if (m_tri) {
        nelem = m_tri->getNumElements();
    } else if (m_quad) {
        nelem = m_quad->getNumElements();
    } else if (m_poly) {
        nelem = m_poly->getNumElements();
    }
    thrust::counting_iterator<Index> first(0), last = first + nelem;
    data.SelectedCellVector().resize(nelem);

    typedef thrust::tuple<typename Data::IndexIterator, typename Data::IndexIterator, typename Data::TypeIterator>
        Iteratortuple;
    typedef thrust::zip_iterator<Iteratortuple> ZipIterator;

    typename Data::VectorIndexIterator end;
    if (m_strbase) {
        end = thrust::copy_if(pol(), first, last, thrust::counting_iterator<Index>(0),
                              data.SelectedCellVector().begin(), SelectCells<Data>(data));
    } else if (m_unstr) {
        ZipIterator ElTupleVec(thrust::make_tuple(&data.El()[0], &data.El()[1], &data.Tl()[0]));
        end =
            thrust::copy_if(pol(), first, last, ElTupleVec, data.SelectedCellVector().begin(), SelectCells<Data>(data));
    } else if (m_poly) {
        ZipIterator ElTupleVec(thrust::make_tuple(&data.El()[0], &data.El()[1], &data.Tl()[0]));
        end = thrust::copy_if(pol(), first, last, ElTupleVec, data.SelectedCellVector().begin(),
                              SelectCells2D<Data>(data));
    } else if (m_tri || m_quad) {
        end = thrust::copy_if(pol(), first, last, thrust::counting_iterator<Index>(0),
                              data.SelectedCellVector().begin(), SelectCells2D<Data>(data));
    }

    size_t numSelectedCells = end - data.SelectedCellVector().begin();
    data.CaseNums().resize(numSelectedCells);
    data.NumVertices().resize(numSelectedCells);
    data.LocationList().resize(numSelectedCells);
    thrust::transform(
        pol(), data.SelectedCellVector().begin(), end,
        thrust::make_zip_iterator(thrust::make_tuple(data.CaseNums().begin(), data.NumVertices().begin())),
        ComputeOutputSizes<Data>(data));
    thrust::exclusive_scan(pol(), data.NumVertices().begin(), data.NumVertices().end(), data.LocationList().begin());
    Index totalNumVertices = 0;
    if (!data.NumVertices().empty())
        totalNumVertices += data.NumVertices().back();
    if (!data.LocationList().empty())
        totalNumVertices += data.LocationList().back();
    for (int i = (m_computeNormals || !m_strbase ? 0 : 3); i < data.NumInVertData(); i++) {
        data.VertData().outVecData[i]->resize(totalNumVertices);
    }
    for (int i = 0; i < data.NumInVertDataI(); i++) {
        data.VertDataIdx().outVecData[i]->resize(totalNumVertices);
    }
    for (int i = 0; i < data.NumInVertDataB(); i++) {
        data.VertDataByte().outVecData[i]->resize(totalNumVertices);
    }
    for (int i = 0; i < data.NumInCellData(); ++i) {
        data.CellData().outVecData[i]->resize(totalNumVertices / 3);
    }
    for (int i = 0; i < data.NumInCellDataI(); ++i) {
        data.CellDataIdx().outVecData[i]->resize(totalNumVertices / 3);
    }
    for (int i = 0; i < data.NumInCellDataB(); ++i) {
        data.CellDataByte().outVecData[i]->resize(totalNumVertices / 3);
    }
    thrust::counting_iterator<Index> start(0), finish(numSelectedCells);
    thrust::for_each(pol(), start, finish, ComputeOutput<Data>(data));

    return totalNumVertices;
}

bool Leveller::process()
{
#ifndef CUTTINGSURFACE
    Vec<Scalar>::const_ptr dataobj = Vec<Scalar>::as(m_data);
    if (!dataobj)
        return false;
    auto bounds = dataobj->getMinMax();
    if (bounds.first[0] <= bounds.second[0]) {
        if (m_isoValue < bounds.first[0] || m_isoValue > bounds.second[0])
            return true;
    }
#else
#endif

    switch (m_processortype) {
    case Host: {
        Index dims[3] = {0, 0, 0};
        if (m_strbase) {
            dims[0] = m_strbase->getNumDivisions(0);
            dims[1] = m_strbase->getNumDivisions(1);
            dims[2] = m_strbase->getNumDivisions(2);
        }
        std::vector<Scalar> unicoords[3];
        const Scalar *coords[3]{nullptr, nullptr, nullptr};
        if (m_uni) {
            for (int i = 0; i < 3; ++i) {
                unicoords[i].resize(dims[i]);
                coords[i] = unicoords[i].data();
                Scalar dist = 0;
                if (dims[i] > 1)
                    dist = (m_uni->max()[i] - m_uni->min()[i]) / (dims[i] - 1);
                Scalar val = m_uni->min()[i];
                for (Index j = 0; j < dims[i]; ++j) {
                    unicoords[i][j] = val;
                    val += dist;
                }
            }
        } else if (m_lg) {
            for (int i = 0; i < 2; ++i) {
                unicoords[i].resize(dims[i]);
                coords[i] = unicoords[i].data();
                Scalar dist = 0;
                if (dims[i] > 1)
                    dist = (m_lg->max()[i] - m_lg->min()[i]) / (dims[i] - 1);
                Scalar val = m_lg->min()[i];
                for (Index j = 0; j < dims[i]; ++j) {
                    unicoords[i][j] = val;
                    val += dist;
                }
            }
            coords[2] = &m_lg->z()[0];
        } else if (m_rect) {
            for (int i = 0; i < 3; ++i)
                coords[i] = &m_rect->coords(i)[0];
        } else if (m_str) {
            for (int i = 0; i < 3; ++i)
                coords[i] = &m_str->x(i)[0];
        }
#ifdef CUTTINGSURFACE
        IsoDataFunctor isofunc =
            m_coord ? m_isocontrol.newFunc(m_grid->getTransform(), &m_coord->x()[0], &m_coord->y()[0], &m_coord->z()[0])
                    : m_isocontrol.newFunc(m_grid->getTransform(), dims, coords[0], coords[1], coords[2]);
#else
        IsoDataFunctor isofunc = m_isocontrol.newFunc(m_grid->getTransform(), &dataobj->x()[0]);
#endif

        HostData hostData =
            m_unstr     ? HostData(m_isoValue, isofunc, m_unstr->el(), m_unstr->tl(), m_unstr->cl(), m_unstr->x(),
                                   m_unstr->y(), m_unstr->z())
            : m_strbase ? HostData(m_isoValue, isofunc, dims[0], dims[1], dims[2], coords[0], coords[1], coords[2])
            : m_poly ? HostData(m_isoValue, isofunc, m_poly->el(), m_poly->cl(), m_poly->x(), m_poly->y(), m_poly->z())
            : m_quad ? HostData(m_isoValue, isofunc, 4, m_quad->cl(), m_quad->x(), m_quad->y(), m_quad->z())
                     : HostData(m_isoValue, isofunc, 3, m_tri->cl(), m_tri->x(), m_tri->y(), m_tri->z());
        hostData.setHaveCoords(m_coord ? true : false);
        if (m_strbase) {
            Index ghost[3][2];
            for (int c = 0; c < 3; ++c) {
                ghost[c][0] = m_strbase->getNumGhostLayers(c, StructuredGridBase::Bottom);
                ghost[c][1] = m_strbase->getNumGhostLayers(c, StructuredGridBase::Top);
            }
            hostData.setGhostLayers(ghost);
        }
        hostData.setComputeNormals(m_computeNormals);

        for (size_t i = 0; i < m_vertexdata.size(); ++i) {
            if (Vec<Scalar, 1>::const_ptr Scal = Vec<Scalar, 1>::as(m_vertexdata[i])) {
                hostData.addMappedData(Scal->x());
            }
            if (Vec<Scalar, 3>::const_ptr Vect = Vec<Scalar, 3>::as(m_vertexdata[i])) {
                hostData.addMappedData(Vect->x());
                hostData.addMappedData(Vect->y());
                hostData.addMappedData(Vect->z());
            }
            if (Vec<Index, 1>::const_ptr Idx = Vec<Index, 1>::as(m_vertexdata[i])) {
                hostData.addMappedData(Idx->x());
            }
            if (auto byte = Vec<Byte>::as(m_vertexdata[i])) {
                hostData.addMappedData(byte->x());
            }
        }
        for (size_t i = 0; i < m_celldata.size(); ++i) {
            if (Vec<Scalar, 1>::const_ptr Scal = Vec<Scalar, 1>::as(m_celldata[i])) {
                hostData.addCellData(Scal->x());
            }
            if (Vec<Scalar, 3>::const_ptr Vect = Vec<Scalar, 3>::as(m_celldata[i])) {
                hostData.addCellData(Vect->x());
                hostData.addCellData(Vect->y());
                hostData.addCellData(Vect->z());
            }
            if (Vec<Index, 1>::const_ptr Idx = Vec<Index, 1>::as(m_celldata[i])) {
                hostData.addCellData(Idx->x());
            }
            if (auto byte = Vec<Byte>::as(m_celldata[i])) {
                hostData.addCellData(byte->x());
            }
        }

        Index totalNumVertices = calculateSurface<HostData, thrust::detail::host_t>(hostData);

        {
            size_t idx = 0;
            if (m_strbase) {
                if (m_computeNormals) {
                    m_normals.reset(new Normals(Object::Initialized));
                    m_normals->d()->x[0] = hostData.VertData().outVecData[idx++];
                    m_normals->d()->x[1] = hostData.VertData().outVecData[idx++];
                    m_normals->d()->x[2] = hostData.VertData().outVecData[idx++];
                    m_normals->setMapping(DataBase::Vertex);
                } else {
                    idx = 3;
                }
            }

            if (m_triangles) {
                m_triangles->d()->x[0] = hostData.VertData().outVecData[idx++];
                m_triangles->d()->x[1] = hostData.VertData().outVecData[idx++];
                m_triangles->d()->x[2] = hostData.VertData().outVecData[idx++];
            } else if (m_lines) {
                m_lines->d()->x[0] = hostData.VertData().outVecData[idx++];
                m_lines->d()->x[1] = hostData.VertData().outVecData[idx++];
                m_lines->d()->x[2] = hostData.VertData().outVecData[idx++];
                std::cerr << "lines with " << totalNumVertices << " vertices" << std::endl;
                auto &cl = m_lines->cl();
                auto &el = m_lines->el();
                for (Index i = 0; i < totalNumVertices; ++i) {
                    cl.emplace_back(i);
                    if (cl.size() % 2 == 0)
                        el.emplace_back(cl.size());
                }
            }

            size_t idxI = 0, idxB = 0;
            for (size_t i = 0; i < m_vertexdata.size(); ++i) {
                if (Vec<Scalar>::as(m_vertexdata[i])) {
                    Vec<Scalar, 1>::ptr out = Vec<Scalar, 1>::ptr(new Vec<Scalar, 1>(Object::Initialized));
                    out->d()->x[0] = hostData.VertData().outVecData[idx++];
                    out->setMeta(m_vertexdata[i]->meta());
                    out->setMapping(DataBase::Vertex);
                    m_outvertData.push_back(out);
                }
                if (Vec<Scalar, 3>::as(m_vertexdata[i])) {
                    Vec<Scalar, 3>::ptr out = Vec<Scalar, 3>::ptr(new Vec<Scalar, 3>(Object::Initialized));
                    out->d()->x[0] = hostData.VertData().outVecData[idx++];
                    out->d()->x[1] = hostData.VertData().outVecData[idx++];
                    out->d()->x[2] = hostData.VertData().outVecData[idx++];
                    out->setMeta(m_vertexdata[i]->meta());
                    out->setMapping(DataBase::Vertex);
                    m_outvertData.push_back(out);
                }
                if (Vec<Index>::as(m_vertexdata[i])) {
                    Vec<Index>::ptr out = Vec<Index>::ptr(new Vec<Index>(Object::Initialized));
                    out->d()->x[0] = hostData.VertDataIdx().outVecData[idxI++];
                    out->setMeta(m_vertexdata[i]->meta());
                    out->setMapping(DataBase::Vertex);
                    m_outvertData.push_back(out);
                }
                if (Vec<Byte>::as(m_vertexdata[i])) {
                    Vec<Byte>::ptr out = Vec<Byte>::ptr(new Vec<Byte>(Object::Initialized));
                    out->d()->x[0] = hostData.VertDataByte().outVecData[idxB++];
                    out->setMeta(m_vertexdata[i]->meta());
                    out->setMapping(DataBase::Vertex);
                    m_outvertData.push_back(out);
                }
            }
        }
        {
            size_t idx = 0;
            size_t idxI = 0;
            size_t idxB = 0;
            for (size_t i = 0; i < m_celldata.size(); ++i) {
                if (Vec<Scalar>::as(m_celldata[i])) {
                    Vec<Scalar, 1>::ptr out = Vec<Scalar, 1>::ptr(new Vec<Scalar, 1>(Object::Initialized));
                    out->d()->x[0] = hostData.CellData().outVecData[idx++];
                    out->setMeta(m_celldata[i]->meta());
                    out->setMapping(DataBase::Element);
                    m_outcellData.push_back(out);
                }
                if (Vec<Scalar, 3>::as(m_celldata[i])) {
                    Vec<Scalar, 3>::ptr out = Vec<Scalar, 3>::ptr(new Vec<Scalar, 3>(Object::Initialized));
                    out->d()->x[0] = hostData.CellData().outVecData[idx++];
                    out->d()->x[1] = hostData.CellData().outVecData[idx++];
                    out->d()->x[2] = hostData.CellData().outVecData[idx++];
                    out->setMeta(m_celldata[i]->meta());
                    out->setMapping(DataBase::Element);
                    m_outcellData.push_back(out);
                }
                if (Vec<Index>::as(m_celldata[i])) {
                    Vec<Index>::ptr out = Vec<Index>::ptr(new Vec<Index>(Object::Initialized));
                    out->d()->x[0] = hostData.CellDataIdx().outVecData[idxI++];
                    out->setMeta(m_celldata[i]->meta());
                    out->setMapping(DataBase::Element);
                    m_outcellData.push_back(out);
                }
                if (Vec<Byte>::as(m_celldata[i])) {
                    Vec<Byte>::ptr out = Vec<Byte>::ptr(new Vec<Byte>(Object::Initialized));
                    out->d()->x[0] = hostData.CellDataByte().outVecData[idxB++];
                    out->setMeta(m_celldata[i]->meta());
                    out->setMapping(DataBase::Element);
                    m_outcellData.push_back(out);
                }
            }
        }
        break;
    }

    case Device: {
        std::cerr << "untested Device code path" << std::endl;
        assert("don't use the Device code path" == 0);

        DeviceData DD(m_isoValue,
#ifndef CUTTINGSURFACE
                      m_isocontrol.newFunc(m_grid->getTransform(), &dataobj->x()[0]),
#else
                      m_isocontrol.newFunc(m_grid->getTransform(), &m_coord->x()[0], &m_coord->y()[0],
                                           &m_coord->z()[0]),
#endif
                      m_unstr->getNumElements(), m_unstr->el(), m_unstr->tl(), m_unstr->getNumCorners(), m_unstr->cl(),
                      m_unstr->getSize(), m_unstr->x(), m_unstr->y(), m_unstr->z());

#if 0
         Index totalNumVertices = calculateSurface<DeviceData, thrust::detail::device_t>(DD);
#else
        Index totalNumVertices = 0;
#endif

        if (m_triangles) {
            m_triangles->x().resize(totalNumVertices);
            Scalar *out_x = m_triangles->x().data();
            thrust::copy(DD.m_outVertData[0]->begin(), DD.m_outVertData[0]->end(), out_x);

            m_triangles->y().resize(totalNumVertices);
            Scalar *out_y = m_triangles->y().data();
            thrust::copy(DD.m_outVertData[1]->begin(), DD.m_outVertData[1]->end(), out_y);

            m_triangles->z().resize(totalNumVertices);
            Scalar *out_z = m_triangles->z().data();
            thrust::copy(DD.m_outVertData[2]->begin(), DD.m_outVertData[2]->end(), out_z);
        } else if (m_lines) {
            m_lines->x().resize(totalNumVertices);
            Scalar *out_x = m_lines->x().data();
            thrust::copy(DD.m_outVertData[0]->begin(), DD.m_outVertData[0]->end(), out_x);

            m_lines->y().resize(totalNumVertices);
            Scalar *out_y = m_lines->y().data();
            thrust::copy(DD.m_outVertData[1]->begin(), DD.m_outVertData[1]->end(), out_y);

            m_lines->z().resize(totalNumVertices);
            Scalar *out_z = m_lines->z().data();
            thrust::copy(DD.m_outVertData[2]->begin(), DD.m_outVertData[2]->end(), out_z);
        }

        if (m_vertexdata.size()) {
            if (Vec<Scalar>::as(m_vertexdata[0])) {
                Vec<Scalar>::ptr out = Vec<Scalar>::ptr(new Vec<Scalar>(totalNumVertices));
                thrust::copy(DD.m_outVertData[3]->begin(), DD.m_outVertData[3]->end(), out->x().data());
                out->setMeta(m_vertexdata[0]->meta());
                m_outvertData.push_back(out);
            }
            if (Vec<Scalar, 3>::as(m_vertexdata[0])) {
                Vec<Scalar, 3>::ptr out = Vec<Scalar, 3>::ptr(new Vec<Scalar, 3>(totalNumVertices));
                thrust::copy(DD.m_outVertData[3]->begin(), DD.m_outVertData[3]->end(), out->x().data());
                thrust::copy(DD.m_outVertData[4]->begin(), DD.m_outVertData[4]->end(), out->y().data());
                thrust::copy(DD.m_outVertData[5]->begin(), DD.m_outVertData[5]->end(), out->z().data());
                out->setMeta(m_vertexdata[0]->meta());
                m_outvertData.push_back(out);
            }
        }
        break;
    }
    }

    return true;
}

#ifndef CUTTINGSURFACE
void Leveller::setIsoData(Vec<Scalar>::const_ptr obj)
{
    m_data = obj;
}
#endif

void Leveller::setComputeNormals(bool value)
{
    m_computeNormals = value;
}

void Leveller::addMappedData(DataBase::const_ptr mapobj)
{
    if (mapobj->mapping() == DataBase::Element)
        m_celldata.push_back(mapobj);
    else
        m_vertexdata.push_back(mapobj);
}

Coords::ptr Leveller::result()
{
    if (m_triangles)
        return m_triangles;
    return m_lines;
}

Normals::ptr Leveller::normresult()
{
    return m_normals;
}

DataBase::ptr Leveller::mapresult() const
{
    if (m_outvertData.size())
        return m_outvertData[0];
    else if (m_outcellData.size())
        return m_outcellData[0];
    else
        return DataBase::ptr();
}

DataBase::ptr Leveller::cellresult() const
{
    if (m_outcellData.size())
        return m_outcellData[0];
    else
        return DataBase::ptr();
}

std::pair<Scalar, Scalar> Leveller::range()
{
    return std::make_pair(gmin, gmax);
}
