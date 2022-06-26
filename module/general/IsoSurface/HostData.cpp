#include "HostData.h"
#include "tables.h"

using namespace vistle;

constexpr Scalar EPSILON(1e-10);

HostData::HostData(Scalar isoValue, IsoDataFunctor isoFunc, const Scalar *x, const Scalar *y, const Scalar *z)
: IsoDataI(isoValue, isoFunc)
, m_vertDataContainer(m_vertData, m_vertDataIdx, m_vertDataByte)
, m_cellDataContainer(m_vertData, m_vertDataIdx, m_vertDataByte)
{
    addMappedData(x);
    addMappedData(y);
    addMappedData(z);
}

//unstructured
HostData::HostData(Scalar isoValue, IsoDataFunctor isoFunc, const Index *el, const Byte *tl, const Index *cl,
                   const Scalar *x, const Scalar *y, const Scalar *z)
: HostData(isoValue, isoFunc, x, y, z)
{
    m_el = el;
    m_cl = cl;
    m_tl = tl;
    m_nvert[0] = 0;
    m_nvert[1] = 0;
    m_nvert[2] = 0;
    m_isUnstructured = el;
    m_haveCoords = true;
}

//structured
HostData::HostData(Scalar isoValue, IsoDataFunctor isoFunc, Index nx, Index ny, Index nz, const Scalar *x,
                   const Scalar *y, const Scalar *z)
: HostData(isoValue, isoFunc, nullptr, nullptr, nullptr, x, y, z)
{
    m_nvert[0] = nx;
    m_nvert[1] = ny;
    m_nvert[2] = nz;
    m_isUnstructured = false;
    m_haveCoords = false;
}

//poly
HostData::HostData(Scalar isoValue, IsoDataFunctor isoFunc, const Index *el, const Index *cl, const Scalar *x,
                   const Scalar *y, const Scalar *z)
: HostData(isoValue, isoFunc, x, y, z)
{
    m_el = el;
    m_cl = cl;
    m_tl = nullptr;
    m_isPoly = true;
}

//triangles/quads
HostData::HostData(Scalar isoValue, IsoDataFunctor isoFunc, int vertPerCell, const Index *cl, const Scalar *x,
                   const Scalar *y, const Scalar *z)
: HostData(isoValue, isoFunc, x, y, z)
{
    m_el = nullptr;
    m_cl = cl;
    m_tl = nullptr;
    m_isTri = vertPerCell == 3;
    m_isQuad = vertPerCell == 4;
    m_numVertPerCell = vertPerCell;
}

void HostData::prepareForOutput()
{
    m_vertDataContainer.prepareContForOutput(m_numInVertData, m_numInVertDataI, m_numInVertDataB);
    m_vertDataContainer.prepareContForOutput(m_numInCellData, m_numInCellDataI, m_numInCellDataB);
}

void HostData::copyCellInToOutAtIdx(Index ValidCellIndex, const Index CellNr)
{
    for (Index idx = 0; idx < m_numVertices[ValidCellIndex] / 3; ++idx) {
        Index outCellIdx = m_LocationList[ValidCellIndex] / 3 + idx;
        m_cellDataContainer.copyContInToOutPerElem(m_numInCellData, m_numInCellDataI, m_numInCellDataB, outCellIdx,
                                                   CellNr);
    }
}

template<size_t xT, size_t yT, size_t xE, size_t yE>
void HostData::computeGeoOutput(const Index Cellbegin, const Index fieldSize, const Index ValidCellIndex,
                                const Index nc, const int (&triTable)[xT][yT], const int (&edgeTable)[xE][yE])
{
    Scalar field[fieldSize];
    prepareField(field, fieldSize, Cellbegin);
    for (Index idx = 0; idx < m_numVertices[ValidCellIndex]; ++idx) {
        const unsigned int edge = triTable[m_caseNums[ValidCellIndex]][idx];
        const unsigned int edgeIdx1 = edgeTable[0][edge];
        const unsigned int edgeIdx2 = edgeTable[1][edge];
        const Scalar t = interpolation_weight<Scalar>(field[edgeIdx1], field[edgeIdx2], m_isovalue);
        Index outVertIdx = m_LocationList[ValidCellIndex] + idx;
        Index interpolEdgeIdx1 = Cellbegin + edgeIdx1;
        Index interpolEdgeIdx2 = Cellbegin + edgeIdx2;
        if (m_cl) { // else flat interpolation based on Cellbegin
            const auto &cl = &m_cl[Cellbegin];
            interpolEdgeIdx1 = cl[edgeIdx1];
            interpolEdgeIdx2 = cl[edgeIdx2];
        }
        interpolate(outVertIdx, field, triTable, edgeTable, interpolEdgeIdx1, interpolEdgeIdx2, t, nc);
    }
}

void HostData::computeStructuredOutput(const Index CellNr, const Index ValidCellIndex)
{
    auto cl = vistle::StructuredGridBase::cellVertices(CellNr, m_nvert);
    assert(cl.size() <= 8);
    Scalar field[8];
    for (unsigned idx = 0; idx < 8; ++idx)
        field[idx] = m_isoFunc(cl[idx]);

    Scalar grad[8][3];
    if (m_computeNormals) {
        const auto &hexahedronIndices = StructuredGridBase::HexahedronIndices;
        auto cellCoords = vistle::StructuredGridBase::cellCoordinates(CellNr, m_nvert);
        for (int idx = 0; idx < 8; ++idx) {
            Index x[3], xl[3], xu[3];
            for (int c = 0; c < 3; ++c)
                x[c] = cellCoords[c] + hexahedronIndices[c][idx];

            for (int c = 0; c < 3; ++c) {
                xl[c] = x[c] > 0 ? x[c] - 1 : x[c];
                xu[c] = x[c] < m_nvert[c] - 1 ? x[c] + 1 : x[c];
            }
            for (int c = 0; c < 3; ++c) {
                Index xx = x[c];
                x[c] = xl[c];
                Index l = StructuredGridBase::vertexIndex(x, m_nvert);
                x[c] = xu[c];
                Index u = StructuredGridBase::vertexIndex(x, m_nvert);
                x[c] = xx;
                grad[idx][c] = m_isoFunc(u) - m_isoFunc(l);
                Scalar diff = 0;

                if (m_haveCoords)
                    diff = (m_vertData.inVecPtr[3 + c][u] - m_vertData.inVecPtr[3 + c][l]);
                else
                    diff = (m_vertData.inVecPtr[3 + c][xu[c]] - m_vertData.inVecPtr[3 + c][xl[c]]);

                if (fabs(diff) > EPSILON)
                    grad[idx][c] /= diff;
                else
                    grad[idx][c] = 0;
            }
        }
    }

    if (m_haveCoords) {
        for (Index idx = 0; idx < m_numVertices[ValidCellIndex]; idx++) {
            interpolate(outCellIdx, field, hexaTriTable, hexaEdgeTable, interpolEdgeIdx1, interpolEdgeIdx2, t, 3);
            /* INTER(3, hexaTriTable, hexaEdgeTable); */
            if (m_computeNormals)
                for (int j = 0; j < 3; j++)
                    m_vertData.outVecPtr[j][outvertexindex] = lerp(grad[v1][j], grad[v2][j], t);
        }
    } else {
        for (Index idx = 0; idx < m_numVertices[ValidCellIndex]; idx++) {
            interpolate(outCellIdx, field, hexaTriTable, hexaEdgeTable, interpolEdgeIdx1, interpolEdgeIdx2, t, 6);
            /* INTER(6, hexaTriTable, hexaEdgeTable); */

            if (m_computeNormals)
                for (int j = 0; j < 3; j++)
                    m_outVertPtr[j][outvertexindex] = lerp(grad[v1][j], grad[v2][j], t);

            auto vc1 = StructuredGridBase::vertexCoordinates(cl[v1], m_nvert);
            auto vc2 = StructuredGridBase::vertexCoordinates(cl[v2], m_nvert);
            for (int j = 0; j < 3; j++)
                m_outVertPtr[3 + j][outvertexindex] =
                    lerp(m_data.m_inVertPtr[3 + j][vc1[j]], m_data.m_inVertPtr[3 + j][vc2[j]], t);
        }
    }
}

/**
 * @brief Find all iso-points on each edge of each face,
 *  build a triangle for each consecutive pair and a center point,
 *  orient outwards towards smaller values
 *
 * @param Cellbegin
 * @param Cellend
 * @param ValidCellIndex
 */
void HostData::computeUnstrPolyhedronOutput(Index Cellbegin, Index Cellend, Index ValidCellIndex)
{
    const auto &cl = m_cl;
    const Index numVert = m_numVertices[ValidCellIndex];
    Index numAvg = 0;
    MiddleData midData;

    Index outIdx = m_LocationList[ValidCellIndex];
    Index facestart = InvalidIndex; // CPOLYHEDRON
    Index term = 0; //CPOLYHEDRON
    // COMMENTS => for VTKPOLYHEDRON
    /* for (Index i = Cellbegin; i < Cellend; i += cl[i] + 1) { */
    for (Index i = Cellbegin; i < Cellend; ++i) {
        if (facestart == InvalidIndex) { // CPOLYHEDRON
            facestart = i;
            term = cl[i];
        } else if (term == cl[i]) { // CPOLYHEDRON
            /* const Index nvert = cl[i]; */
            const Index nvert = i - facestart;
            /* Index c1 = cl[i + nvert]; */
            bool flipped = false, haveIsect = false;
            /* for (Index k = i + 1; k < i + nvert + 1; ++k) { */
            for (Index k = facestart; k < facestart + nvert; ++k) {
                const Index c1 = cl[k];
                /* const Index c2 = cl[k]; */
                const Index c2 = cl[k + 1];

                Scalar d1 = m_isoFunc(c1);
                Scalar d2 = m_isoFunc(c2);

                bool smallToBig = d1 <= m_isovalue && d2 > m_isovalue;
                bool bigToSmall = d1 > m_isovalue && d2 <= m_isovalue;

                if (smallToBig || bigToSmall) {
                    if (!haveIsect) {
                        flipped = bigToSmall;
                        haveIsect = true;
                    }
                    Index out = outIdx;
                    if (flipped) {
                        if (bigToSmall)
                            out += 1;
                        else
                            out -= 1;
                    }
                    Scalar t = tinterp(m_isovalue, d1, d2);
                    m_vertDataContainer.copyContLerpToOutPerElem(m_numInVertData, m_numInVertDataI, m_numInVertDataB,
                                                                 midData, out, c1, c2, t);
                    ++outIdx;
                    if (bigToSmall ^ flipped)
                        ++outIdx;
                    ++numAvg;
                }

                facestart = InvalidIndex; //CPOLYHEDRON
                /* c1 = c2; */
            }
        }
    }
    midData(numAvg);
    for (Index i = 2; i < numVert; i += 3) {
        const Index idx = m_LocationList[ValidCellIndex] + i;
        m_vertDataContainer.copyMidDataToContOutPerElem(m_numInVertData, m_numInVertDataI, m_numInVertDataB, midData,
                                                        idx);
    }
}

template<size_t N, size_t M>
void HostData::setGhostLayers(const Index (&ghost)[N][M])
{
    for (int c = 0; c < 3; ++c)
        for (int i = 0; i < 2; ++i)
            m_nghost[c][i] = ghost[c][i];
}

template<typename T>
void HostData::addData(VecData<T> &vD, const T *data)
{
    vD.inVecPtr.push_back(data);
    vD.outVecData.emplace_back(ShmVector<T>::create(0));
    vD.outVecPtr.push_back(NULL);
    m_numInVertData = vD.inVecPtr.size();
}

template<typename T>
void HostData::addMappedData(const T *data)
{
    auto &vecData = get_vec_vert(*data);
    addData(vecData, data);
}

template<typename T>
void HostData:.addCellData(const T *data)
{
    auto &vecData = get_vec_cell(*data);
    addData(vecData, data);
}

void HostData::prepareField(Scalar *field, const Index fieldSize, const Index Cellbegin)
{
    const bool cl = m_cl == nullptr;

    // if performance is bad => write 2 for loops with m_cl check before loop.
    for (int idx = 0; idx < fieldSize; ++idx)
        field[idx] = m_isoFunc(cl ? m_cl[idx] : Cellbegin + idx);
}

void HostData::interpolate(const Index outVertIdx, const Index edgeIdx1, const Index edgeIdx2, const Scalar &t,
                           const Index nc = 0)
{
    m_vertDataContainer.copyContLerpToOutPerElem(m_numInVertData, m_numInVertDataI, m_numInVertDataB, outVertIdx,
                                                 edgeIdx1, edgeIdx2, t, nc);
}

//***************************************** VecDataContainer *****************************************//

void HostData::VecDataContainer::prepareContForOutput(const int numScalar, const int numIdx, const int numByte)
{
    /* std::apply([&](auto &&...x) { ((x.prepareVecDataForOutput(numScalar, numIdx, numByte)), ...); }, container_data); */
    prepareVecDataForOutput(scalarData, numScalar);
    prepareVecDataForOutput(indexData, numIdx);
    prepareVecDataForOutput(byteData, numByte);
}

void HostData::VecDataContainer::copyContInToOutPerElem(const int numScalar, const int numIdx, const int numByte,
                                                        const Index outCellIdx, const Index inCellIdx)
{
    copyInToOutPerElem(scalarData, numScalar, outCellIdx, inCellIdx);
    copyInToOutPerElem(indexData, numIdx, outCellIdx, inCellIdx);
    copyInToOutPerElem(byteData, numByte, outCellIdx, inCellIdx);
}

void HostData::VecDataContainer::copyContLerpToOutPerElem(const int numScalar, const int numIdx, const int numByte,
                                                          const Index outCellIdx, const Index edgeIdx1,
                                                          const Index edgeIdx2, const Scalar &t, const Index nc = 0)
{
    copyLerpToOutPerElem(scalarData, numScalar, outCellIdx, edgeIdx1, edgeIdx2, t, nc);
    copyLerpToOutPerElem(indexData, numIdx, outCellIdx, edgeIdx1, edgeIdx2, t);
    copyLerpToOutPerElem(byteData, numByte, outCellIdx, edgeIdx1, edgeIdx2, t);
}

void HostData::VecDataContainer::copyContLerpToOutPerElem(const int numScalar, const int numIdx, const int numByte,
                                                          MiddleData &midData, const Index outCellIdx,
                                                          const Index edgeIdx1, const Index edgeIdx2, const Scalar &t)
{
    copyLerpToOutPerElem(scalarData, numScalar, outCellIdx, edgeIdx1, edgeIdx2, t, midData.middleData.data());
    copyLerpToOutPerElem(indexData, numIdx, outCellIdx, edgeIdx1, edgeIdx2, t, midData.middleDataIdx.data());
    copyLerpToOutPerElem(byteData, numByte, outCellIdx, edgeIdx1, edgeIdx2, t, midData.middleDataByte.data());
}

void HostData::VecDataContainer::copyMidDataToContOutPerElem(const int numScalar, const int numIdx, const int numByte,
                                                             MiddleData &midData, const Index outCellIdx)
{
    copyMidDataToOutPerElem(scalarData, numScalar, outCellIdx, midData.middleData.data());
    copyMidDataToOutPerElem(indexData, numIdx, outCellIdx, midData.middleDataIdx.data());
    copyMidDataToOutPerElem(byteData, numByte, outCellIdx, midData.middleDataByte.data());
}

template<typename T>
void HostData::VecDataContainer::prepareVecDataForOutput(VecData<T> &vD, const int numInData)
{
    vD.prepareVecDataForOutput(numInData);
}

template<typename T>
void HostData::VecDataContainer::copyInToOutPerElem(VecData<T> &vD, const int numInData, const Index outCellIdx,
                                                    const Index inCellIdx)
{
    vD.copyInToOutPerElem(numInData, outCellIdx, inCellIdx);
}

template<typename T>
void HostData::VecDataContainer::copyLerpToOutPerElem(VecData<T> &vD, const int numInData, const Index outCellIdx,
                                                      const Index edgeIdx1, const Index edgeIdx2, const Scalar &t,
                                                      const Index nc = 0)
{
    vD.copyLerpToOutPerElem(numInData, outCellIdx, edgeIdx1, edgeIdx2, t, nc);
}

template<typename T>
void HostData::VecDataContainer::copyMidDataToOutPerElem(VecData<T> &vD, const int numInData, const Index outCellIdx,
                                                         T *midData)
{
    vD.copyContToOutPerElem(numInData, outCellIdx, midData);
}

template<typename T>
void HostData::VecDataContainer::copyLerpToOutPerElem(VecData<T> &vD, const int numInData, const Index outCellIdx,
                                                      const Index edgeIdx1, const Index edgeIdx2, const Scalar &t,
                                                      T *middleData)
{
    for (int i = 0; i < numInData; ++i) {
        auto v = vD.lerpIn(i, edgeIdx1, edgeIdx2, t);
        middleData[i] += v;
        vD.outVecPtr[i][outCellIdx] = v;
    }
}
