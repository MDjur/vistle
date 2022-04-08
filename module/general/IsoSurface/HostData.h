#ifndef HOSTDATA_H
#define HOSTDATA_H

#include "IsoDataFunctor.h"
#include "data_interface.h"
#include "tables.h"
#include <algorithm>
#include <array>
#include <tuple>
#include <vector>
#include <vistle/core/index.h>
#include <vistle/core/scalar.h>
#include <vistle/core/unstr.h>
#include <vistle/core/triangles.h>
#include <vistle/core/lines.h>
#include <vistle/core/shm.h>

using namespace vistle;

class HostData: public IsoDataI {
public:
    typedef const Byte *TypeIterator;
    typedef const Index *IndexIterator;
    typedef std::vector<Index>::iterator VectorIndexIterator;

    HostData(Scalar isoValue, IsoDataFunctor isoFunc, const Scalar *x, const Scalar *y, const Scalar *z);

    //unstructured
    HostData(Scalar isoValue, IsoDataFunctor isoFunc, const Index *el, const Byte *tl, const Index *cl, const Scalar *x,
             const Scalar *y, const Scalar *z);

    //structured
    HostData(Scalar isoValue, IsoDataFunctor isoFunc, Index nx, Index ny, Index nz, const Scalar *x, const Scalar *y,
             const Scalar *z);

    //poly
    HostData(Scalar isoValue, IsoDataFunctor isoFunc, const Index *el, const Index *cl, const Scalar *x,
             const Scalar *y, const Scalar *z);

    //triangles/quads
    HostData(Scalar isoValue, IsoDataFunctor isoFunc, int vertPerCell, const Index *cl, const Scalar *x,
             const Scalar *y, const Scalar *z);

    template<size_t xT, size_t yT, size_t xE, size_t yE>
    void computeGeoOutput(const Index Cellbegin, const Index fieldSize, const Index ValidCellIndex, const Index nc,
                          const int (&triTable)[xT][yT], const int (&edgeTable)[xE][yE]);
    void computeStructOutput(const Index CellNr, const Index ValidCellIndex);
    void computeUnstrPolyhedronOutput(Index Cellbegin, Index Cellend, Index ValidCellIndex);
    void copyCellInToOutAtIdx(Index ValidCellIndex, const Index CellNr);
    void prepareForOutput();
    void setComputeNormals(bool val) { m_computeNormals = val; }
    void setHaveCoords(bool val) { m_haveCoords = val; }
    template<size_t N, size_t M>
    void setGhostLayers(const Index (&ghost)[N][M]);
    template<typename T>
    void addMappedData(const T *data);
    template<typename T>
    void addCellData(const T *data);

    auto &NumVertPerCell() { return m_numVertPerCell; }

private:
    static constexpr int MaxNumData = 6;
    struct MiddleData {
        std::array<Scalar, MaxNumData> middleData;
        std::array<Index, MaxNumData> middleDataIdx;
        std::array<Byte, MaxNumData> middleDataByte;

        void operator()(Index numAvg)
        {
            if (numAvg > 0) {
                for (size_t i = 0; i < middleData.size(); ++i)
                    middleData[i] /= numAvg;
                for (int i = 0; i < middleDataIdx.size(); ++i)
                    middleDataIdx[i] /= numAvg;
                for (int i = 0; i < middleDataByte.size(); ++i)
                    middleDataByte[i] /= numAvg;
            }
        }
    };

    struct VecDataContainer {
        VecDataContainer(VecData<Scalar> &scaD, VecData<Index> &idxD, VecData<Byte> &bytD)
        : scalarData(scaD), indexData(idxD), byteData(bytD)
        /* : container_data(scaD, idxD, bytD) */ {}

        void prepareContForOutput(const int numScalar, const int numIdx, const int numByte);
        void copyContInToOutPerElem(const int numScalar, const int numIdx, const int numByte, const Index outCellIdx,
                                    const Index inCellIdx);

        void copyContLerpToOutPerElem(const int numScalar, const int numIdx, const int numByte, const Index outCellIdx,
                                      const Index edgeIdx1, const Index edgeIdx2, const Scalar &t);

        void copyContLerpToOutPerElem(const int numScalar, const int numIdx, const int numByte, MiddleData &midData,
                                      const Index outCellIdx, const Index edgeIdx1, const Index edgeIdx2,
                                      const Scalar &t);

        void copyMidDataToContOutPerElem(const int numScalar, const int numIdx, const int numByte, MiddleData &midData,
                                         const Index outCellIdx);

    private:
        template<typename T>
        void prepareVecDataForOutput(VecData<T> &vD, const int numInData);

        template<typename T>
        void copyInToOutPerElem(VecData<T> &vD, const int numInData, const Index outCellIdx, const Index inCellIdx);

        template<typename T>
        void copyLerpToOutPerElem(VecData<T> &vD, const int numInData, const Index outCellIdx, const Index edgeIdx1,
                                  const Index edgeIdx2, const Scalar &t);

        template<typename T>
        void copyMidDataToOutPerElem(VecData<T> &vD, const int numInData, const Index outCellIdx, T *midData);

        template<typename T>
        void copyLerpToOutPerElem(VecData<T> &vD, const int numInData, const Index outCellIdx, const Index edgeIdx1,
                                  const Index edgeIdx2, const Scalar &t, T *middleData);


        /* const std::tuple<VecData<Scalar> &, VecData<Index> &, VecData<Byte> &> container_data; */
        VecData<Scalar> &scalarData;
        VecData<Index> &indexData;
        VecData<Byte> &byteData;
    };

    template<typename T>
    void addData(VecData<T> &vD, const T *data);

    void prepareField(Scalar *field, const Index fieldSize, const Index Cellbegin);
    void interpolate(const Index outVertexIdx, const Index edgeIdx1, const Index edgeIdx2, const Scalar &t);

    const Index *m_el = nullptr;
    const Index *m_cl = nullptr;
    const Byte *m_tl = nullptr;
    std::vector<Index> m_caseNums;
    std::vector<Index> m_numVertices;
    std::vector<Index> m_LocationList;
    std::vector<Index> m_SelectedCellVector;
    int m_numVertPerCell = 0;
    VecDataContainer m_vertDataContainer;
    VecDataContainer m_cellDataContainer;
};
#endif
