#ifndef DATA_INTERFACE_H
#define DATA_INTERFACE_H

#include <vector>
#include "IsoDataFunctor.h"
#include <vistle/core/shm.h>
#include <vistle/core/index.h>
#include <vistle/core/scalar.h>

using namespace vistle;

class IsoDataI {
public:
    IsoDataI(Scalar isoValue, IsoDataFunctor isoFunc)
    : m_isovalue(isoValue), m_isoFunc(isoFunc), m_nghost{{0, 0}, {0, 0}, {0, 0}}, m_computeNormals(false)
    {}

    //this implementation of VecData is used in HostData
    template<typename T>
    struct VecData {
        std::vector<const T *> inVecPtr;
        std::vector<vistle::shm_array_ref<vistle::shm_array<T, typename vistle::shm<T>::allocator>>> outVecData;
        std::vector<T *> outVecPtr;

        void prepareVecDataForOutput(const int numInData)
        {
            for (int i = 0; i < numInData; i++)
                outVecPtr[i] = outVecData[i]->data();
        }

        void copyInToOutPerElem(const int numInData, const Index outCellIdx, const Index inCellIdx)
        {
            for (int i = 0; i < numInData; ++i)
                outVecPtr[i][outCellIdx] = inVecPtr[i][inCellIdx];
        }

        void copyContToOutPerElem(const int numInData, const Index outCellIdx, T *container)
        {
            for (int i = 0; i < numInData; ++i)
                outVecPtr[i][outCellIdx] = container[i];
        }

        void copyLerpToOutPerElem(const int numInData, const Index outCellIdx, const Index edgeIdx1,
                                  const Index edgeIdx2, const Scalar &t, const Index nc = 0)
        {
            for (int i = nc; i < numInData; ++i)
                outVecPtr[i][outCellIdx] = lerpIn(i, edgeIdx1, edgeIdx2, t);
        }

        auto lerpIn(const int idx, const Index edgeIdx1, const Index edgeIdx2, const Scalar &t) const
        {
            return lerp(inVecPtr[idx][edgeIdx1], inVecPtr[idx][edgeIdx2], t);
        }
    };

    auto &Isovalue() { return m_isovalue; }
    auto &NumInVertData() { return m_numInVertData; }
    auto &NumInVertDataI() { return m_numInVertDataI; }
    auto &NumInVertDataB() { return m_numInVertDataB; }
    auto &NumInCellData() { return m_numInCellData; }
    auto &NumInCellDataI() { return m_numInCellDataI; }
    auto &NumInCellDataB() { return m_numInCellDataB; }
    auto &IsoFunc() { return m_isoFunc; }
    auto &El() { return m_el; }
    auto &Cl() { return m_cl; }
    auto &Tl() { return m_tl; }
    auto &CaseNums() { return m_caseNums; }
    auto &NumVertices() { return m_numVertices; }
    auto &LocationList() { return m_LocationList; }
    auto &SelectedCellVector() { return m_SelectedCellVector; }
    auto &Nvert() { return m_nvert; }
    auto &Nghost() { return m_nghost; }
    auto &VertData() { return m_vertData; }
    auto &CellData() { return m_cellData; }
    auto &VertDataIdx() { return m_vertDataIdx; }
    auto &CellDataIdx() { return m_cellDataIdx; }
    auto &VertDataByte() { return m_vertDataByte; }
    auto &CellDataByte() { return m_cellDataByte; }
    auto &IsUnstructured() { return m_isUnstructured; }
    auto &IsPoly() { return m_isPoly; }
    auto &IsTri() { return m_isTri; }
    auto &IsQuad() { return m_isQuad; }
    auto &HaveCoords() { return m_haveCoords; }
    auto &ComputeNormals() { return m_computeNormals; }

    template<typename T>
    VecData<T> &get_vec_vert(const T &);

    template<typename T>
    VecData<T> &get_vec_cell(const T &);

protected:
    template<>
    VecData<Scalar> &get_vec_vert(const Scalar &)
    {
        return m_vertData;
    }

    template<>
    VecData<Index> &get_vec_vert(const Index &)
    {
        return m_vertDataIdx;
    }

    template<>
    VecData<Byte> &get_vec_vert(const Byte &)
    {
        return m_vertDataByte;
    }

    template<>
    VecData<Scalar> &get_vec_cell(const Scalar &)
    {
        return m_cellData;
    }

    template<>
    VecData<Index> &get_vec_cell(const Index &)
    {
        return m_cellDataIdx;
    }

    template<>
    VecData<Byte> &get_vec_cell(const Byte &)
    {
        return m_cellDataByte;
    }

    Scalar m_isovalue;
    int m_numInVertData = 0, m_numInVertDataI = 0, m_numInVertDataB = 0;
    int m_numInCellData = 0, m_numInCellDataI = 0, m_numInCellDataB = 0;
    IsoDataFunctor m_isoFunc;
    const Index *m_el = nullptr;
    const Index *m_cl = nullptr;
    const Byte *m_tl = nullptr;
    std::vector<Index> m_caseNums;
    std::vector<Index> m_numVertices;
    std::vector<Index> m_LocationList;
    std::vector<Index> m_SelectedCellVector;
    Index m_nvert[3];
    Index m_nghost[3][2];
    VecData<Scalar> m_vertData, m_cellData;
    VecData<Index> m_vertDataIdx, m_cellDataIdx;
    VecData<Byte> m_vertDataByte, m_cellDataByte;
    bool m_isUnstructured = false;
    bool m_isPoly = false;
    bool m_isTri = false;
    bool m_isQuad = false;
    bool m_haveCoords = false;
    bool m_computeNormals = false;
};
#endif
