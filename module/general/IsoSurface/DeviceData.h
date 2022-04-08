#ifndef DEVICEDATA_H
#define DEVICEDATA_H

#include <thrust/device_ptr.h>
#include <thrust/device_vector.h>

#include "IsoDataFunctor.h"
#include "data_interface.h"
#include <vistle/core/index.h>
#include <vistle/core/scalar.h>

using namespace vistle;

class DeviceData: IsoDataI {
public:
    typedef thrust::device_vector<Index>::iterator IndexIterator;
    typedef thrust::device_vector<Byte>::iterator TypeIterator;
    typedef thrust::device_vector<Index>::iterator VectorIndexIterator;

    //reimplement of DataI
    template<typename T>
    struct VecData {
        std::vector<thrust::device_ptr<T>> inVecPtr;
        std::vector<thrust::device_vector<T> *> outVecData;
        std::vector<thrust::device_ptr<T>> outVecPtr;
    };

    DeviceData(Scalar isoValue, IsoDataFunctor isoFunc, Index nelem, const Index *el, const Byte *tl, Index nconn,
               const Index *cl, Index ncoord, const Scalar *x, const Scalar *y, const Scalar *z)
    : IsoDataI(isoValue, isoFunc)
    , m_el(el, el + nelem)
    , m_cl(cl, cl + nconn)
    , m_tl(tl, tl + nelem)
    , m_x(x, x + ncoord)
    , m_y(y, y + ncoord)
    , m_z(z, z + ncoord)
    {
        m_nvert[0] = 0;
        m_nvert[1] = 0;
        m_nvert[2] = 0;
        m_haveCoords = true;
        m_isUnstructured = el;

        m_vertData.inVecPtr.push_back(m_x.data());
        m_vertData.inVecPtr.push_back(m_y.data());
        m_vertData.inVecPtr.push_back(m_z.data());

        initData(m_vertData, m_numInVertData);
        initData(m_vertDataIdx, m_numInVertDataI);
        initData(m_vertDataByte, m_numInVertDataB);
    }

private:
    template<typename T>
    void initData(VecData<T> &vD, int &numInVert)
    {
        auto size = vD.inVecPtr.size();
        for (size_t i = 0; i < size; i++)
            vD.outVecData.push_back(new thrust::device_vector<T>);
        vD.outVecPtr.resize(size);
        numInVert = size;
    }

    thrust::device_vector<Index> m_el;
    thrust::device_vector<Index> m_cl;
    thrust::device_vector<Byte> m_tl;
    thrust::device_vector<Index> m_caseNums;
    thrust::device_vector<Index> m_numVertices;
    thrust::device_vector<Index> m_LocationList;
    thrust::device_vector<Index> m_SelectedCellVector;
    thrust::device_vector<Scalar> m_x;
    thrust::device_vector<Scalar> m_y;
    thrust::device_vector<Scalar> m_z;
    VecData<Scalar> m_vertData, m_cellData;
    VecData<Index> m_vertDataIdx, m_cellDataIdx;
    VecData<Byte> m_vertDataByte, m_cellDataByte;
};
#endif
