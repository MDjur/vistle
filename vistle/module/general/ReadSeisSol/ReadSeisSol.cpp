/**************************************************************************\
 **                                                                      **
 **                                                                      **
 ** Description: Read module for ChEESE Seismic HDF5/XDMF-files  	     **
 **                                                                      **
 **                                                                      **
 **                                                                      **
 **                                                                      **
 **                                                                      **
 ** Author:    Marko Djuric                                              **
 **                                                                      **
 **                                                                      **
 **                                                                      **
 ** Date:  04.03.2021                                                    **
\**************************************************************************/

//*************************** TODO: IMPLEMENT XDMF-READER OUT OF THIS*****************************//

//header
#include "ReadSeisSol.h"

//mpi
#include <mpi.h>

//vistle
#include "vistle/core/object.h"
#include "vistle/core/port.h"
#include "vistle/core/database.h"
#include "vistle/core/parameter.h"
#include "vistle/core/scalar.h"
#include "vistle/core/unstr.h"
#include "vistle/core/vec.h"
#include "vistle/module/module.h"
#include "vistle/module/reader.h"
#include "vistle/util/enum.h"

//boost
#include <boost/mpi/communicator.hpp>
#include <boost/smart_ptr/shared_ptr.hpp>

//xdmf3
#include <XdmfAttribute.hpp>
#include <XdmfTopology.hpp>
#include <XdmfAttributeCenter.hpp>
#include <XdmfAttributeType.hpp>
#include <XdmfTime.hpp>
#include <XdmfGeometry.hpp>
#include <XdmfGridCollection.hpp>
#include <XdmfUnstructuredGrid.hpp>
#include <XdmfDomain.hpp>
#include <XdmfReader.hpp>
#include <XdmfArray.hpp>
#include <XdmfArrayType.hpp>
#include <XdmfHDF5Controller.hpp>

//std
#include <iterator>
#include <array>
#include <functional>
#include <cstring>
#include <iostream>
#include <map>
#include <memory>
#include <numeric>
#include <string>
#include <algorithm>
#include <string_view>
#include <type_traits>
#include <vector>

using namespace vistle;

MODULE_MAIN(ReadSeisSol)

namespace {

constexpr std::string_view ERROR_GRID_GENERATION = "Something went wrong while creation of ugrid! Possible errors:\n"
                                                   "-invalid connectivity list \n"
                                                   "-invalid coordinates\n"
                                                   "-invalid typelist\n"
                                                   "-invalid geo type";

constexpr std::string_view ERROR_ATT_GENERATION =
    "Something went wrong while creation of attribute pointer! Possible errors:\n"
    "-wrong format (hyperslap not supported at the moment) \n"
    "-empty attribute\n"
    "-wrong dimensions defined";

//SeisSol Mode
DEFINE_ENUM_WITH_STRING_CONVERSIONS(SeisSolMode, (XDMF)(HDF5))

//Parallel Mode
DEFINE_ENUM_WITH_STRING_CONVERSIONS(ParallelMode, (BLOCKS)(SERIAL))

/**
  * @brief: Extracts XdmfArrayType from ugrid and returns corresponding vistle UnstructuredGrid::Type.
  *
  * @ugrid: Pointer to XdmfUnstructuredGrid.
  *
  * @return: Corresponding vistle type.
  */
const UnstructuredGrid::Type XdmfUgridToVistleUgridType(const XdmfUnstructuredGrid *ugrid)
{
    const auto &ugridType = ugrid->getTopology()->getType();
    const std::string &topologyType = ugridType->getName();
    if (topologyType == "Tetrahedron")
        return UnstructuredGrid::TETRAHEDRON;
    else if (topologyType == "Triangle")
        return UnstructuredGrid::TRIANGLE;
    else if (topologyType == "Polygon")
        return UnstructuredGrid::POLYGON;
    else if (topologyType == "Pyramid")
        return UnstructuredGrid::PYRAMID;
    else if (topologyType == "Hexahedron")
        return UnstructuredGrid::HEXAHEDRON;
    else if (topologyType == "Polyvertex")
        return UnstructuredGrid::POINT;
    else if (topologyType == "Polyline")
        return UnstructuredGrid::BAR;
    else if (topologyType == "Quadrilateral")
        return UnstructuredGrid::QUAD;
    else if (topologyType == "Wedge")
        return UnstructuredGrid::PRISM;
    else
        return UnstructuredGrid::NONE;
}

/**
 * @brief: Checks if suffix string is the end of main.
 *
 * @param main: main string.
 * @param suffix: suffix string
 *
 * @return: True if given suffix is end of main.
 */
bool checkFileEnding(const std::string &main, const std::string &suffix)
{
    if (main.size() >= suffix.size())
        if (main.compare(main.size() - suffix.size(), suffix.size(), suffix) == 0)
            return true;
    return false;
}

/**
 * @brief: Check for nullptr. Send msgFailure if vObj_ptr is equal nullptr.
 *
 * @param obj: ReadSeisSol object pointer.
 * @param vObj_ptr: Object to check.
 * @param msgFailure: Error message.
 * @brief Release and reset XdmfArray pointer.
 *
 * @return: False if vObj_ptr is a nullptr.
 */
bool checkObjectPtr(ReadSeisSol *obj, vistle::Object::ptr vObj_ptr, const std::string_view &msgFailure)
{
    if (!vObj_ptr) {
        obj->sendInfo("%s", msgFailure.data());
        return false;
    }
    return true;
}

/**
 * @brief Extract unique points from xdmfarray and store them in a set.
 *
 * @tparam T Set generic parameter.
 * @param arr XdmfArray with values to extract.
 * @param refSet Storage std::set.
 */
template<class T>
void extractUniquePoints(XdmfArray *arr, std::set<T> &refSet)
{
    if (arr)
        for (unsigned i = 0; i < arr->getSize(); ++i)
            refSet.insert(arr->getValue<T>(i));
}
} // namespace

/**
  * @brief: Constructor.
  */
ReadSeisSol::ReadSeisSol(const std::string &name, int moduleID, mpi::communicator comm)
: vistle::Reader("Read ChEESE Seismic Data files (xmdf/hdf5)", name, moduleID, comm)
{
    //settings
    m_file = addStringParameter("xfile", "Xdmf File", "/data/ChEESE/SeisSol/LMU_Sulawesi_example/sula-surface.xdmf",
                                Parameter::Filename);
    m_seisMode = addIntParameter("SeisSolMode", "Select file format", (Integer)XDMF, Parameter::Choice);
    V_ENUM_SET_CHOICES(m_seisMode, SeisSolMode);
    m_parallelMode =
        addIntParameter("ParallelMode", "Select ParallelMode (Block or Serial (with and without timestepdistribution))",
                        (Integer)BLOCKS, Parameter::Choice);
    V_ENUM_SET_CHOICES(m_parallelMode, ParallelMode);

    //standard ports
    m_gridOut = createOutputPort("ugrid", "UnstructuredGrid");

    //scalar
    initScalar();

    //partition (maybe relevant for HDF5)
    /* initBlocks(); */

    m_reuseGrid =
        addIntParameter("ReuseGrid", "Reuses grid from first XdmfGrid specified in Xdmf.", 1, Parameter::Boolean);

    //observer
    initObserveParameter();

    //parallel mode
    setParallelizationMode(ParallelizeBlocks);
    setPartitions(size());
}

/**
 * @brief Set parameter to observe.
 */
void ReadSeisSol::initObserveParameter()
{
    observeParameter(m_file);
    observeParameter(m_seisMode);
    observeParameter(m_parallelMode);
    observeParameter(m_reuseGrid);
    observeParameter(m_block);
}


/**
 * @brief Initialize block parameter.
 */
/* void ReadSeisSol::initBlocks() */
/* { */
/*     constexpr auto blockName_constexpr{"Blocks"}; */
/*     constexpr auto blockDescr_constexpr{"Number of blocks for partition."}; */

/*     m_block = addIntParameter(blockName_constexpr, blockDescr_constexpr, 1); */
/*     setParameterRange(m_block, Integer(1), Integer(999999)); */

/* const std::string block_constStr{"Number of blocks "}; */
/* const std::string descr_constStr{"Number of blocks in "}; */

/* char bSpecifier{'x'}; */
/* for (auto &block: m_blocks) { */
/*     const std::string &blockName = block_constStr + bSpecifier; */
/*     const std::string &blockDescription = descr_constStr + (bSpecifier++); */
/*     block = addIntParameter(blockName, blockDescription, 1); */
/*     setParameterRange(block, Integer(1), Integer(999999)); */
/*     observeParameter(block); */
/* } */
/* } */

/**
 * @brief Initialize scalar ports and choice parameters.
 */
void ReadSeisSol::initScalar()
{
    constexpr auto attName_constexpr{"Scalar "};
    constexpr auto portName_constexpr{"Scalar port "};
    constexpr auto portDescr_constexpr{"Scalar port for attribute "};

    for (Index i = 0; i < m_attributes.size(); i++) {
        const std::string i_str{std::to_string(i)};
        const std::string &attName = attName_constexpr + i_str;
        const std::string &portName = portName_constexpr + i_str;
        const std::string &portDescr = portDescr_constexpr + i_str;
        m_attributes[i] = addStringParameter(attName, "Select scalar for scalar port.", "", Parameter::Choice);
        m_scalarsOut[i] = createOutputPort(portName, portDescr);
        observeParameter(m_attributes[i]);
    }
}

/**
 * @brief Destructor.
 */
ReadSeisSol::~ReadSeisSol()
{
    releaseXdmfObjects();
}

/**
 * @brief Clear up allocated xdmf objects by reader.
 */
void ReadSeisSol::releaseXdmfObjects()
{
    if (!xgridCollect) {
        xgridCollect->release();
        xgridCollect.reset();
    }
    if (!m_vugrid_ptr)
        m_vugrid_ptr.reset();
}

/**
  * @brief: Template for switching between XDMF and HDF5 mode. Calls given function with args according to current mode.
  *
  * @xdmfFunc: Function that will be called when in XDMF mode.
  * @hdfFunc: Function that will be called when in HDF5 mode.
  *
  * return: seisFunc(args) -> return type.
  */
template<class Ret, class... Args>
auto ReadSeisSol::switchSeisMode(std::function<Ret(Args...)> xdmfFunc, std::function<Ret(Args...)> hdfFunc,
                                 Args... args)
{
    switch (m_seisMode->getValue()) {
    case XDMF:
        return xdmfFunc(args...);
    case HDF5:
        return hdfFunc(args...);
    }
    return false;
}

/**
  * @brief: Wrapper function for easier call of seismode-function for selected seismode.
  *
  * @xdmfFunc: Function pointer to XdfmFunction.
  * @hdfFunc: Function pointer to hdfFunction.
  * @args: Arguments to pass.
  *
  * @return: Return val of called function.
  */
template<class Ret, class... Args>
auto ReadSeisSol::callSeisModeFunction(Ret (ReadSeisSol::*xdmfFunc)(Args...), Ret (ReadSeisSol::*hdfFunc)(Args...),
                                       Args... args)
{
    using funcType = Ret(ReadSeisSol *, Args...);
    std::function<funcType> xFunc{xdmfFunc};
    std::function<funcType> hFunc{hdfFunc};
    return switchSeisMode<Ret, ReadSeisSol *, Args...>(xFunc, hFunc, this, args...);
}

/**
  * @brief: Helper function to tell user that HDF5 mode is not implemented.
  *
  * @return: false.
  */
bool ReadSeisSol::hdfModeNotImplemented()
{
    sendInfo("HDF5 mode is not implemented.");
    return false;
}
bool ReadSeisSol::hdfModeNotImplementedRead(Token &, int, int)
{
    return hdfModeNotImplemented();
}

/**
  * @brief: Clears attributes choice.
  */
void ReadSeisSol::clearChoice()
{
    for (const auto &att: m_attributes)
        att->setChoices(std::vector<std::string>());

    m_attChoiceStr.clear();
}

/**
  * @brief: Iterate through properties of a xdmfgridcollection.
  *
  * @xgridCol: Constant grid collection pointer.
  */
void ReadSeisSol::inspectXdmfGridCollection(const XdmfGridCollection *xgridCol)
{
    if (xgridCol) {
        const auto &gridColType = xgridCol->getType();
        if (gridColType == XdmfGridCollectionType::Temporal()) {
            //at the moment attribute number and type are the same for all ugrids.
            if (xgridCol->getNumberUnstructuredGrids()) {
                //number of unstructured grids = timestep
                setTimesteps(xgridCol->getNumberUnstructuredGrids());
                inspectXdmfUnstrGrid(xgridCol->getUnstructuredGrid(1).get());
            }
            if (xgridCol->getNumberCurvilinearGrids() | xgridCol->getNumberRectilinearGrids() |
                xgridCol->getNumberRegularGrids())
                sendInfo("Other gridtypes than unstructured have not been implemented yet!");
        } else if (gridColType == XdmfGridCollectionType::Spatial())
            sendInfo("Spatial collectiontype not implemented yet!");
        else
            sendInfo("No implementation without xdmfcollection.");
    }
}

/**
  * @brief: Iterate through unstructured grid attributes stored in an xdmf.
  *
  * @xugrid: const unstructured grid pointer.
  */
void ReadSeisSol::inspectXdmfUnstrGrid(const XdmfUnstructuredGrid *xugrid)
{
    if (xugrid) {
        //TODO: NOT IMPLEMENTED YET => FUNCTIONS DO NOTHING AT THE MOMENT.
        inspectXdmfTopology(xugrid->getTopology().get());
        inspectXdmfGeometry(xugrid->getGeometry().get());
        inspectXdmfTime(xugrid->getTime().get());

        clearChoice();
        for (unsigned i = 0; i < xugrid->getNumberAttributes(); ++i)
            inspectXdmfAttribute(xugrid->getAttribute(i).get());
    }
}

/**
 * @brief NOT IMPLEMENTED
 *
 * @param xgeo
 */
void ReadSeisSol::inspectXdmfGeometry(const XdmfGeometry *xgeo)
{
    if (xgeo) {
        // set geo global ?
    }
}

/**
 * @brief NOT IMPLEMENTED
 *
 * @param xtime
 */
void ReadSeisSol::inspectXdmfTime(const XdmfTime *xtime)
{
    if (xtime) {
        // set timestep defined in xdmf
    }
}

/**
 * @brief NOT IMPLEMENTED
 *
 * @param xtopo
 */
void ReadSeisSol::inspectXdmfTopology(const XdmfTopology *xtopo)
{
    if (xtopo) {
        // set topo global ?
    }
}

/**
  * @brief: Inspect XdmfAttribute and generate a parameter in reader.
  *
  * @xatt: Const attribute pointer.
  */
void ReadSeisSol::inspectXdmfAttribute(const XdmfAttribute *xatt)
{
    if (xatt)
        m_attChoiceStr.push_back(xatt->getName());
}

/**
  * @brief: Inspect XdmfDomain.
  *
  * @xdomain: Const pointer of domain-object.
  */
void ReadSeisSol::inspectXdmfDomain(const XdmfDomain *xdomain)
{
    if (xdomain)
        for (unsigned i = 0; i < xdomain->getNumberGridCollections(); ++i)
            inspectXdmfGridCollection(xdomain->getGridCollection(i).get());
}

/**
  * @brief: Inspect Xdmf.
  *
  * return: False if something went wrong while inspection else true.
  */
bool ReadSeisSol::inspectXdmf()
{
    const std::string &xfile = m_file->getValue();

    //read xdmf
    const auto xreader = XdmfReader::New();

    if (!xreader.get())
        return false;

    const auto &xdomain = shared_dynamic_cast<XdmfDomain>(xreader->read(xfile));
    inspectXdmfDomain(xdomain.get());

    for (auto &scalar: m_attributes)
        setParameterChoices(scalar, m_attChoiceStr);

    m_xAttSelect = DynamicPseudoEnum(m_attChoiceStr);

    return true;
}

//TIP: maybe interesting for HDF5 implementation
/**
 * @brief Checks that total product of blocks are equal mpi size.
 *
 * @return result of "product blocks == mpisize".
 */
/* bool ReadSeisSol::checkBlocks() */
/* { */
/* const int &nBlocks = m_blocks[0]->getValue() * m_blocks[1]->getValue() * m_blocks[2]->getValue(); */

//parallelization over ranks
/* setPartitions(nBlocks); */
/* if (nBlocks == size()) */
/*     return true; */
/* else { */
/*     if (rank() == 0) */
/*         sendInfo("Number of blocks should equal MPISIZE."); */
/*     return false; */
/* } */
/* } */

/**
  * @brief: Called if observed parameter changes.
  *
  * @param: changed parameter.
  *
  * @return: True if xdmf is in correct format else false.
  */
bool ReadSeisSol::examine(const vistle::Parameter *param)
{
    if (!param || param == m_file) {
        if (!checkFileEnding(m_file->getValue(), ".xdmf"))
            return false;

        callSeisModeFunction(&ReadSeisSol::inspectXdmf, &ReadSeisSol::hdfModeNotImplemented);
    } else if (param == m_parallelMode) {
        switch (m_parallelMode->getValue()) {
        case BLOCKS: {
            setParallelizationMode(ParallelizeBlocks);
            setPartitions(size());
            break;
        }
        case SERIAL: {
            setParallelizationMode(Serial);
            setAllowTimestepDistribution(true);
            setPartitions(0);
            break;
        }
        }
    }

    /* return checkBlocks(); */
    return true;
}

/**
  * @brief: Fill typelist of unstructured grid vistle pointer with equal type.
  *
  * @unstr: UnstructuredGrid::ptr that contains type list.
  * @type: UnstructuredGrid::Type.
  */
void ReadSeisSol::fillUnstrGridTypeList(vistle::UnstructuredGrid::ptr unstr, const vistle::UnstructuredGrid::Type &type)
{
    /* const auto &setGhostType{type | vistle::UnstructuredGrid::GHOST_BIT}; */
    std::fill(unstr->tl().begin(), unstr->tl().end(), type);
}

/**
  * @brief: Fill elementlist of unstructured grid vistle pointer.
  *
  * @unstr: UnstructuredGrid::ptr with elementlist to fill.
  * @numCornerPerElem: Number of Corners.
  */
template<class T>
void ReadSeisSol::fillUnstrGridElemList(vistle::UnstructuredGrid::ptr unstr, const T &numCornerPerElem)
{
    //e.g. 4 corner = tetrahedron
    std::generate(unstr->el().begin(), unstr->el().end(),
                  [n = 0, &numCornerPerElem]() mutable { return n++ * numCornerPerElem; });
}

/**
  * @brief: Fill coordinates stored in a XdmfArray into unstructured grid vistle pointer.
  *
  * @unstr: UnstructuredGrid::ptr with coordinates to fill.
  * @xArrGeo: XdmfArray which contains coordinates in an one dimensional array.
  *
  * @return: true if values stored in XdmfArray could be stored in coords array of given unstructured grid.
  */
bool ReadSeisSol::fillUnstrGridCoords(vistle::UnstructuredGrid::ptr unstr, XdmfArray *xArrGeo)
{
    if (!xArrGeo)
        return false;

    auto x = unstr->x().data(), y = unstr->y().data(), z = unstr->z().data();

    //TODO: not the same for all xdmf => hyperslap not working => maybe using XML parser and reading hyperslap to concretize read.
    // -> only 3D at the moment
    constexpr unsigned numCoords{3};
    constexpr unsigned strideVistleArr{1};

    //current order for geo-arrays is contiguous: arrGeo => x1 y1 z1 x2 y2 z2 x3 y3 z3 ... xn yn zn
    xArrGeo->getValues(0, x, xArrGeo->getSize() / numCoords, numCoords, strideVistleArr);
    xArrGeo->getValues(1, y, xArrGeo->getSize() / numCoords, numCoords, strideVistleArr);
    xArrGeo->getValues(2, z, xArrGeo->getSize() / numCoords, numCoords, strideVistleArr);

    return true;
}

/**
 * @brief Fill coordinates stored in a XdmfArray into unstructured grid vistle pointer.
 *
 * @param unstr UnstructuredGrid::ptr with coordinates to fill.
 * @param xArrGeo XdmfArray which contains coordinates in an one dimensional array.
 * @param verticesToRead verticies to read from array.
 *
 * @return true if values stored in XdmfArray could be stroed in coords array of given unstructured grid.
 */
bool ReadSeisSol::fillUnstrGridCoords(vistle::UnstructuredGrid::ptr unstr, XdmfArray *xArrGeo,
                                      const std::set<unsigned> &verticesToRead)
{
    if (!xArrGeo)
        return false;

    if (verticesToRead.empty())
        return fillUnstrGridCoords(unstr, xArrGeo);

    auto x = unstr->x().data(), y = unstr->y().data(), z = unstr->z().data();

    //fill only with point coords from vertToRead
    for (const auto &vert: verticesToRead) {
        x[vert] = xArrGeo->getValue<unsigned>(vert);
        y[vert] = xArrGeo->getValue<unsigned>(vert + 1);
        z[vert] = xArrGeo->getValue<unsigned>(vert + 2);
    }

    return true;
}

/**
  * @brief: Fill connectionlist stored in a XdmfArray into unstructured grid vistle pointer.
  *
  * @unstr: UnstructuredGrid::ptr with connectionlist to fill.
  * @xArrConn: XdmfArray which contains connectionlist in an one dimensional array.
  *
  * @return: true if connectionlist has been updated correctly.
  */
bool ReadSeisSol::fillUnstrGridConnectList(vistle::UnstructuredGrid::ptr unstr, XdmfArray *xArrConn)
{
    if (!xArrConn)
        return false;

    //TODO: Debug invert elementvertices for tetrahedron => vtk vertice order = inverse vistle vertice order
    auto cl = unstr->cl().data();
    xArrConn->getValues(0, cl, xArrConn->getSize());
    return true;
}

/**
  * @brief: Helper function for reading XdmfArray with given XdmfHeavyDataController.
  *
  * @xArr: XdmfArray for read operation.
  * @controller: XdmfHeavyDataController that specifies how to read the XdmfArray.
  */
void ReadSeisSol::readXdmfHeavyController(XdmfArray *xArr, const boost::shared_ptr<XdmfHeavyDataController> &controller)
{
    //has to be empty otherwise potential memoryleak
    if (xArr->isInitialized())
        xArr->release();

    xArr->insert(controller);
    xArr->read();
}

/**
 * @brief Read attributes stored in HDF5 file parallel with blocks into memory.
 *
 * @param xArrAtt XdmfArray pointer for storage.
 * @param attDim Dimension of attribute per timestep.
 * @param defaultController default heavy data controller.
 * @param block current block.
 * @param timestep current timestep.
 */
void ReadSeisSol::readXdmfHDF5AttributeParallel(XdmfArray *xArrAtt, unsigned &attDim,
                                               const boost::shared_ptr<XdmfHeavyDataController> &defaultController,
                                               const int block, const int timestep)
{
    const auto controller = shared_dynamic_cast<XdmfHDF5Controller>(defaultController);
    if (!controller)
        return;

    const auto &h5Path = controller->getFilePath();
    const auto &setPath = controller->getDataSetPath();
    const auto &arrayType = controller->getType();
    const auto &readStride = controller->getStride();

    attDim /= numPartitions();
    const auto &startVal = attDim * block;
    std::vector<unsigned> readStart{startVal};
    std::vector<unsigned> readCount{attDim};
    std::vector<unsigned> readDataSpace{attDim};

    if (controller->getDimensions().size() == 2) {
        readStart.insert(readStart.begin(), timestep);
        readCount.insert(readCount.begin(), 1);
        readDataSpace.insert(readDataSpace.begin(), 1);
    }

    auto hdfController =
        XdmfHDF5Controller::New(h5Path, setPath, arrayType, readStart, readStride, readCount, readDataSpace);

    xArrAtt->setHeavyDataController(hdfController);
    xArrAtt->read();

    hdfController->closeFiles();
}

/**
 * @brief Reads topology parallel and return unique vertices list.
 *
 * @param xArrTopo Array to store topology data.
 * @param defaultControllerTopo default HeavyDataController for topology.
 * @param block current block.
 */
void ReadSeisSol::readXdmfHDF5TopologyParallel(XdmfArray *xArr,
                                              const boost::shared_ptr<XdmfHeavyDataController> &defaultController,
                                              const int block)
{
    /*************** read topology parallel **************/
    auto controller = shared_dynamic_cast<XdmfHDF5Controller>(defaultController);
    if (!controller)
        return;

    const auto &h5Path = controller->getFilePath();
    const auto &setPath = controller->getDataSetPath();
    const auto &arrayType = controller->getType();
    const auto &xtopo = controller->getDimensions();
    const auto &readStride = controller->getStride();

    const auto &[numCorner, numElem] = std::minmax_element(xtopo.begin(), xtopo.end());
    const unsigned &countTopo{*numElem / numPartitions()};
    const std::vector<unsigned> readCount{countTopo, *numCorner};
    const std::vector<unsigned> readStart{countTopo * block, 0};
    const std::vector<unsigned> readDataSpace{countTopo, *numCorner};

    auto hdfController = XdmfHDF5Controller::New(h5Path, setPath, arrayType, readStart, readStride,
                                                 readCount, readDataSpace);

    xArr->setHeavyDataController(hdfController);
    xArr->read();

    hdfController->closeFiles();
}

/**
 * @brief: Generate a scalar pointer with values from XdmfAttribute.
 *
 * @param xattribute: XdmfAttribute pointer.
 * @param timestep: Current timestep.
 * @param block: Current block number.
 *
 * @return: generated vistle scalar pointer.
 */
vistle::Vec<Scalar>::ptr ReadSeisSol::generateScalarFromXdmfAttribute(XdmfAttribute *xattribute, const int &timestep,
                                                                      const int &block)
{
    if (!xattribute)
        return nullptr;

    const auto &xattDimVec = xattribute->getDimensions();
    const auto &numDims = xattDimVec.size();
    unsigned attDim = *std::max_element(xattDimVec.begin(), xattDimVec.end());

    unsigned startRead{0};
    constexpr unsigned arrStride{1};
    constexpr unsigned valStride{1};
    if (numDims == 2)
        startRead = attDim * timestep;

    const shared_ptr<XdmfArray> xArrAtt(XdmfArray::New());
    const auto &mode = m_parallelMode->getValue();
    const auto &controller = xattribute->getHeavyDataController();

    if (mode == BLOCKS) {
        if (XdmfAttributeCenter::Cell() == xattribute->getCenter())
            readXdmfHDF5AttributeParallel(xArrAtt.get(), attDim, controller, block, timestep);
        else  {
            sendInfo("Other options than centering attributes according to cells in block mode are not implemented");
            return nullptr;
        }
        startRead = 0;
    } else if (mode == SERIAL)
        readXdmfHeavyController(xArrAtt.get(), controller);
    else
        return nullptr;

    if (xArrAtt->getSize() == 0)
        return nullptr;

    Vec<Scalar>::ptr att(new Vec<Scalar>(attDim));
    auto vattDataX = att->x().data();
    xArrAtt->getValues<Scalar>(startRead, vattDataX, attDim, arrStride, valStride);
    return att;
}

/**
  * @brief: Generate an unstructured grid from XdmfUnstructuredGrid.
  *
  * @param unstr: XdmfUnstructuredGrid pointer.
  * @param timestep: Current timestep.
  * @param block: Current block number.
  *
  * @return: Created UnstructuredGrid pointer.
  */
vistle::UnstructuredGrid::ptr ReadSeisSol::generateUnstrGridFromXdmfGrid(XdmfUnstructuredGrid *unstr, const int block)
{
    if (!unstr)
        return nullptr;

    constexpr auto numCoords{3};

    const auto &xtopology = unstr->getTopology();
    const auto &xgeometry = unstr->getGeometry();
    const auto &xtopologyDimVec = xtopology->getDimensions();

    //get lowest val from topology dimension vector = corner
    const int &numCornerPerElem = *std::min_element(xtopologyDimVec.begin(), xtopologyDimVec.end());
    const auto xArrConn(XdmfArray::New());
    const auto xArrGeo(XdmfArray::New());
    const auto &geoContr = xgeometry->getHeavyDataController();
    const auto &topoContr = xtopology->getHeavyDataController();
    std::set<unsigned> uniqueVerts;

    if (m_parallelMode->getValue() == SERIAL) {
        readXdmfHeavyController(xArrConn.get(), topoContr);
        readXdmfHeavyController(xArrGeo.get(), geoContr);
    } else if (m_parallelMode->getValue() == BLOCKS) {
        readXdmfHDF5TopologyParallel(xArrConn.get(), topoContr, block);
        /* extractUniquePoints(xArrConn.get(), uniqueVerts); */

        //TODO: fetches at the moment all points => should fetch only relevant points => need unique verts
        // -> problem:
        //      - data fetched from topology is unordered => means: first element can contain last vertice
        //          => vistle tries to read last vert from coords array
        // -> solutions:
        //      1. create an vector with the size to store all geometry data, but read only relevant parameter.
        readXdmfHeavyController(xArrGeo.get(), geoContr);
    } else
        return nullptr;

    if (xArrGeo->getSize() == 0 || xArrConn->getSize() == 0)
        return nullptr;

    auto numVertices{xArrGeo->getSize() / numCoords};

    UnstructuredGrid::ptr unstr_ptr(
        new UnstructuredGrid(xArrConn->getSize() / numCornerPerElem, xArrConn->getSize(), numVertices));

    fillUnstrGridConnectList(unstr_ptr, xArrConn.get());
    fillUnstrGridCoords(unstr_ptr, xArrGeo.get(), uniqueVerts);
    fillUnstrGridElemList(unstr_ptr, numCornerPerElem);
    fillUnstrGridTypeList(unstr_ptr, XdmfUgridToVistleUgridType(unstr));

    //return copy of shared_ptr
    return unstr_ptr;
}

/**
  * @brief: TODO: implement format trancision between XdmfArrayType and vistle array type. []
  * Int8, Int16, Int32 = Integer
  * Float32 = float32
  * Float64 = double/Float
  * UInt8, UInt16, UInt32 = Unsigned
  */
void ReadSeisSol::setArrayType(boost::shared_ptr<const XdmfArrayType> type)
{
    if (type == XdmfArrayType::Int8()) {
        sendInfo("int8");
    } else if (type == XdmfArrayType::Int16()) {
        sendInfo("int16");
    } else if (type == XdmfArrayType::Int32()) {
        sendInfo("int32");
    } else if (type == XdmfArrayType::Int64()) {
        sendInfo("int64");
    } else if (type == XdmfArrayType::Float32()) {
        sendInfo("float32");
    } else if (type == XdmfArrayType::Float64()) {
        sendInfo("float64");
    } else if (type == XdmfArrayType::UInt8()) {
        sendInfo("UInt8");
    } else if (type == XdmfArrayType::UInt16()) {
        sendInfo("UInt16");
    } else if (type == XdmfArrayType::UInt32()) {
        sendInfo("UInt32");
    } else if (type == XdmfArrayType::String()) {
        sendInfo("String");
    }
}

/**
  * TODO: Implement trancision between XdmfAttributeCenter and vistle cell center. []
  * => relevant for color module not this reader.
  */
void ReadSeisSol::setGridCenter(const boost::shared_ptr<const XdmfAttributeCenter> &attCenter)
{
    if (attCenter == XdmfAttributeCenter::Grid())
        sendInfo("grid");
    else if (attCenter == XdmfAttributeCenter::Cell())
        sendInfo("cell");
    else if (attCenter == XdmfAttributeCenter::Node())
        sendInfo("node");
    else if (attCenter == XdmfAttributeCenter::Other())
        sendInfo("other");
}

/**
  * @brief: Initialize xdmf variables before read will be called.
  *
  * return: true if every xdmf variable could be initialized.
  */
bool ReadSeisSol::prepareReadXdmf()
{
    const std::string xfile = m_file->getValue();
    const auto xreader = XdmfReader::New();
    if (!xreader.get())
        return false;

    const auto xdomain = shared_dynamic_cast<XdmfDomain>(xreader->read(xfile));
    xgridCollect = xdomain->getGridCollection(0);
    const auto &xugrid = xgridCollect->getUnstructuredGrid(0);
    const auto &xattribute = xugrid->getAttribute(m_xAttSelect[m_attributes[0]->getValue()]);
    if (xattribute->getReadMode() == XdmfArray::Reference) {
        sendInfo("Hyperslab is not working in current version of this reader, the attributes have a "
                 "format that is not supported or "
                 "you aren't using HDF5 files. Chosen Param: %s",
                 m_attributes[0]->getValue().c_str());
        releaseXdmfObjects();
        return false;
    }
    return true;
}

/**
  * @brief: Call prepare function to chosen seismode.
  *
  * @return: true if everything went well in called prepare function.
  */
bool ReadSeisSol::prepareRead()
{
    return callSeisModeFunction(&ReadSeisSol::prepareReadXdmf, &ReadSeisSol::hdfModeNotImplemented);
}

/**
 * @brief Reuses grid from first timestep if m_reuseGrid is enabled or create a new one for each timestep otherwise.
 *
 * @param xugrid XdmfUnstructuredGrid pointer.
 * @param block current blockNum.
 *
 * @return Clone of m_unstr_grid if m_reuseGrid is true else the function creates a new one.
 */
vistle::UnstructuredGrid::ptr ReadSeisSol::reuseGrid(XdmfUnstructuredGrid *xugrid, int block)
{
    UnstructuredGrid::ptr ugrid_ptr;
    if (m_reuseGrid->getValue())
        ugrid_ptr = m_vugrid_ptr->clone();
    else {
        ugrid_ptr = generateUnstrGridFromXdmfGrid(xugrid, block);
        checkObjectPtr(this, ugrid_ptr, ERROR_GRID_GENERATION);
    }
    return ugrid_ptr;
}

/**
  * @brief: Called by read(). Reads internal data links to hdf5 files stored in xdmf.
  *
  * @token: Token for adding objects to ports.
  * @timestep: current timestep.
  * @block: current block.
  *
  * @return: true if everything has been added correctly to the ports.
  */
bool ReadSeisSol::readXdmfHDF5Data(Token &token, int timestep, int block)
{
    sendInfo("rank: %d timestep: %d block: %d", rank(), timestep, block);

    /*************** Create Unstructured Grid **************/
    const auto &xugrid = xgridCollect->getUnstructuredGrid(timestep);
    vistle::UnstructuredGrid::ptr vugrid_ptr = reuseGrid(xugrid.get(), block);

    /*************** Add ugrid to port & set meta data **************/
    if (m_gridOut->isConnected() && vugrid_ptr.get()) {
        vugrid_ptr->setBlock(block);
        vugrid_ptr->setTimestep(timestep);
        vugrid_ptr->updateInternals();
        token.addObject(m_gridOut, vugrid_ptr);
    }

    /*************** Create Scalars **************/
    for (size_t i = 0; i < NUM_ATTRIBUTES; ++i) {
        const auto &att_port = m_scalarsOut[i];
        if (!att_port->isConnected())
            continue;

        const auto &att_idx{m_xAttSelect[m_attributes[i]->getValue()]};
        const auto &xatt{xugrid->getAttribute(att_idx)};
        auto vatt_ptr = generateScalarFromXdmfAttribute(xatt.get(), timestep, block);
        if (checkObjectPtr(this, vatt_ptr, ERROR_ATT_GENERATION)) {
            /*************** Add attribute to port & set meta data **************/
            vatt_ptr->setGrid(vugrid_ptr);
            vatt_ptr->setBlock(block);
            vatt_ptr->setMapping(DataBase::Element);
            vatt_ptr->addAttribute("_species", xatt->getName());
            vatt_ptr->setTimestep(timestep);
            vatt_ptr->updateInternals();

            token.addObject(att_port, vatt_ptr);
        } else
            return false;
    }
    return true;
}

/**
  * @brief: Called every timestep and block.
  *
  * @token: Token for adding objects to ports.
  * @timestep: current timestep.
  * @block: current block.
  *
  * @return: true if everything has been added correctly to ports.
  */
bool ReadSeisSol::read(Token &token, int timestep, int block)
{
    if (timestep == -1) {
        if (m_reuseGrid->getValue()) {
            const auto &xugrid = xgridCollect->getUnstructuredGrid(0);
            m_vugrid_ptr = generateUnstrGridFromXdmfGrid(xugrid.get(), block);

            bool cOP{checkObjectPtr(this, m_vugrid_ptr, ERROR_GRID_GENERATION)};
            if (!cOP)
                releaseXdmfObjects();
            return cOP;
        }
        return true;
    }

    return callSeisModeFunction<bool, Token &, int, int>(
        &ReadSeisSol::readXdmfHDF5Data, &ReadSeisSol::hdfModeNotImplementedRead, token, timestep, block);
}

/**
  * @brief: Delete xdmf element variables cached.
  *
  * @return: no proper error handling at the moment.
  */
bool ReadSeisSol::finishReadXdmf()
{
    releaseXdmfObjects();
    sendInfo("Released reader objects rank: %d", rank());
    return true;
}

/**
  * @brief: Called after final read call. Calls function according chosen seismode.
  *
  * @return: Return value from called function.
  */
bool ReadSeisSol::finishRead()
{
    return callSeisModeFunction(&ReadSeisSol::finishReadXdmf, &ReadSeisSol::hdfModeNotImplemented);
}

/*************** DynamicPseudoEnum **************/
/**
  * @brief: Template for PseudoEnum which initializes map with an index like enum.
  *
  * @params: vector with values to store.
  */
template<class T>
void ReadSeisSol::DynamicPseudoEnum<T>::init(const std::vector<T> &params)
{
    std::for_each(params.begin(), params.end(), [n = 0, &indices = indices_map](auto &param) mutable {
        indices.insert({param, n++});
    });
}
