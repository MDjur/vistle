//header
#include "ReadTsunami.h"

//vistle
#include "vistle/core/database.h"
#include "vistle/core/index.h"
#include "vistle/core/layergrid.h"
#include "vistle/core/object.h"
#include "vistle/core/parameter.h"
#include "vistle/core/scalar.h"
#include "vistle/core/vec.h"
#include "vistle/module/module.h"

//vistle-module-util
#include "vistle/module/utils/ghost.h"

//mpi
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/intercommunicator.hpp>
#include <mpi.h>

//std
#include <algorithm>
#include <array>
#include <iterator>
#include <memory>
#include <omp.h>
#include <string>
#include <cstddef>
#include <vector>

using namespace vistle;
using namespace netCDF;

MODULE_MAIN_THREAD(ReadTsunami, boost::mpi::threading::multiple)
namespace {
constexpr auto lat{"lat"};
constexpr auto lon{"lon"};
constexpr auto bathy{"bathy"};
constexpr auto ETA{"eta"};
constexpr auto _species{"_species"};
constexpr auto fillValue{"fillValue"};
constexpr auto fillValueNew{"fillValueNew"};
constexpr auto NONE{"None"};
} // namespace

ReadTsunami::ReadTsunami(const string &name, int moduleID, mpi::communicator comm)
: vistle::Reader(name, moduleID, comm), m_needSea(false)
{
    // file-browser
    m_filedir = addStringParameter("file_dir", "NC File directory", "", Parameter::Filename);

    //ghost
    m_ghost = addIntParameter("ghost", "Show ghostcells.", 1, Parameter::Boolean);
    m_fill = addIntParameter("fill", "Replace filterValue.", 1, Parameter::Boolean);
    m_verticalScale = addFloatParameter("VerticalScale", "Vertical Scale parameter sea", 1.0);

    // define ports
    m_seaSurface_out = createOutputPort("Sea surface", "Grid Sea (Heightmap/LayerGrid)");
    m_groundSurface_out = createOutputPort("Ground surface", "Sea floor (Heightmap/LayerGrid)");

    // block size
    m_blocks[0] = addIntParameter("blocks latitude", "number of blocks in lat-direction", 2);
    m_blocks[1] = addIntParameter("blocks longitude", "number of blocks in lon-direction", 2);

    //fillvalue
    addFloatParameter(fillValue, "ncFile fillValue offset for eta", -9999.f);
    addFloatParameter(fillValueNew, "set new fillValue offset for eta", 0.0f);

    //bathymetryname
    m_bathy = addStringParameter("bathymetry ", "Select bathymetry stored in netCDF", "", Parameter::Choice);

    //scalar
    initScalarParamReader();

    // timestep built-in params
    setParameterRange(m_first, Integer(0), Integer(9999999));
    setParameterRange(m_last, Integer(-1), Integer(9999999));
    setParameterRange(m_increment, Integer(1), Integer(9999999));

    // observer these parameters
    observeParameter(m_filedir);
    observeParameter(m_blocks[0]);
    observeParameter(m_blocks[1]);
    observeParameter(m_verticalScale);

    setParallelizationMode(ParallelizeBlocks);
}

/**
 * @brief Destructor.
 */
ReadTsunami::~ReadTsunami()
{}

/**
 * @brief Initialize scalar choice parameter and ports.
 */
void ReadTsunami::initScalarParamReader()
{
    constexpr auto scalar_name{"Scalar "};
    constexpr auto scalarPort_name{"Scalar port "};
    constexpr auto scalarPort_descr{"Port for scalar number "};

    for (int i = 0; i < NUM_SCALARS; ++i) {
        const string i_str{to_string(i)};
        const string &scName = scalar_name + i_str;
        const string &portName = scalarPort_name + i_str;
        const string &portDescr = scalarPort_descr + i_str;
        m_scalars[i] = addStringParameter(scName, "Select scalar.", "", Parameter::Choice);
        m_scalarsOut[i] = createOutputPort(portName, portDescr);
        observeParameter(m_scalars[i]);
    }
}

/**
 * @brief Print string only on rank 0.
 *
 * @param str Format str to print.
 */
template<class... Args>
void ReadTsunami::printRank0(const string &str, Args... args) const
{
    if (rank() == 0)
        sendInfo(str, args...);
}

/**
  * @brief Prints current rank and the number of all ranks to the console.
  */
void ReadTsunami::printMPIStats() const
{
    printRank0("Current Rank: " + to_string(rank()) + " Processes (MPISIZE): " + to_string(size()));
}

/**
  * @brief Called when any of the reader parameter changing.
  *
  * @param: Parameter that got changed.
  * @return: true if all essential parameters could be initialized.
  */
bool ReadTsunami::examine(const vistle::Parameter *param)
{
    if (!param || param == m_filedir) {
        printMPIStats();

        if (!inspectNetCDF())
            return false;
    }

    const int &nBlocks = m_blocks[0]->getValue() * m_blocks[1]->getValue();
    if (nBlocks % size() != 0) {
        printRank0("Number of blocks = x * MPISIZE!");
        return false;
    }
    setPartitions(nBlocks);
    return true;
}


/**
 * @brief Open NcmpiFile and catch failure.
 *
 * @return Return open ncFile as unique_ptr if there is no error, else a nullptr.
 */
/* unique_ptr<NcmpiFile> ReadTsunami::openNcmpiFile() */
unique_ptr<NcFile> ReadTsunami::openNcmpiFile()
{
    const auto &fileName = m_filedir->getValue();
    try {
        return make_unique<NcFile>(fileName, NcFile::read);
    } catch (...) {
        sendInfo("Couldn't open NetCDF file!");
        return nullptr;
    }
}

/**
 * @brief Inspect netCDF dimension (all independent variables.)
 *
 * @param ncFile Open nc file pointer.
 *
 * @return tru if everything is initialized.
 */
bool ReadTsunami::inspectDims(const NcFile *ncFile)
{
    const int &maxTime = ncFile->getDim("time").getSize();
    setTimesteps(maxTime);

    const Integer &maxlatDim = ncFile->getDim(lat).getSize();
    const Integer &maxlonDim = ncFile->getDim(lon).getSize();
    setParameterRange(m_blocks[0], Integer(1), maxlatDim);
    setParameterRange(m_blocks[1], Integer(1), maxlonDim);
    return inspectScalars(ncFile);
}

/**
 * @brief Helper for checking if name includes contains.
 *
 * @param name name of attribute.
 * @param contains string which will be checked for.
 *
 * @return true if contains is in name.
 */
inline auto strContains(const string &name, const string &contains)
{
    return name.find(contains) != string::npos;
}

/**
 * @brief Helper which adds NONE to vector if its empty.
 *
 * @param vec ref to string vector.
 */
inline void ifEmptyAddNone(vector<string> &vec)
{
    if (vec.empty())
        vec.push_back(NONE);
}

/**
 * @brief Inspect netCDF variables which depends on dim variables.
 *
 * @param ncFile Open nc file pointer.
 *
 * @return true if everything is initialized.
 */
bool ReadTsunami::inspectScalars(const NcFile *ncFile)
{
    vector<string> scalarChoiceVec;
    vector<string> bathyChoiceVec;
    auto latLonContainsGrid = [&](auto &name, int i) {
        if (strContains(name, "grid"))
            m_latLon_Ground[i] = name;
        else
            m_latLon_Sea[i] = name;
    };

    //delete previous choicelists.
    m_bathy->setChoices(vector<string>());
    for (const auto &scalar: m_scalars)
        scalar->setChoices(vector<string>());

    //read names of scalars
    for (auto &[name, val]: ncFile->getVars()) {
        if (strContains(name, lat))
            latLonContainsGrid(name, 0);
        else if (strContains(name, lon))
            latLonContainsGrid(name, 1);
        else if (strContains(name, bathy)) // || strContains(name, "deformation"))
            bathyChoiceVec.push_back(name);
        else if (val.getDimCount() == 2) // for now: only scalars with 2 dim depend on lat lon.
            scalarChoiceVec.push_back(name);
    }
    ifEmptyAddNone(scalarChoiceVec);
    ifEmptyAddNone(bathyChoiceVec);

    //init choice param with scalardata
    setParameterChoices(m_bathy, bathyChoiceVec);

    for (auto &scalar: m_scalars)
        setParameterChoices(scalar, scalarChoiceVec);

    if (m_latLon_Ground.empty() || m_latLon_Sea.empty()) {
        sendInfo("No parameter lat, lon, grid_lat or grid_lon found. Reader not able to read tsunami.");
        return false;
    }
    return true;
}


/**
 * @brief Inspect netCDF variables stored in file.
 */
bool ReadTsunami::inspectNetCDF()
{
    auto ncFile = openNcmpiFile();

    if (!ncFile)
        return false;

    return inspectDims(ncFile.get());
}

/**
 * @brief Fill z coordinate from layergrid.
 *
 * @surface: LayerGrid pointer which will be modified.
 * @dim: Dimension in lat (y) and lon (x).
 * @zCalc: Function for computing z-coordinate.
 */
template<class T>
void ReadTsunami::fillHeight(LayerGrid::ptr surface, const Dim<T> &dim, const ZCalcFunc &zCalc)
{
    int n = 0;
    auto z = surface->z().begin();
    for (T i = 0; i < dim.X(); ++i)
        for (T j = 0; j < dim.Y(); ++j, ++n)
            z[n] = zCalc(i, j);
}

/**
  * @brief Generates NcVarParams struct which contains start, count and stride values computed based on given parameters.
  *
  * @dim: Current dimension of NcVar.
  * @ghost: Number of ghost cells to add.
  * @numDimBlocks: Total number of blocks for this direction.
  * @partition: Partition scalar.
  * @return: NcVarParams object.
  */
template<class T, class PartionIdx>
auto ReadTsunami::generateNcVarExt(const NcVar &ncVar, const T &dim, const T &ghost, const T &numDimBlock,
                                   const PartionIdx &partition) const
{
    T count = dim / numDimBlock;
    T start = partition * count;
    structured_ghost_addition(start, count, dim, ghost);
    return PNcVarExt(ncVar, start, count);
}

/**
  * @brief Called for each timestep and for each block (MPISIZE).
  *
  * @token: Ref to internal vistle token.
  * @timestep: current timestep.
  * @block: current block number of parallelization.
  * @return: true if all data is set and valid.
  */
bool ReadTsunami::read(Token &token, int timestep, int block)
{
    return computeBlock(token, block, timestep);
}

/**
 * @brief Called before read and used to initialize steering parameters.
 *
 * @return True if everything is initialized.
 */
bool ReadTsunami::prepareRead()
{
    const auto &part = numPartitions();

    m_block_etaIdx = vector<int>(part);
    m_block_dimSea = vector<moffDim>(part);
    m_block_etaVec = vector<vector<float>>(part);
    m_block_VecScalarPtr = VecArrVecScalarPtrs(part);
    m_block_min = VecLatLon(part);
    m_block_max = VecLatLon(part);

    m_needSea = m_seaSurface_out->isConnected();
    for (auto &val: m_scalarsOut)
        m_needSea = m_needSea || val->isConnected();
    return true;
}

/**
  * @brief Computing per block.
  *
  * @token: Ref to internal vistle token.
  * @block: current block number of parallelization.
  * @timestep: current timestep.
  * @return: true if all data is set and valid.
  */
bool ReadTsunami::computeBlock(Reader::Token &token, const int block, const int timestep)
{
    if (timestep == -1)
        return computeConst(token, block);
    else if (m_needSea)
        return computeTimestep(token, block, timestep);
    return true;
}

/**
 * @brief Compute the unique block partition index for current block.
 *
 * @tparam Iter Iterator.
 * @param block Current block.
 * @param nLatBlocks Storage for number of blocks for latitude.
 * @param nLonBlocks Storage for number of blocks for longitude.
 * @param blockPartitionIterFirst Start iterator for storage partition indices.
 */
template<class Iter>
void ReadTsunami::computeBlockPartition(const int block, vistle::Index &nLatBlocks, vistle::Index &nLonBlocks,
                                        Iter blockPartitionIterFirst)
{
    array<Index, NUM_BLOCKS> blocks;
    for (int i = 0; i < NUM_BLOCKS; ++i)
        blocks[i] = m_blocks[i]->getValue();

    nLatBlocks = blocks[0];
    nLonBlocks = blocks[1];

    structured_block_partition(blocks.begin(), blocks.end(), blockPartitionIterFirst, block);
}

/**
 * @brief Initialize eta values (wave amplitude per timestep).
 *
 * @param ncFile open ncfile pointer.
 * @param ncExtSea array which contains latitude and longitude NcVarExtended objects.
 * @param time ReadTimer object for current run.
 * @param verticesSea number of vertices for sea surface of the current block.
 * @param block current block.
 */
void ReadTsunami::initETA(const NcFile *ncFile, const array<PNcVarExt, 2> &ncExtSea, const ReaderTime &time,
                          const size_t &verticesSea, int block)
{
    const NcVar &etaVar = ncFile->getVar(ETA);
    const auto &latSea = ncExtSea[0];
    const auto &lonSea = ncExtSea[1];
    const auto &nTimesteps = time.calc_numtime();
    m_block_etaVec[block] = vector<float>(nTimesteps * verticesSea);
    const vector<size_t> vecStartEta{static_cast<size_t>(time.first()), latSea.Start(), lonSea.Start()};
    const vector<size_t> vecCountEta{static_cast<size_t>(nTimesteps), latSea.Count(), lonSea.Count()};
    const vector<ptrdiff_t> vecStrideEta{static_cast<ptrdiff_t>(time.inc()), latSea.Stride(), lonSea.Stride()};
    auto &vecEta = m_block_etaVec[block];
    etaVar.getVar(vecStartEta, vecCountEta, vecStrideEta, vecEta.data());

    //filter fillvalue
    if (m_fill->getValue()) {
        const float &fillNew = getFloatParameter(fillValueNew);
        const float &fill = getFloatParameter(fillValue);
        //TODO: Bad! needs rework.
        replace(vecEta.begin(), vecEta.end(), fill, fillNew);
    }
}

/**
 * @brief Initialize scalars for each attribute and each block.
 *
 * @param ncFile open ncfile pointer.
 * @param ncExtSea array which contains latitude and longitude NcVarExtended objects.
 * @param verticesSea number of vertices for the sea surface.
 * @param block current block.
 */
void ReadTsunami::initScalars(const NcFile *ncFile, const array<PNcVarExt, 2> &ncExtSea, const size_t &verticesSea,
                              int block)
{
    const auto &latSea = ncExtSea[0];
    const auto &lonSea = ncExtSea[1];
    auto &vecScalarPtrArr = m_block_VecScalarPtr[block];
    for (int i = 0; i < NUM_SCALARS; ++i) {
        if (!m_scalarsOut[i]->isConnected())
            continue;
        const auto &scName = m_scalars[i]->getValue();
        if (scName == NONE)
            continue;
        vector<float> vecScalar(verticesSea);
        const vector<size_t> vecScalarStart{latSea.Start(), lonSea.Start()};
        const vector<size_t> vecScalarCount{latSea.Count(), lonSea.Count()};
        const vector<ptrdiff_t> vecScalarStride{1, 1};
        const vector<ptrdiff_t> vecScalarImap{1, static_cast<ptrdiff_t>(latSea.Count())}; // transponse
        const auto &val = ncFile->getVar(scName);
        VisVecScalarPtr scalarPtr(new Vec<Scalar>(verticesSea));
        vecScalarPtrArr[i] = scalarPtr;
        val.getVar(vecScalarStart, vecScalarCount, vecScalarStride, vecScalarImap, scalarPtr->x().begin());

        //set some meta data
        scalarPtr->addAttribute(_species, scName);
        scalarPtr->setBlock(block);
    }
}

/**
 * @brief Initialize sea surface variables for current block.
 *
 * @param ncFile open ncfile pointer.
 * @param nBlocks number of blocks in latitude and longitude direction.
 * @param nBlockPartIdx partition index for latitude and longitude.
 * @param time ReaderTime object which holds time parameters.
 * @param ghost ghostcells to add.
 * @param block current block.
 */
void ReadTsunami::initSea(const NcFile *ncFile, const array<Index, 2> &nBlocks,
                          const array<Index, NUM_BLOCKS> &nBlockPartIdx, const ReaderTime &time, int ghost, int block)
{
    const NcVar &latVar = ncFile->getVar(m_latLon_Sea[0]);
    const NcVar &lonVar = ncFile->getVar(m_latLon_Sea[1]);
    const sztDim dimSeaTotal(latVar.getDim(0).getSize(), lonVar.getDim(0).getSize(), 1);

    array<PNcVarExt, 2> latLonSea;
    latLonSea[0] = generateNcVarExt<MPI_Offset, Index>(latVar, dimSeaTotal.X(), ghost, nBlocks[0], nBlockPartIdx[0]);
    latLonSea[1] = generateNcVarExt<MPI_Offset, Index>(lonVar, dimSeaTotal.Y(), ghost, nBlocks[1], nBlockPartIdx[1]);

    // lat = Y and lon = X
    const auto &latSea = latLonSea[0];
    const auto &lonSea = latLonSea[1];

    const size_t &vertSea = latSea.Count() * lonSea.Count();
    initETA(ncFile, latLonSea, time, vertSea, block);

    vector<float> vecLat(latSea.Count()), vecLon(lonSea.Count());

    latSea.readNcVar(vecLat.data());
    lonSea.readNcVar(vecLon.data());
    m_block_dimSea[block] = moffDim(lonSea.Count(), latSea.Count(), 1);
    m_block_min[block][0] = *vecLon.begin();
    m_block_min[block][1] = *vecLat.begin();
    m_block_max[block][0] = *(vecLon.end() - 1);
    m_block_max[block][1] = *(vecLat.end() - 1);

    initScalars(ncFile, latLonSea, vertSea, block);
}

/**
 * @brief Initialize and create ground surface.
 *
 * @param token Token object ref.
 * @param ncFile open ncfile pointer.
 * @param nBlocks number of blocks in latitude and longitude.
 * @param blockPartIdx partition index for latitude and longitude.
 * @param ghost ghostcells to add.
 * @param block current block.
 */
void ReadTsunami::createGround(Token &token, const NcFile *ncFile, const array<vistle::Index, 2> &nBlocks,
                               const array<vistle::Index, NUM_BLOCKS> &blockPartIdx, int ghost, int block)
{
    const NcVar &gridLatVar = ncFile->getVar(m_latLon_Ground[0]);
    const NcVar &gridLonVar = ncFile->getVar(m_latLon_Ground[1]);

    const sztDim dimGroundTotal(gridLatVar.getDim(0).getSize(), gridLonVar.getDim(0).getSize(), 0);
    const auto latGround =
        generateNcVarExt<MPI_Offset, Index>(gridLatVar, dimGroundTotal.X(), ghost, nBlocks[0], blockPartIdx[0]);
    const auto lonGround =
        generateNcVarExt<MPI_Offset, Index>(gridLonVar, dimGroundTotal.Y(), ghost, nBlocks[1], blockPartIdx[1]);
    const size_t &verticesGround = latGround.Count() * lonGround.Count();
    vector<float> vecDepth(verticesGround);
    const NcVar &bathymetryVar = ncFile->getVar(m_bathy->getValue());
    const vector<size_t> vecStartBathy{latGround.Start(), lonGround.Start()};
    const vector<size_t> vecCountBathy{latGround.Count(), lonGround.Count()};
    bathymetryVar.getVar(vecStartBathy, vecCountBathy, vecDepth.data());
    if (!vecDepth.empty()) {
        vector<float> vecLatGrid(latGround.Count()), vecLonGrid(lonGround.Count());
        latGround.readNcVar(vecLatGrid.data());
        lonGround.readNcVar(vecLonGrid.data());

        const auto &scale = m_verticalScale->getValue();
        LayerGrid::ptr grndPtr(new LayerGrid(lonGround.Count(), latGround.Count(), 1));
        const auto grndDim = moffDim(lonGround.Count(), latGround.Count(), 0);
        grndPtr->min()[0] = *vecLonGrid.begin();
        grndPtr->min()[1] = *vecLatGrid.begin();
        grndPtr->max()[0] = *(vecLonGrid.end() - 1);
        grndPtr->max()[1] = *(vecLatGrid.end() - 1);
        fillHeight(grndPtr, grndDim, [&vecDepth, &lonGround, &scale](const auto &j, const auto &k) {
            return -vecDepth[k * lonGround.Count() + j] * scale;
        });

        // add ground data to port
        grndPtr->setBlock(block);
        grndPtr->setTimestep(-1);
        grndPtr->updateInternals();
        token.addObject(m_groundSurface_out, grndPtr);
    }
}

/**
  * @brief Generates the initial surfaces for sea and ground and adds only ground to scene.
  *
  * @token: Ref to internal vistle token.
  * @block: current block number of parallel process.
  * @return: true if all data could be initialized.
  */
bool ReadTsunami::computeConst(Token &token, const int block)
{
    m_mtx.lock();
    auto ncFile = openNcmpiFile();
    m_mtx.unlock();
    if (!ncFile)
        return false;

    array<Index, 2> nBlocks;
    array<Index, NUM_BLOCKS> blockPartitionIdx;
    computeBlockPartition(block, nBlocks[0], nBlocks[1], blockPartitionIdx.begin());

    int ghost{0};
    if (m_ghost->getValue() == 1 && !(nBlocks[0] == 1 && nBlocks[1] == 1))
        ghost++;

    if (m_needSea) {
        auto last = m_last->getValue();
        if (last < 0)
            return false;
        auto inc = m_increment->getValue();
        auto first = m_first->getValue();
        last = last - (last % inc);
        const auto &time = ReaderTime(first, last, inc);
        initSea(ncFile.get(), nBlocks, blockPartitionIdx, time, ghost, block);
    }

    if (m_groundSurface_out->isConnected()) {
        if (m_bathy->getValue() == NONE)
            printRank0("File doesn't provide bathymetry data");
        else
            createGround(token, ncFile.get(), nBlocks, blockPartitionIdx, ghost, block);
    }
    return true;
}

/**
  * @brief Generates heightmap (layergrid) for corresponding timestep and adds object to scene.
  *
  * @token: Ref to internal vistle token.
  * @block: current block number of parallel process.
  * @timestep: current timestep.
  * @return: true. TODO: add proper error-handling here.
  */
bool ReadTsunami::computeTimestep(Token &token, const int block, const int timestep)
{
    auto &blockSeaDim = m_block_dimSea[block];
    LayerGrid::ptr gridPtr(new LayerGrid(blockSeaDim.X(), blockSeaDim.Y(), blockSeaDim.Z()));
    copy_n(m_block_min[block].begin(), 2, gridPtr->min());
    copy_n(m_block_max[block].begin(), 2, gridPtr->max());

    // getting z from vecEta and copy to z()
    // verticesSea * timesteps = total Count() vecEta
    const auto &verticesSea = gridPtr->getNumVertices();
    vector<float> &etaVec = m_block_etaVec[block];
    auto count = m_block_etaIdx[block]++ * verticesSea;
    fillHeight(gridPtr, blockSeaDim, [&etaVec, &blockSeaDim, &count](const auto &j, const auto &k) {
        return etaVec[count + k * blockSeaDim.X() + j];
    });
    gridPtr->updateInternals();
    gridPtr->setTimestep(timestep);
    gridPtr->setBlock(block);
    if (m_seaSurface_out->isConnected())
        token.addObject(m_seaSurface_out, gridPtr);

    //add scalar to ports
    ArrVecScalarPtrs &arrVecScaPtrs = m_block_VecScalarPtr[block];
    for (int i = 0; i < NUM_SCALARS; ++i) {
        if (!m_scalarsOut[i]->isConnected())
            continue;

        auto vecScalarPtr = arrVecScaPtrs[i]->clone();
        vecScalarPtr->setGrid(gridPtr);
        vecScalarPtr->addAttribute(_species, vecScalarPtr->getAttribute(_species));
        vecScalarPtr->setBlock(block);
        vecScalarPtr->setTimestep(timestep);
        vecScalarPtr->updateInternals();

        token.addObject(m_scalarsOut[i], vecScalarPtr);
    }

    return true;
}

/**
 * @brief Called after last read operation for clearing reader parameter.
 *
 * @return true
 */
bool ReadTsunami::finishRead()
{
    //reset block partition variables
    m_block_VecScalarPtr.clear();
    m_block_dimSea.clear();
    m_block_etaVec.clear();
    m_block_max.clear();
    m_block_min.clear();
    m_block_etaIdx.clear();
    sendInfo("Cleared Cache for rank: " + to_string(rank()));
    return true;
}
