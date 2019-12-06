#ifndef VISIT_VISTLE_ENGINE_H
#define VISIT_VISTLE_ENGINE_H

#include <mpi.h>
#include <boost/asio.hpp>

#include "ModuleInterface.h"
#include "MetaData.h"

#include <string>
#include <vector>
#include <memory>
#include <functional>
#include <mutex>
#include <thread>

#include <manager/manager.h>
#include <module/module.h>


namespace in_situ {
enum class SimulationDataTyp {
     mesh
    ,variable
    ,material
    ,curve
    ,expression
    ,species
    ,genericCommand
    ,customCommand
    ,message

};

class V_VISITXPORT Engine {
public:
    static Engine* createEngine();
    static void DisconnectSimulation();
    bool initialize(int argC, char** argV);
    bool isInitialized() const noexcept;
    bool setMpiComm(void* newConn);


    //********************************
    //***functions called by module***
    //********************************
    void setModule(vistle::Module* module);
    void setDoReadMutex(std::mutex* m);
    bool getNumObjects(SimulationDataTyp type, int& num);
    bool getNthObject(SimulationDataTyp type, int n, visit_handle&);
    bool getDataNames(SimulationDataTyp type, std::vector<std::string>& names);
    //set callbacks (called from module)
    void SetTimestepChangedCb(std::function<bool(void)> cb);
    void SetDisconnectCb(std::function<void(void)> cb);
    //********************************
    //****functions called by sim****
    //********************************
    //adds all available data to the according outputs to execute the pipeline
    bool sendData();
    void SimulationTimeStepChanged();
    bool getMetaData(Metadata& md);
    void SimulationInitiateCommand(const char* command);
    void DeleteData();
    //set callbacks (called from sim)
    void SetSimulationCommandCallback(void(*sc)(const char*, const char*, void*), void* scdata);


private:
    static Engine* instance;
    bool m_initialized = false;
    MPI_Comm comm = MPI_COMM_WORLD;
    vistle::Module* m_module = nullptr;
    std::mutex* m_doReadMutex = nullptr;
    std::thread managerThread;

    std::map<std::string, vistle::Port*> m_portsList;
    std::vector<vistle::Object::ptr> m_dataObjects;
    //callbacks from ConnectLibSim module
    std::function<bool(void)> timestepChangedCb; //returns true, if module is ready to receive data;
    std::function<void(void)> disconnectCb;


    //callbacks from simulation
    void (*simulationCommandCallback)(const char*, const char*, void*) = nullptr;
    void* simulationCommandCallbackData = nullptr;


    bool addPorts();

    bool makeCurvilinearMesh(visit_handle h);
    bool makeUntructuredMesh(visit_handle h);
    bool makeAmrMesh(visit_handle h);
    bool sendDataToModule();



    Engine();
    ~Engine();


    void printToConsole(const std::string& msg) const;
};





}

#endif // !VISIT_VISTLE_ENGINE_H
