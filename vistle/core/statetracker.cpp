#include <boost/foreach.hpp>

#include "message.h"
#include "messages.h"
#include "messagerouter.h"
#include "parameter.h"
#include "port.h"
#include "porttracker.h"
#include <cassert>

#include "statetracker.h"

#include <vistle/util/vecstreambuf.h>
#include "archives.h"

#define CERR \
   std::cerr << m_name << ": "

//#define DEBUG

namespace bi = boost::interprocess;

namespace vistle {

namespace {
const std::string unknown("(unknown)");
}

using message::Id;

int StateTracker::Module::state() const {

   int s = StateObserver::Known;
   if (initialized)
      s |= StateObserver::Initialized;
   if (busy)
      s |= StateObserver::Busy;
   if (killed)
      s |= StateObserver::Killed;
   return s;
}

StateTracker::StateTracker(const std::string &name, std::shared_ptr<PortTracker> portTracker)
: m_portTracker(portTracker)
, m_traceType(message::INVALID)
, m_traceId(Id::Invalid)
, m_name(name)
{
   if (!m_portTracker) {
      m_portTracker.reset(new PortTracker());
      m_portTracker->setTracker(this);
   }

   Module session(Id::Vistle, Id::MasterHub); // for session parameters
   session.name = "Vistle Session";
   runningMap.emplace(Id::Vistle, session);
}

StateTracker::mutex &StateTracker::getMutex() {

   return m_replyMutex;
}

int StateTracker::getMasterHub() const {

   return Id::MasterHub;
}

std::vector<int> StateTracker::getHubs() const {
    mutex_locker guard(m_stateMutex);
    std::vector<int> hubs;
    for (const auto &h: m_hubs)
       hubs.push_back(h.id);
    return hubs;
}

std::vector<int> StateTracker::getSlaveHubs() const {
    mutex_locker guard(m_stateMutex);
    std::vector<int> hubs;
    for (const auto &h: m_hubs)
       if (h.id != Id::MasterHub)
          hubs.push_back(h.id);
    return hubs;
}

const std::string &StateTracker::hubName(int id) const {

    mutex_locker guard(m_stateMutex);
    for (const auto &h: m_hubs) {
       if (h.id == id)
          return h.name;
    }
    return unknown;
}

std::vector<int> StateTracker::getRunningList() const {

   mutex_locker guard(m_stateMutex);
   std::vector<int> result;
   for (RunningMap::const_iterator it = runningMap.begin();
         it != runningMap.end(); 
         ++it) {
      if (Id::isModule(it->first))
          result.push_back(it->first);
   }
   return result;
}

std::vector<int> StateTracker::getBusyList() const {

   mutex_locker guard(m_stateMutex);
   std::vector<int> result;
   for (ModuleSet::const_iterator it = busySet.begin();
         it != busySet.end();
         ++it) {
      result.push_back(*it);
   }
   return result;
}

int StateTracker::getHub(int id) const {

   mutex_locker guard(m_stateMutex);
   if (Id::isHub(id)) {
       return id;
   }

   RunningMap::const_iterator it = runningMap.find(id);
   if (it == runningMap.end()) {
      it = quitMap.find(id);
      if (it == quitMap.end())
         return Id::Invalid;
   }

   if (it->second.hub > Id::MasterHub) {
      CERR << "getHub for " << id << " failed - invalid value " << it->second.hub  << std::endl;
   }
   //assert(it->second.hub <= Id::MasterHub);
   return it->second.hub;
}

const HubData &StateTracker::getHubData(int id) const {
    static HubData invalidHub(Id::Invalid, "");

    mutex_locker guard(m_stateMutex);
    for (const auto &hub: m_hubs) {
        if (hub.id == id)
            return hub;
    }

    return invalidHub;
}

std::string StateTracker::getModuleName(int id) const {

   mutex_locker guard(m_stateMutex);
   RunningMap::const_iterator it = runningMap.find(id);
   if (it == runningMap.end())
      return std::string();

   return it->second.name;
}

int StateTracker::getModuleState(int id) const {

   mutex_locker guard(m_stateMutex);
   RunningMap::const_iterator it = runningMap.find(id);
   if (it == runningMap.end()) {
      it = quitMap.find(id);
      if (it == quitMap.end())
         return StateObserver::Unknown;
      else
         return StateObserver::Quit;
   }

   return it->second.state();
}

namespace {

void appendMessage(std::vector<StateTracker::MessageWithPayload> &v, const message::Message &msg, std::shared_ptr<const buffer> payload = std::shared_ptr<const buffer>()) {

   v.emplace_back(msg, payload);
}

}

StateTracker::VistleState StateTracker::getState() const {

   mutex_locker guard(m_stateMutex);
   using namespace vistle::message;
   VistleState state;

   for (const auto &slave: m_hubs) {
      AddHub msg(slave.id, slave.name);
      msg.setNumRanks(slave.numRanks);
      msg.setPort(slave.port);
      msg.setDataPort(slave.dataPort);
      if (!slave.address.is_unspecified())
          msg.setAddress(slave.address);
      appendMessage(state, msg);
   }

   // available modules
   auto avail = availableModules();
   for(const auto &keymod: avail) {
      const auto &mod = keymod.second;
      appendMessage(state, ModuleAvailable(mod.hub, mod.name, mod.path));
   }

   // loaded map
   appendMessage(state, UpdateStatus(UpdateStatus::LoadedFile, m_loadedWorkflowFile));

   // modules with parameters and ports
   for (auto &it: runningMap) {
      const int id = it.first;
      const Module &m = it.second;

      if (Id::isModule(id)) {
          Spawn spawn(m.hub, m.name);
          spawn.setSpawnId(id);
          appendMessage(state, spawn);

          if (m.initialized) {
              Started s(m.name);
              s.setSenderId(id);
              s.setReferrer(spawn.uuid());
              appendMessage(state, s);
          }

          if (m.busy) {
              Busy b;
              b.setSenderId(id);
              appendMessage(state, b);
          }

          if (m.killed) {
              Kill k(id);
              appendMessage(state, k);
          }
      }

      const ParameterMap &pmap = m.parameters;
      for (const auto &it2: m.paramOrder) {
         //CERR << "module " << id << ": " << it2.first << " -> " << it2.second << std::endl;
         const std::string &name = it2.second;
         const auto it3 = pmap.find(name);
         assert(it3 != pmap.end());
         const auto param = it3->second;

         AddParameter add(*param, getModuleName(id));
         add.setSenderId(id);
         appendMessage(state, add);

         SetParameter setDef(id, name, param, Parameter::Value, true);
         setDef.setSenderId(id);
         appendMessage(state, setDef);

         if (param->presentation() == Parameter::Choice) {
            SetParameterChoices choices(name, param->choices().size());
            choices.setSenderId(id);
            SetParameterChoices::Payload pl(param->choices());
            auto vec = addPayload(choices, pl);
            auto shvec = std::make_shared<buffer>(vec);
            appendMessage(state, choices, shvec);
         }

         SetParameter setV(id, name, param, Parameter::Value);
		 setV.setSenderId(id);
		 appendMessage(state, setV);
         SetParameter setMin(id, name, param, Parameter::Minimum);
		 setMin.setSenderId(id);
		 appendMessage(state, setMin);
         SetParameter setMax(id, name, param, Parameter::Maximum);
		 setMax.setSenderId(id);
		 appendMessage(state, setMax);
      }

      if (portTracker()) {
         for (auto &portname: portTracker()->getInputPortNames(id)) {
            AddPort cp(*portTracker()->getPort(id, portname));
            cp.setSenderId(id);
            appendMessage(state, cp);
         }

         for (auto &portname: portTracker()->getOutputPortNames(id)) {
            AddPort cp(*portTracker()->getPort(id, portname));
            cp.setSenderId(id);
            appendMessage(state, cp);
         }
      }
   }

   // connections
   for (auto &it: runningMap) {
      const int id = it.first;

      if (portTracker()) {
         for (auto &portname: portTracker()->getOutputPortNames(id)) {
            const Port::ConstPortSet *connected = portTracker()->getConnectionList(id, portname);
            for (auto &dest: *connected) {
               Connect c(id, portname, dest->getModuleID(), dest->getName());
               appendMessage(state, c);
            }
         }

         for (auto &paramname: getParameters(id)) {
            const Port::ConstPortSet *connected = portTracker()->getConnectionList(id, paramname);
            for (auto &dest: *connected) {
               Connect c(id, paramname, dest->getModuleID(), dest->getName());
               appendMessage(state, c);
            }
         }
      }
   }

   for (const auto &m: m_queue)
      appendMessage(state, m.message, m.payload);

   // finalize
   appendMessage(state, ReplayFinished());

   return state;
}

const std::map<AvailableModule::Key, AvailableModule> &StateTracker::availableModules() const {

    mutex_locker guard(m_stateMutex);
    return m_availableModules;
}

bool StateTracker::handle(const message::Message &msg, const buffer *payload, bool track) {
    return handle(msg, payload?payload->data():nullptr, payload?payload->size():0, track);
}

bool StateTracker::handle(const message::Message &msg, const char *payload, size_t payloadSize, bool track) {

   using namespace vistle::message;

   mutex_locker guard(m_stateMutex);

   ++m_numMessages;
   m_aggregatedPayload += msg.payloadSize();

#ifndef NDEBUG
   if (msg.type() != message::ADDOBJECT && msg.uuid() != msg.referrer()) {
       if (m_alreadySeen.find(msg.uuid()) != m_alreadySeen.end()) {
           CERR << "duplicate message: " << msg << std::endl;
       }
       m_alreadySeen.insert(msg.uuid());
   }
#endif

   if (m_traceId != Id::Invalid && m_traceType != INVALID) {

      if (msg.type() == m_traceType || m_traceType == ANY) {

         if (msg.senderId() == m_traceId || msg.destId() == m_traceId || m_traceId == Id::Broadcast) {
            std::cout << m_name << ": " << msg << std::endl << std::flush;
         }
      }
   }

   if (!track)
      return true;

   bool handled = true;

   buffer pl;
   if (payload) {
       std::copy(payload, payload+payloadSize, std::back_inserter(pl));
   }

   mutex_locker locker(getMutex());
   switch (msg.type()) {
      case IDENTIFY: {
         break;
      }
      case ADDHUB: {
         const AddHub &slave = static_cast<const AddHub &>(msg);
         handled = handlePriv(slave);
         break;
      }
      case SPAWN: {
         const Spawn &spawn = static_cast<const Spawn &>(msg);
         registerReply(msg.uuid(), msg);
         handled = handlePriv(spawn);
         break;
      }
      case SPAWNPREPARED: {
         break;
      }
      case STARTED: {
         const Started &started = static_cast<const Started &>(msg);
         handled = handlePriv(started);
         break;
      }
      case KILL: {
         const Kill &kill = static_cast<const Kill &>(msg);
         handled = handlePriv(kill);
         break;
      }
      case DEBUG: {
         break;
      }
      case QUIT: {
         const Quit &quit = static_cast<const Quit &>(msg);
         handled = handlePriv(quit);
         break;
      }
      case MODULEEXIT: {
         const ModuleExit &modexit = static_cast<const ModuleExit &>(msg);
         handled = handlePriv(modexit);
         break;
      }
      case EXECUTE: {
         break;
      }
      case CANCELEXECUTE: {
         break;
      }
      case ADDOBJECT: {
         break;
      }
      case ADDOBJECTCOMPLETED: {
         break;
      }
      case ADDPORT: {
         const AddPort &cp = static_cast<const AddPort &>(msg);
         handled = handlePriv(cp);
         break;
      }
      case REMOVEPORT: {
         const RemovePort &dp = static_cast<const RemovePort &>(msg);
         handled = handlePriv(dp);
         break;
      }
      case ADDPARAMETER: {
         const AddParameter &add = static_cast<const AddParameter &>(msg);
         handled = handlePriv(add);
         break;
      }
      case REMOVEPARAMETER: {
         const RemoveParameter &rem = static_cast<const RemoveParameter &>(msg);
         handled = handlePriv(rem);
         break;
      }
      case CONNECT: {
         const Connect &conn = static_cast<const Connect &>(msg);
         handled = handlePriv(conn);
         break;
      }
      case DISCONNECT: {
         const Disconnect &disc = static_cast<const Disconnect &>(msg);
         handled = handlePriv(disc);
         break;
      }
      case SETPARAMETER: {
         const SetParameter &set = static_cast<const SetParameter &>(msg);
         handled = handlePriv(set);
         break;
      }
      case SETPARAMETERCHOICES: {
         const SetParameterChoices &choice = static_cast<const SetParameterChoices &>(msg);
         handled = handlePriv(choice, pl);
         break;
      }
      case PING: {
         const Ping &ping = static_cast<const Ping &>(msg);
         handled = handlePriv(ping);
         break;
      }
      case PONG: {
         const Pong &pong = static_cast<const Pong &>(msg);
         handled = handlePriv(pong);
         break;
      }
      case TRACE: {
         const Trace &trace = static_cast<const Trace &>(msg);
         handled = handlePriv(trace);
         break;
      }
      case BUSY: {
         const Busy &busy = static_cast<const Busy &>(msg);
         handled = handlePriv(busy);
         break;
      }
      case IDLE: {
         const Idle &idle = static_cast<const Idle &>(msg);
         handled = handlePriv(idle);
         break;
      }
      case BARRIER: {
         const Barrier &barrier = static_cast<const Barrier &>(msg);
         handled = handlePriv(barrier);
         break;
      }
      case BARRIERREACHED: {
         const BarrierReached &reached = static_cast<const BarrierReached &>(msg);
         handled = handlePriv(reached);
         registerReply(msg.uuid(), msg);
         break;
      }
      case SETID: {
         const SetId &setid = static_cast<const SetId &>(msg);
         (void)setid;
         break;
      }
      case REPLAYFINISHED: {
         const ReplayFinished &fin = static_cast<const ReplayFinished &>(msg);
         handled = handlePriv(fin);
         break;
      }
      case SENDTEXT: {
         const SendText &info = static_cast<const SendText &>(msg);
         handled = handlePriv(info, pl);
         break;
      }
      case UPDATESTATUS: {
         const UpdateStatus &status = static_cast<const UpdateStatus &>(msg);
         handled = handlePriv(status);
         break;
      }
      case MODULEAVAILABLE: {
         const ModuleAvailable &mod = static_cast<const ModuleAvailable &>(msg);
         handled = handlePriv(mod);
         break;
      }
      case EXECUTIONPROGRESS: {
         break;
      }
      case LOCKUI: {
         break;
      }
      case OBJECTRECEIVEPOLICY: {
         const ObjectReceivePolicy &m = static_cast<const ObjectReceivePolicy &>(msg);
         handled = handlePriv(m);
         break;
      }
      case SCHEDULINGPOLICY: {
         const SchedulingPolicy &m = static_cast<const SchedulingPolicy &>(msg);
         handled = handlePriv(m);
         break;
      }
      case REDUCEPOLICY: {
         const ReducePolicy &m = static_cast<const ReducePolicy &>(msg);
         handled = handlePriv(m);
         break;
      }
      case REQUESTTUNNEL: {
         const RequestTunnel &m = static_cast<const RequestTunnel &>(msg);
         handled = handlePriv(m);
         break;
      }
      case CLOSECONNECTION: {
         const CloseConnection &m = static_cast<const CloseConnection &>(msg);
         handled = handlePriv(m);
         break;
      }

      case FILEQUERY:
      case FILEQUERYRESULT:
      case DATATRANSFERSTATE:
         break;

      default:
         CERR << "message type not handled: " << msg << std::endl;
         assert("message type not handled" == 0);
         break;
   }

   if (handled) {
      if (msg.typeFlags() & TriggerQueue) {
         processQueue();
      }
   } else {
      if (msg.typeFlags() & QueueIfUnhandled) {
          if (payload) {
              auto pl = std::make_shared<const buffer>(payload, payload+payloadSize);
              m_queue.emplace_back(msg, pl);

          } else {
              m_queue.emplace_back(msg, nullptr);
          }
#ifndef NDEBUG
         m_alreadySeen.erase(msg.uuid());
#endif
      }
   }

   return true;
}

void StateTracker::processQueue() {

   if (m_processingQueue)
      return;
   m_processingQueue = true;

   VistleState queue;
   std::swap(m_queue, queue);

   for (auto &m: queue) {
      handle(m.message, m.payload.get());
   }

   m_processingQueue = false;
}

void StateTracker::cleanQueue(int id) {

   using namespace message;

   VistleState queue;
   std::swap(m_queue, queue);

   for (auto &m: queue) {
      auto &msg = m.message;
      if (msg.destId() == id)
          continue;
      switch(msg.type()) {
      case CONNECT: {
          const auto &m = msg.as<Connect>();
          if (m.getModuleA() == id || m.getModuleB() == id)
              continue;
          break;
      }
      case DISCONNECT: {
          const auto &m = msg.as<Disconnect>();
          if (m.getModuleA() == id || m.getModuleB() == id)
              continue;
          break;
      }
      default:
          break;
      }
      m_queue.emplace_back(m);
   }
}

bool StateTracker::handlePriv(const message::AddHub &slave) {
   std::lock_guard<mutex> locker(m_slaveMutex);
   for (auto &h: m_hubs) {
       if (h.id == slave.id()) {
           m_slaveCondition.notify_all();
           return true;
       }
   }
   m_hubs.emplace_back(slave.id(), slave.name());
   m_hubs.back().numRanks = slave.numRanks();
   m_hubs.back().port = slave.port();
   m_hubs.back().dataPort = slave.dataPort();
   if (slave.hasAddress())
      m_hubs.back().address = slave.address();

   // for per-hub parameters
   Module hub(slave.id(), slave.id());
   hub.name = slave.name();
   runningMap.emplace(slave.id(), hub);

   m_slaveCondition.notify_all();
   return true;
}

bool StateTracker::handlePriv(const message::Ping &ping) {

   //CERR << "Ping [" << ping.senderId() << " " << ping.getCharacter() << "]" << std::endl;
   return true;
}

bool StateTracker::handlePriv(const message::Pong &pong) {

   CERR << "Pong [" << pong.senderId() << " " << pong.getCharacter() << "]" << std::endl;
   return true;
}

bool StateTracker::handlePriv(const message::Trace &trace) {

   if (trace.on()) {
      m_traceType = trace.messageType();
      m_traceId = trace.module();
      CERR << "tracing " << m_traceType << " from/to " << m_traceId << std::endl;
   } else {
      CERR << "disabling tracing of " << m_traceType << " from/to " << m_traceId << std::endl;
      m_traceId = Id::Invalid;
      m_traceType = message::INVALID;
   }
   return true;
}

bool StateTracker::handlePriv(const message::Spawn &spawn) {

   ++m_graphChangeCount;

   int moduleId = spawn.spawnId();
   if (moduleId == Id::Invalid) {
      // don't track when master hub has not yet provided a module id
      return true;
   }

   int hub = spawn.hubId();

   auto result = runningMap.emplace(moduleId, Module(moduleId, hub));
   Module &mod = result.first->second;
   assert(hub <= Id::MasterHub);
   mod.hub = hub;
   mod.name = spawn.getName();

   mutex_locker guard(m_stateMutex);
   for (StateObserver *o: m_observers) {
      o->incModificationCount();
      o->newModule(moduleId, spawn.uuid(), mod.name);
   }

   return true;
}

bool StateTracker::handlePriv(const message::Started &started) {

   ++m_graphChangeCount;

   int moduleId = started.senderId();
   auto it = runningMap.find(moduleId);
   if (it == runningMap.end()) {
       CERR << "did not find " << moduleId << " in runningMap, contents are: ";
       for (auto it = runningMap.begin(); it != runningMap.end(); ++it) {
           if (it != runningMap.begin())
               std::cerr << ", ";
           std::cerr << it->first;
       }
       std::cerr << std::endl;
   }
   assert(it != runningMap.end());
   auto &mod = it->second;
   mod.initialized = true;

   mutex_locker guard(m_stateMutex);
   for (StateObserver *o: m_observers) {
      o->moduleStateChanged(moduleId, mod.state());
   }

   return true;
}

bool StateTracker::handleConnect(const message::Connect &connect) {

    if (!handlePriv(connect)) {
        // to be queued by caller
        //m_queue.emplace_back(connect);
        return false;
    }
    return true;
}

bool StateTracker::handleDisconnect(const message::Disconnect &disconnect) {

    if (!handlePriv(disconnect)) {
        // to be queued by caller
        //m_queue.emplace_back(disconnect);
        return false;
    }
    return true;
}

bool StateTracker::handlePriv(const message::Connect &connect) {

   ++m_graphChangeCount;

   bool ret = true;
   if (portTracker()) {
      ret = portTracker()->addConnection(connect.getModuleA(),
            connect.getPortAName(),
            connect.getModuleB(),
            connect.getPortBName());
   }

   computeHeights();

   return ret;
}

bool StateTracker::handlePriv(const message::Disconnect &disconnect) {

   ++m_graphChangeCount;

   bool ret = true;
   if (portTracker()) {
      ret = portTracker()->removeConnection(disconnect.getModuleA(),
            disconnect.getPortAName(),
            disconnect.getModuleB(),
            disconnect.getPortBName());
   }

   computeHeights();

   return ret;
}

bool StateTracker::handlePriv(const message::ModuleExit &moduleExit) {

   ++m_graphChangeCount;

   const int mod = moduleExit.senderId();
   portTracker()->removeModule(mod);

   //CERR << " Module [" << mod << "] quit" << std::endl;

   { 
      RunningMap::iterator it = runningMap.find(mod);
      if (it != runningMap.end()) {
         quitMap.insert(*it);
         runningMap.erase(it);
      } else {
         it = quitMap.find(mod);
         if (it == quitMap.end())
            CERR << " ModuleExit [" << mod << "] not found in map" << std::endl;
      }
   }
   {
      ModuleSet::iterator it = busySet.find(mod);
      if (it != busySet.end())
         busySet.erase(it);
   }

   cleanQueue(mod);

   mutex_locker guard(m_stateMutex);
   for (StateObserver *o: m_observers) {
      o->incModificationCount();
      o->deleteModule(mod);
   }

   return true;
}

bool StateTracker::handlePriv(const message::Execute &execute) {

   return true;
}

bool StateTracker::handlePriv(const message::ExecutionProgress &prog) {

   return true;
}

bool StateTracker::handlePriv(const message::Busy &busy) {

   const int id = busy.senderId();
   if (busySet.find(id) != busySet.end()) {
      //CERR << "module " << id << " sent Busy twice" << std::endl;
   } else {
      busySet.insert(id);
   }
   auto it = runningMap.find(id);
   if (it == runningMap.end()) {
       assert(quitMap.find(id) != quitMap.end());
       return false;
   }
   auto &mod = it->second;
   mod.busy = true;

   mutex_locker guard(m_stateMutex);
   for (StateObserver *o: m_observers) {
      o->moduleStateChanged(id, mod.state());
   }

   return true;
}

bool StateTracker::handlePriv(const message::Idle &idle) {

   const int id = idle.senderId();
   ModuleSet::iterator it = busySet.find(id);
   if (it != busySet.end()) {
      busySet.erase(it);
   } else {
      //CERR << "module " << id << " sent Idle, but was not busy" << std::endl;
   }
   auto rit = runningMap.find(id);
   if (rit == runningMap.end()) {
       assert(quitMap.find(id) != quitMap.end());
       return false;
   }
   auto &mod = rit->second;
   mod.busy = false;

   mutex_locker guard(m_stateMutex);
   for (StateObserver *o: m_observers) {
      o->moduleStateChanged(id, mod.state());
   }

   return true;
}

bool StateTracker::handlePriv(const message::AddParameter &addParam) {

#ifdef DEBUG
   CERR << "AddParameter: module=" << addParam.moduleName() << "(" << addParam.senderId() << "), name=" << addParam.getName() << std::endl;
#endif

   auto mit = runningMap.find(addParam.senderId());
   if (mit == runningMap.end()) {
      CERR << addParam << ": did not find sending module" << std::endl;
      return true;
   }
   assert(mit != runningMap.end());
   auto &mod = mit->second;
   ParameterMap &pm = mod.parameters;
   ParameterOrder &po = mod.paramOrder;
   ParameterMap::iterator it = pm.find(addParam.getName());
   if (it != pm.end()) {
      if (addParam.senderId() == Id::Vistle)
          return true;
      CERR << "duplicate parameter " << addParam.moduleName() << ":" << addParam.getName() << std::endl;
   } else {
      pm[addParam.getName()] = addParam.getParameter();
      int maxIdx = 0;
      const auto rit = po.rbegin();
      if (rit != po.rend())
         maxIdx = rit->first;
      po[maxIdx+1] = addParam.getName();
   }

   mutex_locker guard(m_stateMutex);
   for (StateObserver *o: m_observers) {
      o->incModificationCount();
      o->newParameter(addParam.senderId(), addParam.getName());
   }

   if (portTracker()) {
      const Port *p = portTracker()->addPort(addParam.senderId(), addParam.getName(), addParam.description(), Port::PARAMETER);

      for (StateObserver *o: m_observers) {
         o->newPort(p->getModuleID(), p->getName());
      }
   }

   return true;
}

bool StateTracker::handlePriv(const message::RemoveParameter &removeParam) {

#ifdef DEBUG
   CERR << "RemoveParameter: module=" << removeParam.moduleName() << "(" << removeParam.senderId() << "), name=" << removeParam.getName() << std::endl;
#endif

   auto mit = runningMap.find(removeParam.senderId());
   if (mit == runningMap.end())
       return false;
   assert(mit != runningMap.end());
   auto &mod = mit->second;
   ParameterMap &pm = mod.parameters;
   ParameterOrder &po = mod.paramOrder;
   ParameterMap::iterator it = pm.find(removeParam.getName());
   if (it == pm.end()) {
      CERR << "parameter to be removed not found: " << removeParam.moduleName() << ":" << removeParam.senderId() << ": " << removeParam.getName() << std::endl;
      return false;
   } else {
      pm.erase(it);
      for (auto rit = po.begin(); rit != po.end(); ++rit) {
          if (rit->second == removeParam.getName()) {
              po.erase(rit);
              break;
          }
      }
   }

   mutex_locker guard(m_stateMutex);
   if (portTracker()) {
      for (StateObserver *o: m_observers) {
         o->deletePort(removeParam.senderId(), removeParam.getName());
      }

      portTracker()->removePort(Port(removeParam.senderId(), removeParam.getName(), Port::PARAMETER));
   }

   for (StateObserver *o: m_observers) {
      o->incModificationCount();
      o->deleteParameter(removeParam.senderId(), removeParam.getName());
   }

   return true;
}

bool StateTracker::handlePriv(const message::SetParameter &setParam) {

#ifdef DEBUG
   CERR << "SetParameter: sender=" << setParam.senderId() << ", module=" << setParam.getModule() << ", name=" << setParam.getName() << std::endl;
#endif

   bool handled = false;

   const int senderId = setParam.senderId();
   if (setParam.getModule()==senderId && runningMap.find(senderId) != runningMap.end()) {
      auto param = getParameter(setParam.getModule(), setParam.getName());
      if (param) {
         setParam.apply(param);
         handled = true;
      }
   }

   if (handled) {
      mutex_locker guard(m_stateMutex);
      for (StateObserver *o: m_observers) {
         o->incModificationCount();
         o->parameterValueChanged(setParam.senderId(), setParam.getName());
      }
   }

   return handled;
}

bool StateTracker::handlePriv(const message::SetParameterChoices &choices, const buffer &payload) {

   const int senderId = choices.senderId();
   if (runningMap.find(senderId) == runningMap.end())
       return false;

   auto p = getParameter(choices.senderId(), choices.getName());
   if (!p)
      return false;

   auto pl = message::getPayload<message::SetParameterChoices::Payload>(payload);

   choices.apply(p, pl);

   //CERR << "choices changed for " << choices.getModule() << ":" << choices.getName() << ": #" << p->choices().size() << std::endl;

   mutex_locker guard(m_stateMutex);
   for (StateObserver *o: m_observers) {
      o->incModificationCount();
      o->parameterChoicesChanged(choices.senderId(), choices.getName());
   }

   return true;
}

bool StateTracker::handlePriv(const message::Quit &quit) {

   mutex_locker guard(m_stateMutex);
   for (StateObserver *o: m_observers) {
      o->quitRequested();
   }

   return true;
}

bool StateTracker::handlePriv(const message::Kill &kill) {

    const int destId = kill.getModule();
    std::set<int> ids;
    if (destId == message::Id::Broadcast) {
        for (const auto &p: runningMap) {
            ids.insert(p.first);
        }
        for (const auto &p: quitMap) {
            ids.insert(p.first);
        }
    } else {
        ids.insert(destId);
    }

    for (auto id: ids) {
        auto it = runningMap.find(id);
        if (it == runningMap.end()) {
            it = quitMap.find(id);
            assert(it != quitMap.end());
        }

        auto &mod = it->second;
        mod.killed = true;

        mutex_locker guard(m_stateMutex);
        for (StateObserver *o: m_observers) {
            o->moduleStateChanged(id, mod.state());
        }
    }

    return true;
}

bool StateTracker::handlePriv(const message::AddObject &addObj) {

   ++m_numObjects;
   return true;
}

bool StateTracker::handlePriv(const message::Barrier &barrier) {

   return true;
}

bool StateTracker::handlePriv(const message::BarrierReached &barrReached) {

   return true;
}

bool StateTracker::handlePriv(const message::AddPort &createPort) {

   if (portTracker()) {
      const Port * p = portTracker()->addPort(createPort.getPort());

      if (!p)
         return false;

      mutex_locker guard(m_stateMutex);
      for (StateObserver *o: m_observers) {
         o->incModificationCount();
         o->newPort(p->getModuleID(), p->getName());
      }
   }

   return true;
}

bool StateTracker::handlePriv(const message::RemovePort &destroyPort) {

   if (portTracker()) {
      Port p = destroyPort.getPort();
      int id = p.getModuleID();
      std::string name = p.getName();

      if (portTracker()->findPort(p)) {

          mutex_locker guard(m_stateMutex);
          for (StateObserver *o: m_observers) {
              o->incModificationCount();
              o->deletePort(id, name);
          }
          portTracker()->removePort(p);
      }
   }

   return true;
}

bool StateTracker::handlePriv(const message::ReplayFinished &reset)
{
   mutex_locker guard(m_stateMutex);
   for (StateObserver *o: m_observers) {
      o->resetModificationCount();
   }
   return true;
}

bool StateTracker::handlePriv(const message::SendText &info, const buffer &payload)
{
    auto pl = message::getPayload<message::SendText::Payload>(payload);
    mutex_locker guard(m_stateMutex);
    for (StateObserver *o: m_observers) {
        o->info(pl.text, info.textType(), info.senderId(), info.rank(), info.referenceType(), info.referenceUuid());
    }

    return true;
}

bool StateTracker::handlePriv(const message::UpdateStatus &status) {

    if (status.statusType() == message::UpdateStatus::LoadedFile) {
        m_loadedWorkflowFile = status.text();
        mutex_locker guard(m_stateMutex);
        for (StateObserver *o: m_observers) {
            o->loadedWorkflowChanged(m_loadedWorkflowFile);
        }

        return true;
    }

    auto it = runningMap.find(status.senderId());
    if (it == runningMap.end())
        return false;

    auto &mod = it->second;
    mod.statusText = status.text();
    mod.statusImportance = status.importance();
    mod.statusTime = m_statusTime;
    ++m_statusTime;

    mutex_locker guard(m_stateMutex);
    for (StateObserver *o: m_observers) {
        o->status(mod.id, mod.statusText, mod.statusImportance);
    }

    if (mod.statusText.empty() || mod.statusImportance >= m_currentStatusImportance) {
        auto oid = m_currentStatusId;
        auto otext = m_currentStatus;
        auto oprio = m_currentStatusImportance;
        updateStatus();
        if (oid != m_currentStatusId || otext != m_currentStatus || oprio != m_currentStatusImportance) {
            for (StateObserver *o: m_observers) {
                o->updateStatus(m_currentStatusId, m_currentStatus, m_currentStatusImportance);
            }
        }
    }

    return true;
}

bool StateTracker::handlePriv(const message::ModuleAvailable &avail) {

    if (avail.hub() == Id::Invalid)
        return true;

    AvailableModule mod;
    mod.hub = avail.hub();
    mod.name = avail.name();
    mod.path = avail.path();

    AvailableModule::Key key(mod.hub, mod.name);

    mutex_locker guard(m_stateMutex);
    if (m_availableModules.emplace(key, mod).second) {
        for (StateObserver *o: m_observers) {
            o->moduleAvailable(mod.hub, mod.name, mod.path);
        }
    } else {
        CERR << "Duplicate module: " << mod.hub << " " << mod.name << std::endl;
    }

    return true;
}

bool StateTracker::handlePriv(const message::ObjectReceivePolicy &receivePolicy)
{
   const int id = receivePolicy.senderId();
   RunningMap::iterator it = runningMap.find(id);
   if (it == runningMap.end()) {
      CERR << " Module [" << id << "] changed ObjectReceivePolicy, but not found in running map" << std::endl;
      return false;
   }
   Module &mod = it->second;
   mod.objectPolicy = receivePolicy.policy();
   return true;
}

bool StateTracker::handlePriv(const message::SchedulingPolicy &schedulingPolicy)
{
   const int id = schedulingPolicy.senderId();
   RunningMap::iterator it = runningMap.find(id);
   if (it == runningMap.end()) {
      CERR << " Module [" << id << "] changed SchedulingPolicy, but not found in running map" << std::endl;
      return false;
   }
   Module &mod = it->second;
   mod.schedulingPolicy = schedulingPolicy.policy();
   return true;
}

bool StateTracker::handlePriv(const message::ReducePolicy &reducePolicy)
{
   const int id = reducePolicy.senderId();
   RunningMap::iterator it = runningMap.find(id);
   if (it == runningMap.end()) {
      CERR << " Module [" << id << "] changed ReducePolicy, but not found in running map" << std::endl;
      return false;
   }
   Module &mod = it->second;
   mod.reducePolicy = reducePolicy.policy();
   return true;
}

bool StateTracker::handlePriv(const message::RequestTunnel &tunnel)
{
   return true;
}

bool StateTracker::handlePriv(const message::CloseConnection &close)
{
   CERR << "socket shutdown requested: " << close.reason() << std::endl;
   return true;
}


StateTracker::~StateTracker() {

    if (m_portTracker) {
        m_portTracker->setTracker(nullptr);
    }
}

std::shared_ptr<PortTracker> StateTracker::portTracker() const {

   return m_portTracker;
}

std::vector<std::string> StateTracker::getParameters(int id) const {

   std::vector<std::string> result;

   mutex_locker guard(m_stateMutex);

   RunningMap::const_iterator rit = runningMap.find(id);
   if (rit == runningMap.end())
      return result;

   const ParameterOrder &po = rit->second.paramOrder;
   BOOST_FOREACH (ParameterOrder::value_type val, po) {
      const auto &name = val.second;
      result.push_back(name);
   }

   return result;
}

std::shared_ptr<Parameter> StateTracker::getParameter(int id, const std::string &name) const {

   mutex_locker guard(m_stateMutex);

   RunningMap::const_iterator rit = runningMap.find(id);
   if (rit == runningMap.end())
      return std::shared_ptr<Parameter>();

   ParameterMap::const_iterator pit = rit->second.parameters.find(name);
   if (pit == rit->second.parameters.end())
      return std::shared_ptr<Parameter>();

   return pit->second;
}

bool StateTracker::registerRequest(const message::uuid_t &uuid) {

   std::lock_guard<mutex> locker(m_replyMutex);

   auto it = m_outstandingReplies.find(uuid);
   if (it != m_outstandingReplies.end()) {
      CERR << "duplicate attempt to wait for reply" << std::endl;
      return false;
   }

   //CERR << "waiting for " << uuid  << std::endl;
   m_outstandingReplies[uuid] = std::shared_ptr<message::Buffer>();
   return true;
}

std::shared_ptr<message::Buffer> StateTracker::waitForReply(const message::uuid_t &uuid) {

   std::unique_lock<mutex> locker(m_replyMutex);
   std::shared_ptr<message::Buffer> ret = removeRequest(uuid);
   while (!ret) {
      m_replyCondition.wait(locker);
      ret = removeRequest(uuid);
   }
   return ret;
}

std::shared_ptr<message::Buffer> StateTracker::removeRequest(const message::uuid_t &uuid) {

   //CERR << "remove request try: " << uuid << std::endl;
   std::shared_ptr<message::Buffer> ret;
   auto it = m_outstandingReplies.find(uuid);
   if (it != m_outstandingReplies.end() && it->second) {
      ret = it->second;
      //CERR << "remove request success: " << uuid << std::endl;
      m_outstandingReplies.erase(it);
   }
   return ret;
}

bool StateTracker::registerReply(const message::uuid_t &uuid, const message::Message &msg) {

   std::lock_guard<mutex> locker(m_replyMutex);
   auto it = m_outstandingReplies.find(uuid);
   if (it == m_outstandingReplies.end()) {
      return false;
   }
   if (it->second) {
      CERR << "attempt to register duplicate reply for " << uuid << std::endl;
      assert(!it->second);
      return false;
   }

   it->second.reset(new message::Buffer(msg));

   //CERR << "notifying all for " << uuid  << " and " << m_outstandingReplies.size() << " others" << std::endl;

   m_replyCondition.notify_all();

   return true;
}

std::vector<int> StateTracker::waitForSlaveHubs(size_t count) {

    std::unique_lock<mutex> locker(m_slaveMutex);
   auto hubIds = getSlaveHubs();
   while (hubIds.size() < count) {
      m_slaveCondition.wait(locker);
      hubIds = getSlaveHubs();
   }
   return hubIds;
}

std::vector<int> StateTracker::waitForSlaveHubs(const std::vector<std::string> &names) {

   auto findAll = [this](const std::vector<std::string> &names, std::vector<int> &ids) -> bool {
      const auto hubIds = getSlaveHubs();
      std::vector<std::string> available;
      for (int id: hubIds)
         available.push_back(hubName(id));
      
      ids.clear();
      size_t found=0;
      for (const auto &name: names) {
         for (const auto &slave: m_hubs) {
            if (slave.id == Id::MasterHub)
               continue;
            if (name == slave.name) {
               ++found;
               ids.push_back(slave.id);
            } else {
               ids.push_back(Id::Invalid);
            }
         }
      }
      return found == names.size();
   };

   std::unique_lock<mutex> locker(m_slaveMutex);
   std::vector<int> ids;
   while (!findAll(names, ids)) {
      m_slaveCondition.wait(locker);
   }
   return ids;
}

void StateTracker::registerObserver(StateObserver *observer) {

   mutex_locker guard(m_stateMutex);
   m_observers.insert(observer);
}

ParameterSet StateTracker::getConnectedParameters(const Parameter &param) const {

   mutex_locker guard(m_stateMutex);

   std::function<ParameterSet (const Port *, ParameterSet)> findAllConnectedPorts;
   findAllConnectedPorts = [this, &findAllConnectedPorts] (const Port *port, ParameterSet conn) -> ParameterSet {
      if (const Port::ConstPortSet *list = portTracker()->getConnectionList(port)) {
         for (auto port: *list) {
            auto param = getParameter(port->getModuleID(), port->getName());
            if (param && conn.find(param) == conn.end()) {
               conn.insert(param);
               const Port *port = portTracker()->getPort(param->module(), param->getName());
               conn = findAllConnectedPorts(port, conn);
            }
         }
      }
      return conn;
   };

   if (!portTracker())
      return ParameterSet();
   const Port *port = portTracker()->findPort(param.module(), param.getName());
   if (!port)
      return ParameterSet();
   if (port->getType() != Port::PARAMETER)
      return ParameterSet();
   return findAllConnectedPorts(port, ParameterSet());
}

void StateTracker::computeHeights() {

   std::set<Module *> modules;
   for (auto &mod: runningMap) {
      mod.second.height = -1;
      modules.insert(&mod.second);
   }

   while (!modules.empty()) {
      for (auto mod: modules) {
         int id = mod->id;
         auto outputs = portTracker()->getOutputPorts(id);

         int height = -1;
         bool isSink = true;
         for (auto &output: outputs) {
            for (auto &port: output->connections()) {
               const int otherId = port->getModuleID();
               auto it = runningMap.find(otherId);
               if (it == runningMap.end()) {
                  if (quitMap.find(otherId) == quitMap.end())
                     CERR << "did not find module " << otherId << ", connected to " << id << " at port " << output->getName() << std::endl;
                  continue;
               }
               isSink = false;
               assert(it != runningMap.end());
               const auto &otherMod = it->second;
               if (otherMod.height != -1 && (height == -1 || otherMod.height+1 < height)) {
                  height = otherMod.height + 1;
               }
            }
         }
         if (isSink) {
            height = 0;
         }
         if (height != -1) {
            mod->height = height;
            modules.erase(mod);
            break;
         }
      }
   }
}

int StateTracker::graphChangeCount() const {

    return m_graphChangeCount;
}

std::string StateTracker::loadedWorkflowFile() const {

    return m_loadedWorkflowFile;
}

std::string StateTracker::statusText() const {

    return m_currentStatus;
}

void StateTracker::updateStatus() {

    using namespace message;

    m_currentStatusImportance = UpdateStatus::Bulk;
    m_currentStatus.clear();
    m_currentStatusId = Id::Invalid;
    unsigned long time = 0;
    bool system = false;

    for (auto &p: runningMap) {
        auto &mod = p.second;
        if (!mod.statusText.empty()) {
            if (mod.statusImportance > m_currentStatusImportance) {
                m_currentStatusImportance = mod.statusImportance;
                time = mod.statusTime;
                m_currentStatusId = mod.id;
                system = !Id::isModule(mod.id);
            }
            if (mod.statusImportance == m_currentStatusImportance) {
                if ((mod.statusTime>=time && (!system || !Id::isModule(mod.id)))
                        || (!system && !Id::isModule(mod.id))) {
                    time = mod.statusTime;
                    m_currentStatus = mod.statusText;
                    m_currentStatusId = mod.id;
                    system = !Id::isModule(mod.id);
                }
            }
        }
    }
}

void StateObserver::quitRequested() {

}

void StateObserver::incModificationCount() {
   ++m_modificationCount;
}

long StateObserver::modificationCount() const {
    return m_modificationCount;
}

void StateObserver::loadedWorkflowChanged(const std::string &filename) {

}

void StateObserver::resetModificationCount() {
   m_modificationCount = 0;
}

} // namespace vistle
