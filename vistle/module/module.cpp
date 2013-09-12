#include <mpi.h>

#include <cstdlib>
#include <cstdio>

#include <sys/types.h>
#ifndef _WIN32
#include <unistd.h>
#else
#include<Winsock2.h>
#pragma comment(lib, "Ws2_32.lib")
#endif

#include <sstream>
#include <iostream>
#include <iomanip>
#include <algorithm>

#include <boost/interprocess/shared_memory_object.hpp>
#include <boost/interprocess/mapped_region.hpp>
#include <boost/mpl/vector.hpp>
#include <boost/mpl/for_each.hpp>
#include <boost/mpl/transform.hpp>

#include <core/object.h>
#include <core/message.h>
#include <core/messagequeue.h>
#include <core/parameter.h>
#include <core/shm.h>
#include <core/objectcache.h>
#include <core/port.h>

#include "module.h"

using namespace boost::interprocess;

namespace vistle {

template<typename CharT, typename TraitsT = std::char_traits<CharT> >
class msgstreambuf: public std::basic_streambuf<CharT, TraitsT> {

 public:
   msgstreambuf(Module *mod)
   : m_module(mod)
   {}

   int overflow(int ch) {
      if (ch != EOF) {
         m_buf.push_back(ch);
         if (ch == '\n')
            flush();
         return 0;
      } else {
         return EOF;
      }
   }

   std::streamsize xsputn (const CharT *s, std::streamsize num) {
      size_t end = m_buf.size();
      m_buf.resize(end+num);
      memcpy(m_buf.data()+end, s, num);
      flush();
      return num;
   }

   void flush() {
      m_module->sendInfo("%s", std::string(m_buf.data(), m_buf.size()).c_str());
      m_buf.clear();
   }

 private:
   Module *m_module = nullptr;
   std::vector<char> m_buf;
};



Module::Module(const std::string &n, const std::string &shmname,
      const unsigned int r,
      const unsigned int s, const int m)
: m_name(n)
, m_rank(r)
, m_size(s)
, m_id(m)
, m_executionCount(0)
, m_defaultCacheMode(ObjectCache::CacheNone)
, m_syncMessageProcessing(false)
{
#ifdef _WIN32
    WSADATA wsaData;
    WSAStartup(MAKEWORD(2, 2), &wsaData);
#endif

   message::DefaultSender::init(m_id, m_rank);

   const int HOSTNAMESIZE = 64;
   char hostname[HOSTNAMESIZE];
   gethostname(hostname, HOSTNAMESIZE - 1);

   try {
      Shm::attach(shmname, id(), rank(), sendMessageQueue);
   } catch (interprocess_exception &ex) {
      std::cerr << "module " << id() << " [" << rank() << "/" << size() << "] "
                << ex.what() << std::endl;
      exit(2);
   }

   // names are swapped relative to communicator
   std::string smqName = message::MessageQueue::createName("rmq", id(), rank());
   try {
      sendMessageQueue = message::MessageQueue::open(smqName);
   } catch (interprocess_exception &ex) {
      std::cerr << "module " << id() << " [" << rank() << "/" << size() << "] "
         << ex.what() << std::endl;
      exit(2);
   }

   std::string rmqName = message::MessageQueue::createName("smq", id(), rank());
   try {
      receiveMessageQueue = message::MessageQueue::open(rmqName);
   } catch (interprocess_exception &ex) {
      std::cerr << "module " << id() << " [" << rank() << "/" << size() << "] "
         << ex.what() << std::endl;
      exit(2);
   }

   std::cerr << "  module [" << name() << "] [" << id() << "] [" << rank()
             << "/" << size() << "] started as " << hostname << ":"
#ifndef _WIN32
             << getpid() << std::endl;
#else
             << std::endl;
#endif

   Parameter *cm = addIntParameter("_cache_mode", "input object caching", ObjectCache::CacheDefault, Parameter::Choice);
   std::vector<std::string> modes;
   assert(ObjectCache::CacheDefault == 0);
   modes.push_back("default");
   assert(ObjectCache::CacheNone == 1);
   modes.push_back("none");
   assert(ObjectCache::CacheAll == 2);
   modes.push_back("all");
   setParameterChoices(cm, modes);

   addVectorParameter("_position", "position in GUI", ParamVector(0., 0.));
}

void Module::initDone() {

   m_streambuf = new msgstreambuf<char>(this);
   m_origStreambuf = std::cerr.rdbuf(m_streambuf);

   for (auto &pair: parameters) {
      parameterChanged(pair.second);
   }

   sendMessage(message::Started(name()));
}

const std::string &Module::name() const {

   return m_name;
}

int Module::id() const {

   return m_id;
}

unsigned int Module::rank() const {

   return m_rank;
}

unsigned int Module::size() const {

   return m_size;
}

bool Module::syncMessageProcessing() const {

   return m_syncMessageProcessing;
}

void Module::setSyncMessageProcessing(bool sync) {

   m_syncMessageProcessing = sync;
}

ObjectCache::CacheMode Module::setCacheMode(ObjectCache::CacheMode mode, bool updateParam) {

   if (mode == ObjectCache::CacheDefault)
      m_cache.setCacheMode(m_defaultCacheMode);
   else
      m_cache.setCacheMode(mode);

   if (updateParam)
      setIntParameter("_cache_mode", mode);

   return m_cache.cacheMode();
}

void Module::setDefaultCacheMode(ObjectCache::CacheMode mode) {

   assert(mode != ObjectCache::CacheDefault);
   m_defaultCacheMode = mode;
   setCacheMode(m_defaultCacheMode, false);
}


ObjectCache::CacheMode Module::cacheMode(ObjectCache::CacheMode mode) const {

   return m_cache.cacheMode();
}

Port *Module::createInputPort(const std::string &name, const std::string &description, const int flags) {

   std::map<std::string, Port*>::iterator i = inputPorts.find(name);

   if (i == inputPorts.end()) {

      Port *p = new Port(id(), name, Port::INPUT, flags);
      inputPorts[name] = p;

      message::CreatePort message(p);
      sendMessageQueue->getMessageQueue().send(&message, sizeof(message), 0);
      return p;
   }

   return NULL;
}

Port *Module::createOutputPort(const std::string &name, const std::string &description, const int flags) {

   std::map<std::string, Port *>::iterator i = outputPorts.find(name);

   if (i == outputPorts.end()) {

      Port *p = new Port(id(), name, Port::OUTPUT, flags);
      outputPorts[name] = p;

      message::CreatePort message(p);
      sendMessage(message);
      return p;
   }
   return NULL;
}

Port *Module::findInputPort(const std::string &name) const {

   std::map<std::string, Port *>::const_iterator i = inputPorts.find(name);

   if (i == inputPorts.end())
      return NULL;

   return i->second;
}

Port *Module::findOutputPort(const std::string &name) const {

   std::map<std::string, Port *>::const_iterator i = outputPorts.find(name);

   if (i == outputPorts.end())
      return NULL;

   return i->second;
}


bool Module::addParameterGeneric(const std::string &name, Parameter *param, Parameter::Presentation presentation) {

   std::map<std::string, Parameter *>::iterator i =
      parameters.find(name);

   if (i != parameters.end())
      return false;

   parameters[name] = param;

   message::AddParameter add(name, param->description(), param->type(), presentation, m_name);
   sendMessage(add);
   message::SetParameter set(id(), name, param);
   set.setInit();
   set.setUuid(add.uuid());
   sendMessage(set);

   return true;
}

bool Module::updateParameter(const std::string &name, const Parameter *param, const message::SetParameter *inResponseTo, Parameter::RangeType rt) {

   std::map<std::string, Parameter *>::iterator i = parameters.find(name);

   if (i == parameters.end()) {
      std::cerr << "setParameter: " << name << " not found" << std::endl;
      return false;
   }

   if (i->second->type() != param->type()) {
      std::cerr << "setParameter: type mismatch for " << name << " " << i->second->type() << " != " << param->type() << std::endl;
      return false;
   }

   if (&*i->second != param) {
      std::cerr << "setParameter: pointer mismatch for " << name << std::endl;
      return false;
   }

   message::SetParameter set(id(), name, param, rt);
   if (inResponseTo) {
      set.setReply();
      set.setUuid(inResponseTo->uuid());
   }
   sendMessage(set);

   return true;
}

void Module::setParameterChoices(const std::string &name, const std::vector<std::string> &choices)
{
   Parameter *p = findParameter(name);
   if (p)
      setParameterChoices(p, choices);
}

void Module::setParameterChoices(Parameter *param, const std::vector<std::string> &choices)
{
   if (choices.size() <= message::param_num_choices) {
      message::SetParameterChoices sc(id(), param->getName(), choices);
      sendMessage(sc);
   }
}

template<class T>
Parameter *Module::addParameter(const std::string &name, const std::string &description, const T &value, Parameter::Presentation pres) {

   Parameter *p = new ParameterBase<T>(id(), name, value);
   p->setDescription(description);
   if (!addParameterGeneric(name, p, pres)) {
      delete p;
      return NULL;
   }
   return p;
}

Parameter *Module::findParameter(const std::string &name) const {

   std::map<std::string, Parameter *>::const_iterator i = parameters.find(name);

   if (i == parameters.end())
      return NULL;

   return i->second;
}

template<class T>
bool Module::setParameter(const std::string &name, const T &value, const message::SetParameter *inResponseTo) {

   ParameterBase<T> *p = dynamic_cast<ParameterBase<T> *>(findParameter(name));
   if (!p)
      return false;

   return setParameter(p, value, inResponseTo);
}

template<class T>
bool Module::setParameter(ParameterBase<T> *param, const T &value, const message::SetParameter *inResponseTo) {

   param->setValue(value);
   parameterChanged(param);
   return updateParameter(param->getName(), param, inResponseTo);
}

template<class T>
bool Module::setParameterMinimum(ParameterBase<T> *param, const T &minimum) {

   return Module::setParameterRange(param, minimum, param->maximum());
}

template<class T>
bool Module::setParameterMaximum(ParameterBase<T> *param, const T &maximum) {

   return Module::setParameterRange(param, param->minimum(), maximum);
}

template<class T>
bool Module::setParameterRange(const std::string &name, const T &minimum, const T &maximum) {

   ParameterBase<T> *p = dynamic_cast<ParameterBase<T> *>(findParameter(name));
   if (!p)
      return false;

   return setParameterRange(p, minimum, maximum);
}

template<class T>
bool Module::setParameterRange(ParameterBase<T> *param, const T &minimum, const T &maximum) {

   bool ok = true;

   if (maximum > minimum) {
      param->setMinimum(minimum);
      param->setMaximum(maximum);
   } else {
      std::cerr << "range for parameter " << param->getName() << " swapped" << std::endl;
      param->setMinimum(maximum);
      param->setMaximum(minimum);
   }
   T value = param->getValue();
   if (value < param->minimum()) {
      std::cerr << "value for parameter " << param->getName() << " increased to minimum: " << param->minimum() << std::endl;
      param->setValue(param->minimum());
      ok &= updateParameter(param->getName(), param, nullptr);
   }
   if (value > param->maximum()) {
      std::cerr << "value for parameter " << param->getName() << " decreased to maximum: " << param->maximum() << std::endl;
      param->setValue(param->maximum());
      ok &= updateParameter(param->getName(), param, nullptr);
   }
   ok &= updateParameter(param->getName(), param, nullptr, Parameter::Minimum);
   ok &= updateParameter(param->getName(), param, nullptr, Parameter::Maximum);
   return ok;
}

template<class T>
bool Module::getParameter(const std::string &name, T &value) const {

   if (ParameterBase<T> *p = dynamic_cast<ParameterBase<T> *>(findParameter(name))) {
      value = p->getValue();
   } else {
      std::cerr << "Module::getParameter(" << name << "): type failure" << std::endl;
      assert("dynamic_cast failed" == 0);
      return false;
   }

   return true;
}

StringParameter *Module::addStringParameter(const std::string & name, const std::string &description,
                              const std::string & value, Parameter::Presentation p) {

   return dynamic_cast<StringParameter *>(addParameter(name, description, value, p));
}

bool Module::setStringParameter(const std::string & name,
                              const std::string & value, const message::SetParameter *inResponseTo) {

   return setParameter(name, value, inResponseTo);
}

std::string Module::getStringParameter(const std::string & name) const {

   std::string value = "";
   getParameter(name, value);
   return value;
}

FloatParameter *Module::addFloatParameter(const std::string & name, const std::string &description,
                               const Float value) {

   return dynamic_cast<FloatParameter *>(addParameter(name, description, value));
}

bool Module::setFloatParameter(const std::string & name,
                               const Float value, const message::SetParameter *inResponseTo) {

   return setParameter(name, value, inResponseTo);
}

Float Module::getFloatParameter(const std::string & name) const {

   Float value = 0.;
   getParameter(name, value);
   return value;
}

IntParameter *Module::addIntParameter(const std::string & name, const std::string &description,
                             const Integer value, Parameter::Presentation p) {

   return dynamic_cast<IntParameter *>(addParameter(name, description, value, p));
}

bool Module::setIntParameter(const std::string & name,
                             Integer value, const message::SetParameter *inResponseTo) {

   if (name == "_cache_mode") {
      switch (value) {
      case ObjectCache::CacheNone:
      case ObjectCache::CacheAll:
      case ObjectCache::CacheDefault:
         break;
      default:
         value = ObjectCache::CacheDefault;
         break;
      }

      setCacheMode(ObjectCache::CacheMode(value), false);
   }

   return setParameter(name, value, inResponseTo);
}

Integer Module::getIntParameter(const std::string & name) const {

   Integer value = 0;
   getParameter(name, value);
   return value;
}

VectorParameter *Module::addVectorParameter(const std::string & name, const std::string &description,
                                const ParamVector & value) {

   return dynamic_cast<VectorParameter *>(addParameter(name, description, value));
}

bool Module::setVectorParameter(const std::string & name,
                                const ParamVector & value, const message::SetParameter *inResponseTo) {

   return setParameter(name, value, inResponseTo);
}

ParamVector Module::getVectorParameter(const std::string & name) const {

   ParamVector value;
   getParameter(name, value);
   return value;
}

void Module::updateMeta(vistle::Object::ptr obj) const {

   obj->setCreator(id());
   obj->setExecutionCounter(m_executionCount);
}

bool Module::addObject(const std::string &portName, vistle::Object::ptr object) {

   std::map<std::string, Port *>::iterator i = outputPorts.find(portName);
   if (i != outputPorts.end()) {
      return addObject(i->second, object);
   }
   std::cerr << "Module::addObject: output port " << portName << " not found" << std::endl;
   assert(i != outputPorts.end());
   return false;
}

bool Module::addObject(Port *port, vistle::Object::ptr object) {

   updateMeta(object);
   vistle::Object::const_ptr cobj = object;
   return passThroughObject(port, cobj);
}

bool Module::passThroughObject(const std::string &portName, vistle::Object::const_ptr object) {
   std::map<std::string, Port *>::iterator i = outputPorts.find(portName);
   if (i != outputPorts.end()) {
      return passThroughObject(i->second, object);
   }
   std::cerr << "Module::passThroughObject: output port " << portName << " not found" << std::endl;
   assert(i != outputPorts.end());
   return false;
}

bool Module::passThroughObject(Port *port, vistle::Object::const_ptr object) {

   if (!object)
      return false;

   assert(object->check());

   message::AddObject message(port->getName(), object);
   sendMessageQueue->getMessageQueue().send(&message, sizeof(message), 0);
   return true;
}

ObjectList Module::getObjects(const std::string &portName) {

   ObjectList objects;
   std::map<std::string, Port *>::iterator i = inputPorts.find(portName);

   if (i != inputPorts.end()) {

      std::list<shm_handle_t>::iterator shmit;
      ObjectList &olist = i->second->objects();
      for (ObjectList::const_iterator it = olist.begin(); it != olist.end(); it++) {
         Object::const_ptr object = *it;
         if (object.get())
            assert(object->check());
            objects.push_back(object);
      }
   } else {
      std::cerr << "Module::getObjects: input port " << portName << " not found" << std::endl;
      assert(i != inputPorts.end());
   }

   return objects;
}

void Module::removeObject(const std::string &portName, vistle::Object::const_ptr object) {

   bool erased = false;
   shm_handle_t handle = object->getHandle();
   std::map<std::string, Port *>::iterator i = inputPorts.find(portName);

   if (i != inputPorts.end()) {
      ObjectList &olist = i->second->objects();

      for (ObjectList::iterator it = olist.begin(); it != olist.end(); ) {
         if (handle == (*it)->getHandle()) {
            erased = true;
            //object->unref(); // XXX: doesn't erasing the it handle that already?
            it = olist.erase(it);
         } else
            ++it;
      }
      if (!erased)
         std::cerr << "Module " << id() << " removeObject didn't find"
            " object [" << object->getName() << "]" << std::endl;
   } else {
      std::cerr << "Module " << id() << " removeObject didn't find port ["
                << portName << "]" << std::endl;

      assert(i != inputPorts.end());
   }
}

bool Module::hasObject(const std::string &portName) const {

   std::map<std::string, Port *>::const_iterator i = inputPorts.find(portName);

   if (i != inputPorts.end()) {

      return !i->second->objects().empty();
   }

   std::cerr << "Module::hasObject: input port " << portName << " not found" << std::endl;
   assert(i != inputPorts.end());

   return false;
}

vistle::Object::const_ptr Module::takeFirstObject(const std::string &portName) {

   std::map<std::string, Port *>::iterator i = inputPorts.find(portName);

   if (i == inputPorts.end()) {
      std::cerr << "Module::takeFirstObject: input port " << portName << " not found" << std::endl;
      assert(i != inputPorts.end());
      return vistle::Object::ptr();
   }

   if (!i->second->objects().empty()) {

      Object::const_ptr obj = i->second->objects().front();
      assert(obj->check());
      i->second->objects().pop_front();
      return obj;
   }

   return vistle::Object::ptr();
}

bool Module::addInputObject(const std::string & portName,
                            Object::const_ptr object) {

   if (!object)
      return false;

   if (object)
      assert(object->check());

   if (m_executionCount < object->getExecutionCounter())
      m_executionCount = object->getExecutionCounter();

   Port *p = findInputPort(portName);

   if (p) {
      m_cache.addObject(portName, object);
      p->objects().push_back(object);
      return true;
   }

   std::cerr << "Module::addInputObject: input port " << portName << " not found" << std::endl;
   assert(p);

   return false;
}

bool Module::isConnected(const std::string &portname) const {

   Port *p = findInputPort(portname);
   if (!p)
      p = findOutputPort(portname);
   if (!p)
      return false;

   return !p->connections().empty();
}

bool Module::parameterChanged(Parameter *p) {

   (void)p;
   return true;
}

bool Module::parameterAdded(const int senderId, const std::string &name, const message::AddParameter &msg, const std::string &moduleName) {

   (void)senderId;
   (void)name;
   (void)msg;
   (void)moduleName;
   return false;
}

bool Module::parameterChanged(const int senderId, const std::string &name, const message::SetParameter &msg) {

   (void)senderId;
   (void)name;
   (void)msg;
   return false;
}


bool Module::dispatch() {

   size_t msgSize;
   unsigned int priority;
   char msgRecvBuf[message::Message::MESSAGE_SIZE];

   bool again = true;
   receiveMessageQueue->getMessageQueue().receive(
                                               (void *) msgRecvBuf,
                                               message::Message::MESSAGE_SIZE,
                                               msgSize, priority);
   vistle::message::Message *message = (vistle::message::Message *) msgRecvBuf;


   if (syncMessageProcessing()) {
      int sync = 0, allsync = 0;

      switch (message->type()) {
         case vistle::message::Message::OBJECTRECEIVED:
         case vistle::message::Message::QUIT:
            sync = 1;
            break;
         default:
            break;
      }

      MPI_Allreduce(&sync, &allsync, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

      do {
         vistle::message::Message *message = (vistle::message::Message *) msgRecvBuf;

         switch (message->type()) {
            case vistle::message::Message::OBJECTRECEIVED:
            case vistle::message::Message::QUIT:
               sync = 1;
               break;
            default:
               break;
         }

         again &= handleMessage(message);

         if (allsync && !sync) {
            receiveMessageQueue->getMessageQueue().receive(
                  (void *) msgRecvBuf,
                  message::Message::MESSAGE_SIZE,
                  msgSize, priority);
         }

      } while(allsync && !sync);
   } else {

      again &= handleMessage(message);
   }

   return again;
}


void Module::sendMessage(const message::Message &message) const {

   sendMessageQueue->getMessageQueue().send(&message, message.size(), 0);
}


bool Module::handleMessage(const vistle::message::Message *message) {

   switch (message->type()) {

      case vistle::message::Message::PING: {

         const vistle::message::Ping *ping =
            static_cast<const vistle::message::Ping *>(message);

         std::cerr << "    module [" << name() << "] [" << id() << "] ["
                   << rank() << "/" << size() << "] ping ["
                   << ping->getCharacter() << "]" << std::endl;
         vistle::message::Pong m(ping->getCharacter(), ping->senderId());
         m.setUuid(ping->uuid());
         sendMessage(m);
         break;
      }

      case vistle::message::Message::PONG: {

         const vistle::message::Pong *pong =
            static_cast<const vistle::message::Pong *>(message);

         std::cerr << "    module [" << name() << "] [" << id() << "] ["
                   << rank() << "/" << size() << "] pong ["
                   << pong->getCharacter() << "]" << std::endl;
         break;
      }

      case message::Message::QUIT: {

         const message::Quit *quit =
            static_cast<const message::Quit *>(message);
         //TODO: uuid should be included in coresponding ModuleExit message
         (void) quit;
         return false;
         break;
      }

      case message::Message::KILL: {

         const message::Kill *kill =
            static_cast<const message::Kill *>(message);
         //TODO: uuid should be included in coresponding ModuleExit message
         if (kill->getModule() == id()) {
            return false;
         } else {
            std::cerr << "module [" << name() << "] [" << id() << "] ["
               << rank() << "/" << size() << "]" << ": received invalid Kill message" << std::endl;
         }
         break;
      }

      case message::Message::CREATEPORT: {

         const message::CreatePort *cp =
            static_cast<const message::CreatePort *>(message);
         Port *port = cp->getPort();
         std::string name = port->getName();
         std::string::size_type p = name.find('[');
         std::string basename = name;
         size_t idx = 0;
         if (p != std::string::npos) {
            basename = name.substr(0, p-1);
            std::stringstream idxstr(name.substr(p+1));
            idxstr >> idx;
         }
         Port *existing = NULL;
         Port *parent = NULL;
         Port *newport = NULL;
         switch (port->getType()) {
            case Port::INPUT:
               existing = findInputPort(name);
               if (!existing)
                  parent = findInputPort(basename);
               if (parent) {
                  newport = parent->child(idx, true);
                  inputPorts[name] = newport;
               }
               break;
            case Port::OUTPUT:
               existing = findOutputPort(name);
               if (!existing)
                  parent = findInputPort(basename);
               if (parent) {
                  newport = parent->child(idx, true);
                  outputPorts[name] = newport;
               }
               break;
            case Port::PARAMETER:
               // does not really happen - handled in AddParameter
               break;
            case Port::ANY:
               break;
         }
         if (newport) {
            message::CreatePort np(newport);
            np.setUuid(cp->uuid());
            sendMessage(np);
            const Port::PortSet &links = newport->linkedPorts();
            for (Port::PortSet::iterator it = links.begin();
                  it != links.end();
                  ++it) {
               message::CreatePort linked(*it);
               linked.setUuid(cp->uuid());
               sendMessage(linked);
            }
         }
         break;
      }

      case message::Message::CONNECT: {

         const message::Connect *conn =
            static_cast<const message::Connect *>(message);
         Port *port = NULL;
         Port *other = NULL;
         const Port::PortSet *ports = NULL;
         std::string ownPortName;
         //std::cerr << name() << " receiving connection: " << conn->getModuleA() << ":" << conn->getPortAName() << " -> " << conn->getModuleB() << ":" << conn->getPortBName() << std::endl;
         if (conn->getModuleA() == id()) {
            port = findOutputPort(conn->getPortAName());
            ownPortName = conn->getPortAName();
            if (port) {
               other = new Port(conn->getModuleB(), conn->getPortBName(), Port::INPUT);
               ports = &port->connections();
            }
            
         } else if (conn->getModuleB() == id()) {
            ownPortName = conn->getPortBName();
            port = findInputPort(conn->getPortBName());
            if (port) {
               other = new Port(conn->getModuleA(), conn->getPortAName(), Port::OUTPUT);
               ports = &port->connections();
            }
         }

         if (ports && port && other) {
            if (ports->find(other) == ports->end())
               port->addConnection(other);
            else
               delete other;
         } else {
            if (!findParameter(ownPortName))
               std::cerr << name() << " did not find port" << std::endl;
         }
         break;
      }

      case message::Message::DISCONNECT: {

         const message::Disconnect *disc =
            static_cast<const message::Disconnect *>(message);
         Port *port = NULL;
         Port *other = NULL;
         const Port::PortSet *ports = NULL;
         if (disc->getModuleA() == id()) {
            port = findOutputPort(disc->getPortAName());
            if (port) {
               other = new Port(disc->getModuleB(), disc->getPortBName(), Port::INPUT);
               ports = &port->connections();
            }
            
         } else if (disc->getModuleB() == id()) {
            port = findInputPort(disc->getPortBName());
            if (port) {
               other = new Port(disc->getModuleA(), disc->getPortAName(), Port::OUTPUT);
               ports = &port->connections();
            }
         }

         if (ports && port && other) {
            Port *p = port->removeConnection(other);
            delete other;
            delete p;
         }
         break;
      }

      case message::Message::COMPUTE: {

         const message::Compute *comp =
            static_cast<const message::Compute *>(message);
         if (m_executionCount < comp->getExecutionCount())
            m_executionCount = comp->getExecutionCount();

         if (comp->getExecutionCount() > 0) {
            // Compute not triggered by adding an object
            for (std::map<std::string, Port *>::iterator pit = inputPorts.begin();
                  pit != inputPorts.end();
                  ++pit) {
               pit->second->objects() = m_cache.getObjects(pit->first);
            }
         }
         /*
         std::cerr << "    module [" << name() << "] [" << id() << "] ["
                   << rank() << "/" << size << "] compute" << std::endl;
         */
         message::Busy busy;
         busy.setUuid(comp->uuid());
         sendMessage(busy);
         bool ret = false;
         try {
            ret = compute();
         } catch (std::exception &e) {
            std::cerr << name() << ": exception - " << e.what() << std::endl;
         }
         message::Idle idle;
         idle.setUuid(comp->uuid());
         sendMessage(idle);
         return ret;
         break;
      }

      case message::Message::ADDOBJECT: {

         const message::AddObject *add =
            static_cast<const message::AddObject *>(message);
         addInputObject(add->getPortName(), add->takeObject());
         break;
      }

      case message::Message::SETPARAMETER: {

         const message::SetParameter *param =
            static_cast<const message::SetParameter *>(message);

         if (param->getModule() == id()) {

            // sent by controller
            switch (param->getParameterType()) {
               case Parameter::Integer:
                  setIntParameter(param->getName(), param->getInteger(), param);
                  break;
               case Parameter::Float:
                  setFloatParameter(param->getName(), param->getFloat(), param);
                  break;
               case Parameter::Vector:
                  setVectorParameter(param->getName(), param->getVector(), param);
                  break;
               case Parameter::String:
                  setStringParameter(param->getName(), param->getString(), param);
                  break;
               default:
                  std::cerr << "Module::handleMessage: unknown parameter type " << param->getParameterType() << std::endl;
                  assert("unknown parameter type" == 0);
                  break;
            }

            if (Parameter *p = findParameter(param->getName())) {
               parameterChanged(p);
            }

            // notify controller about current value
#if 0
            if (Parameter *p = findParameter(param->getName())) {
               sendMessage(message::SetParameter(id(), param->getName(), p));
            }
#endif
         } else {

            parameterChanged(param->senderId(), param->getName(), *param);
         }
         break;
      }

      case message::Message::SETPARAMETERCHOICES: {
         const message::SetParameterChoices *choices = static_cast<const message::SetParameterChoices *>(message);
         if (choices->senderId() != id()) {
            //FIXME: handle somehow
            //parameterChanged(choices->senderId(), choices->getName());
         }
         break;
      }

      case message::Message::ADDPARAMETER: {

         const message::AddParameter *param =
            static_cast<const message::AddParameter *>(message);

         parameterAdded(param->senderId(), param->getName(), *param, param->moduleName());
         break;
      }

      case message::Message::SPAWN: {
         const message::Spawn *spawn =
            static_cast<const message::Spawn *>(message);
         m_otherModuleMap[spawn->spawnId()] = spawn->getName();
         break;
      }

      case message::Message::MODULEEXIT: {
         const message::ModuleExit *exit =
            static_cast<const message::ModuleExit *>(message);
         m_otherModuleMap.erase(exit->senderId());
         break;
      }

      case message::Message::BARRIER: {

         const message::Barrier *barrier = static_cast<const message::Barrier *>(message);
         message::BarrierReached reached(barrier->getBarrierId());
         reached.setUuid(barrier->uuid());
         sendMessage(reached);
         break;
      }

      case message::Message::OBJECTRECEIVED:
         // currently only relevant for renderers
         break;

      default:
         std::cerr << "    module [" << name() << "] [" << id() << "] ["
                   << rank() << "/" << size() << "] unknown message type ["
                   << message->type() << "]" << std::endl;

         break;
   }

   return true;
}

std::string Module::getModuleName(int id) const {

   auto it = m_otherModuleMap.find(id);
   if (it == m_otherModuleMap.end())
      return std::string();

   return it->second;
}

Module::~Module() {

   m_cache.clear();

   if (m_origStreambuf)
      std::cerr.rdbuf(m_origStreambuf);
   delete m_streambuf;

   vistle::message::ModuleExit m;
   sendMessage(m);

   Shm::the().detach();

   std::cerr << "  module [" << name() << "] [" << id() << "] [" << rank()
             << "/" << size() << "] quit" << std::endl;

   MPI_Barrier(MPI_COMM_WORLD);
}

namespace {

using namespace boost;

class InstModule: public Module {
public:
   InstModule()
   : Module("inst", "anything", 0, 1, 1)
   {}
   bool compute() { return true; }
};

struct instantiator {
   template<typename P> P operator()(P) {
      InstModule m;
      ParameterBase<P> p(m.id(), "param");
      m.setParameterMinimum(&p, P());
      m.setParameterMaximum(&p, P());
      return P();
   }
};
}

void instantiate_parameter_functions() {

   mpl::for_each<Parameters>(instantiator());
}

void vistle::Module::sendInfo(const char *fmt, ...)
{
   if(!fmt) {
      fmt = "(empty message)";
   }
   char *text = new char[strlen(fmt)+500];

   va_list args;
   va_start(args, fmt);
   vsnprintf(text, strlen(fmt)+500, fmt, args);
   va_end(args);

   message::SendInfo info(text);
   sendMessage(info);
}

void Module::sendInfo(const message::Message &msg, const char *fmt, ...)
{
   if(!fmt) {
      fmt = "(empty message)";
   }
   char *text = new char[strlen(fmt)+500];

   va_list args;
   va_start(args, fmt);
   vsnprintf(text, strlen(fmt)+500, fmt, args);
   va_end(args);

   message::SendInfo info(text, &msg);
   sendMessage(info);
}

} // namespace vistle
