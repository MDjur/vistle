#include <cstdio>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <string>

#include <util/hostname.h>
#include <util/sleep.h>
#include <core/message.h>
#include <core/tcpmessage.h>
#include <core/parameter.h>
#include <core/port.h>
#include "userinterface.h"

namespace asio = boost::asio;

namespace vistle {

UserInterface::UserInterface(const std::string &host, const unsigned short port, StateObserver *observer)
: m_id(-1)
, m_remoteHost(host)
, m_remotePort(port)
, m_isConnected(false)
, m_stateTracker("UI state")
, m_socket(m_ioService)
, m_locked(false)
{
   message::DefaultSender::init(message::Id::UI, 0);

   if (observer)
      m_stateTracker.registerObserver(observer);

   std::cerr << "  userinterface ["  << id() << "] started as " << hostname() << ":"
#ifndef _WIN32
             << getpid() << std::endl;
#else
             << std::endl;
#endif

   m_hostname = hostname();

   tryConnect();
}

void UserInterface::stop() {

   if (isConnected()) {
      vistle::message::ModuleExit m;
      m.setDestId(vistle::message::Id::LocalHub);
      sendMessage(m);
   }

   cancel();
}

void UserInterface::cancel() {

   m_socket.cancel();
   if (isConnected()) {
      m_socket.shutdown(asio::ip::tcp::socket::shutdown_both);
   }
   m_ioService.stop();
}

int UserInterface::id() const {

   return m_id;
}

const std::string &UserInterface::host() const {

   return m_hostname;
}

boost::asio::ip::tcp::socket &UserInterface::socket() {

   return m_socket;
}

const boost::asio::ip::tcp::socket &UserInterface::socket() const {

   return m_socket;
}

bool UserInterface::tryConnect() {

   assert(!isConnected());
   std::string host = m_remoteHost;
   if (host == m_hostname)
       host = "localhost";

   asio::ip::tcp::resolver resolver(m_ioService);
   asio::ip::tcp::resolver::query query(host, boost::lexical_cast<std::string>(m_remotePort));
   asio::ip::tcp::resolver::iterator endpoint_iterator = resolver.resolve(query);
   boost::system::error_code ec;
   asio::connect(socket(), endpoint_iterator, ec);
   if (ec) {
      std::cerr << "  userinterface [" << id() << "]: could not establish connection to "
      << host << ":" << m_remotePort << std::endl;
      m_isConnected = false;
      return false;
   }
   m_isConnected = true;
   return true;
}

bool UserInterface::isConnected() const {

   return m_isConnected && socket().is_open();
}

StateTracker &UserInterface::state() {

   return m_stateTracker;
}

bool UserInterface::dispatch() {

   bool work = false;
   while (isConnected()) {

      work = true;
      bool received = false;

      message::Buffer buf;
      std::vector<char> payload;
      if (!message::recv(socket(), buf, received, true /* blocking */, &payload)) {
         return false;
      }
      
      if (!received)
         break;

      if (!handleMessage(&buf, payload))
         return false;

   }

   vistle::adaptive_wait(work, this);

   return true;
}


bool UserInterface::sendMessage(const message::Message &message, const std::vector<char> *payload) {

   if (m_locked) {
      assert(!payload);
      m_sendQueue.emplace_back(message);
      return true;
   }

   return message::send(socket(), message, payload);
}


bool UserInterface::handleMessage(const vistle::message::Message *message, const std::vector<char> &payload) {

   bool ret = m_stateTracker.handle(*message);

   {
      std::lock_guard<std::mutex> lock(m_messageMutex);
      MessageMap::iterator it = m_messageMap.find(const_cast<message::uuid_t &>(message->uuid()));
      if (it != m_messageMap.end()) {
         it->second->buf.resize(message->size());
         memcpy(it->second->buf.data(), message, message->size());
         it->second->received = true;
         it->second->cond.notify_all();
      }
   }

   switch (message->type()) {
      case message::IDENTIFY: {
         const message::Identify *id = static_cast<const message::Identify *>(message);
         if (id->identity() == message::Identify::REQUEST) {
            const message::Identify reply(message::Identify::UI);
            sendMessage(reply);
         }
         return true;
         break;
      }

      case message::SETID: {
         const message::SetId *id = static_cast<const message::SetId *>(message);
         m_id = id->getId();
         assert(m_id > 0);
         message::DefaultSender::init(id->senderId(), -m_id);
         //std::cerr << "received new UI id: " << m_id << std::endl;
         break;
      }

      case message::LOCKUI: {
         auto lock = static_cast<const message::LockUi *>(message);
         m_locked = lock->locked();
         if (!m_locked) {
            for (auto &m: m_sendQueue) {
               sendMessage(m);
            }
            m_sendQueue.clear();
         }
         break;
      }

      case message::QUIT: {
         const message::Quit *quit = static_cast<const message::Quit *>(message);
         (void)quit;
         return false;
         break;
      }

      case message::FILEQUERYRESULT: {
          auto &fq = message->as<message::FileQueryResult>();
          bool found = false;
          for (auto b: m_fileBrowser) {
              if (b->id() == fq.filebrowserId()) {
                  b->handleMessage(fq, payload);
                  found = true;
              }
              if (!found) {
                  std::cerr << "userinterface [" << host() << "] [" << id() << "] "
                            << "did not find filebrowser with id " << fq.filebrowserId() << std::endl;
              }
          }
          break;
      }

      default:
         break;
   }

   return ret;
}

bool UserInterface::getLockForMessage(const message::uuid_t &uuid) {

   std::lock_guard<std::mutex> lock(m_messageMutex);
   MessageMap::iterator it = m_messageMap.find(uuid);
   if (it == m_messageMap.end()) {
      it = m_messageMap.insert(std::make_pair(uuid, std::shared_ptr<RequestedMessage>(new RequestedMessage()))).first;
   }
   it->second->mutex.lock();
   //m_messageMap[const_cast<message::uuid_t &>(uuid)]->mutex.lock();
   return true;
}

bool UserInterface::getMessage(const message::uuid_t &uuid, message::Message &msg) {

   m_messageMutex.lock();
   MessageMap::iterator it = m_messageMap.find(uuid);
   if (it == m_messageMap.end()) {
      m_messageMutex.unlock();
      return false;
   }

   if (!it->second->received) {
      std::mutex &mutex = it->second->mutex;
      std::condition_variable &cond = it->second->cond;
      std::unique_lock<std::mutex> lock(mutex, std::adopt_lock_t());

      m_messageMutex.unlock();
      cond.wait(lock);
      m_messageMutex.lock();
   }

   if (!it->second->received) {
      m_messageMutex.unlock();
      return false;
   }

   memcpy(&msg, &*it->second->buf.data(), it->second->buf.size());
   m_messageMap.erase(it);
   m_messageMutex.unlock();
   return true;
}

void UserInterface::registerObserver(StateObserver *observer) {

    m_stateTracker.registerObserver(observer);
}

void UserInterface::registerFileBrowser(FileBrowser *browser)
{
    assert(browser->m_ui == nullptr);

    ++m_fileBrowserCount;
    browser->m_ui = this;
    browser->m_id = m_fileBrowserCount;
    m_fileBrowser.push_back(browser);
}

void UserInterface::removeFileBrowser(FileBrowser *browser)
{
    assert(browser->m_ui == this);

    auto it = std::find(m_fileBrowser.begin(), m_fileBrowser.end(), browser);
    if (it != m_fileBrowser.end()) {
        m_fileBrowser.erase(it);
        browser->m_ui = nullptr;
    }
}

UserInterface::~UserInterface() {

   std::cerr << "  userinterface [" << host() << "] [" << id()
             << "] quit" << std::endl;
}

FileBrowser::~FileBrowser() {
}

int FileBrowser::id() const {
    return m_id;
}

bool FileBrowser::sendMessage(const message::Message &message, const std::vector<char> *payload)
{
    assert(m_ui);
    message::Buffer buf(message);
    if (buf.type() == message::FILEQUERY) {
        auto &fq = buf.as<message::FileQuery>();
        fq.setFilebrowserId(m_id);
    }
    return m_ui->sendMessage(buf, payload);
}

} // namespace vistle
