#ifndef TCPMESSAGE_H
#define TCPMESSAGE_H

#include <boost/asio/ip/tcp.hpp>

#include "export.h"

namespace vistle {

namespace message {

class Message;

typedef boost::asio::ip::tcp::socket socket_t;

bool V_COREEXPORT send(socket_t &sock, const message::Message &msg);
void V_COREEXPORT async_send(socket_t &sock, const message::Message &msg, const std::function<void(boost::system::error_code ec)> &handler);

bool V_COREEXPORT recv(socket_t &sock, message::Message &msg, bool &received, bool block=false);
void V_COREEXPORT async_recv(socket_t &sock, const message::Message &msg, const std::function<void(boost::system::error_code ec)> &handler);

} // namespace message
} // namespace vistle
#endif
