#ifndef INSITU_MESSAGE_H
#define INSITU_MESSAGE_H

#include "export.h"

#include <string>
#include <array>

#include <boost/asio/ip/tcp.hpp>
#include <boost/asio/streambuf.hpp>

#include <boost/mpi.hpp>
#include <core/message.h>
#include <core/messages.h>
#include <core/messagepayload.h>
#include <core/tcpmessage.h>
#include <core/archives_config.h>
#include <core/archives.h>

#include <util/vecstreambuf.h>

#include <boost/interprocess/shared_memory_object.hpp>
#include <boost/interprocess/mapped_region.hpp>
#include <boost/interprocess/sync/scoped_lock.hpp>


namespace insitu{
namespace message {

enum class V_INSITUMESSAGEEXPORT InSituMessageType {
    Invalid
    , ShmInit
    , AddObject
    , SetPorts //detected ports from sim to module <--> connected ports from module to Engine
    , SetCommands
    , Ready
    , ExecuteCommand
    , GoOn
    , ConstGrids
    , NthTimestep
    , ConnectionClosed
    , VTKVariables
    , CombineGrids
    , KeepTimesteps

};


struct V_INSITUMESSAGEEXPORT InSituMessageBase {
    InSituMessageBase(InSituMessageType t) :m_type(t) {};
    InSituMessageType type() const;
protected:

    InSituMessageType m_type;

};
#define COMMA ,

#define DECLARE_ENGINE_MESSAGE_WITH_PARAM(messageType,  payloadType)\
struct V_INSITUMESSAGEEXPORT messageType : public InSituMessageBase\
{\
    typedef payloadType value_type;\
    friend class InSituTcpMessage; \
    messageType(const payloadType& payloadName):InSituMessageBase(type), value(payloadName){}\
    static const InSituMessageType type = InSituMessageType::messageType;\
    payloadType value; \
    ARCHIVE_ACCESS\
    template<class Archive>\
    void serialize(Archive& ar) {\
       ar& value;\
    }\
private:\
    messageType():InSituMessageBase(type){}\
};\



#define DECLARE_ENGINE_MESSAGE(messageType)\
struct V_INSITUMESSAGEEXPORT messageType : public InSituMessageBase {\
    static const InSituMessageType type = InSituMessageType::messageType;\
    messageType() :InSituMessageBase(type) {}\
    ARCHIVE_ACCESS\
    template<class Archive>\
    void serialize(Archive& ar) {}\
};\

DECLARE_ENGINE_MESSAGE(Invalid)
DECLARE_ENGINE_MESSAGE(GoOn)
DECLARE_ENGINE_MESSAGE(ConnectionClosed)

DECLARE_ENGINE_MESSAGE_WITH_PARAM(ShmInit,  std::string)
DECLARE_ENGINE_MESSAGE_WITH_PARAM(AddObject,  std::string)
DECLARE_ENGINE_MESSAGE_WITH_PARAM(SetPorts,  std::vector<std::vector<std::string>>)
DECLARE_ENGINE_MESSAGE_WITH_PARAM(SetCommands,  std::vector<std::string>)
DECLARE_ENGINE_MESSAGE_WITH_PARAM(Ready,  bool)
DECLARE_ENGINE_MESSAGE_WITH_PARAM(ExecuteCommand,  std::string)
DECLARE_ENGINE_MESSAGE_WITH_PARAM(ConstGrids,  bool)
DECLARE_ENGINE_MESSAGE_WITH_PARAM(CombineGrids,  bool)
DECLARE_ENGINE_MESSAGE_WITH_PARAM(NthTimestep,  size_t)
DECLARE_ENGINE_MESSAGE_WITH_PARAM(VTKVariables,  bool)
DECLARE_ENGINE_MESSAGE_WITH_PARAM(KeepTimesteps, bool)


struct V_INSITUMESSAGEEXPORT InSituMessage : public vistle::message::MessageBase<InSituMessage, vistle::message::INSITU> {
    InSituMessage(InSituMessageType t) :m_ismType(t) {}
    InSituMessageType ismType() const {
        return m_ismType;
    }
private:
    InSituMessageType m_ismType;
};
static_assert(sizeof(InSituMessage) <= vistle::message::Message::MESSAGE_SIZE, "message too large");




//after being initialized, sends and receives messages in a blocking manner. 
//When the connection is closed returns EngineMEssageType::ConnectionClosed and becomes uninitialized.
//while uninitialized calls to send and received are ignored.
//Received Messages are broadcasted to all ranks so make sure they all call receive together.
class V_INSITUMESSAGEEXPORT InSituTcpMessage {
public:

    InSituMessageType type() const;

    static void initialize(std::shared_ptr< boost::asio::ip::tcp::socket> socket, boost::mpi::communicator comm);
    static bool isInitialized();
    template<typename SomeMessage>
    static bool send(const SomeMessage& msg) {
        if (!m_initialized) {
            std::cerr << "InSituTcpMessage uninitialized: can not send message!" << std::endl;
            return false;
        }
        if (msg.type == InSituMessageType::Invalid) {
            std::cerr << "InSituTcpMessage : can not send invalid message!" << std::endl;
            return false;
        }
        bool error = false;
        if (m_comm.rank() == 0) {
            boost::system::error_code err;
            InSituMessage ism(msg.type);
            vistle::vecostreambuf<vistle::buffer> buf;
            vistle::oarchive ar(buf);
            ar& msg;
            vistle::buffer vec = buf.get_vector();

            ism.setPayloadSize(vec.size());
            vistle::message::send(*m_socket, ism, err, &vec);
            if (err) {
                error = true;
            }
        }
        return error;
    }

    static InSituTcpMessage recv();

    template<typename SomeMessage>
    SomeMessage& unpackOrCast() {


        assert(SomeMessage::type != type());
        if (!m_msg) {
            vistle::vecistreambuf<vistle::buffer> buf(m_payload);
            m_msg.reset(new SomeMessage{});
            try {
                vistle::iarchive ar(buf);
                ar&* static_cast<SomeMessage*>(m_msg.get());
            } catch (yas::io_exception & ex) {
                std::cerr << "ERROR: failed to get InSituTcpMessage payload from " << m_payload.size() << " bytes: " << ex.what() << std::endl;
            }
        }
        return *static_cast<SomeMessage*>(m_msg.get());

    }



private:
    InSituTcpMessage(InSituMessageType type, vistle::buffer&& payload);
    InSituTcpMessage();
    InSituMessageType m_type;
    vistle::buffer m_payload;
    std::unique_ptr<InSituMessageBase> m_msg;

    static bool m_initialized;
    static boost::mpi::communicator m_comm;
    static std::shared_ptr< boost::asio::ip::tcp::socket> m_socket;




};

class V_INSITUMESSAGEEXPORT SyncShmIDs {
public:
    enum class Mode {
          Create
        , Attach
    };
#ifdef MODULE_THREAD
    typedef std::lock_guard<std::mutex> Guard;
#else
    typedef boost::interprocess::scoped_lock<boost::interprocess::interprocess_mutex> Guard;
#endif

    static void initialize(int moduleID, int rank, int instance, Mode mode);
    static void close();
    static bool isInitialized();
    //we might need to use a mutex for this?
    static void set(int objID, int arrayID); 
    static int objectID(); 
    static int arrayID(); 

    struct ShmData {
        int objID = -1;
        int arrayID = -1;
#ifdef MODULE_THREAD
        std::mutex mutex;
#else
        boost::interprocess::interprocess_mutex mutex;
#endif

    };
private:
#ifndef MODULE_THREAD
    class ShmSegment{
    public:
        ShmSegment() {}
        ShmSegment(ShmSegment&& other);
        ShmSegment(const std::string& name, Mode mode) throw(vistle::exception);
        ~ShmSegment();
        const ShmData* data() const;
        ShmData* data();
        ShmSegment& operator=(ShmSegment&& other) noexcept;
    private:
        std::string m_name;

        boost::interprocess::mapped_region m_region;
    };
    static ShmSegment m_segment;
#else
    class Data { //wrapper if we dont need shm
public:
    ShmData* data(){
        return &m_data;
    }
    private:
        ShmData m_data;
};
    static Data m_segment;
#endif
    static int m_rank;
    static int m_moduleID;

    int m_objectID;
    int m_array_ID;

    static bool m_initialized;
};
} //message
} //insitu



#endif // !ENGINE_MESSAGE_H
