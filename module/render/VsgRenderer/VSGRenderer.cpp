#include "VSGRenderer.h"
#include "vistle/core/message.h"
#include "vistle/core/messages.h"

//std
#include <iostream>
#include <memory>
#include <vector>
#include <vulkan/vulkan_core.h>

//vsg
/* #include <vsg/core/Object.h> */
/* #include <vsg/app/Trackball.h> */
/* #include <vsg/app/Viewer.h> */
/* #include <vsg/app/WindowTraits.h> */
/* #include <vsg/io/FileSystem.h> */

#ifdef vsgXchange_FOUND

#include <vsgXchange/all.h>

#endif

#define CERR std::cerr << "VsgRenderer: "

namespace {

void printVSGMetaData(vsg::ref_ptr<vsg::Object> obj)
{
    CERR << "Show Metadata";

    if (auto auxiliary = obj->getAuxiliary()) {
        for (auto &[key, object]: auxiliary->userObjects) {
            if (auto s = dynamic_cast<vsg::stringValue *>(object.get()))
                CERR << "metadata key = " << key << ", stringValue = " << s->value() << std::endl;
            else if (auto d = dynamic_cast<vsg::doubleValue *>(object.get()))
                CERR << "metadata key = " << key << ", doubleValue = " << d->value() << std::endl;
            else
                CERR << "metadata key = " << key << ", object = " << object << std::endl;
        }
    } else {
        CERR << "No vsg::Auxiliary assigned to Object." << std::endl;
    }
}

vsg::ref_ptr<vsg::Camera> createCameraForScene(vsg::Node *scene, int32_t x, int32_t y, uint32_t width, uint32_t height)
{
    // compute the bounds of the scene graph to help position camera
    vsg::ComputeBounds computeBounds;
    scene->accept(computeBounds);
    vsg::dvec3 centre = (computeBounds.bounds.min + computeBounds.bounds.max) * 0.5;
    double radius = vsg::length(computeBounds.bounds.max - computeBounds.bounds.min) * 0.6;
    double nearFarRatio = 0.001;

    // set up the camera
    auto lookAt = vsg::LookAt::create(centre + vsg::dvec3(0.0, -radius * 3.5, 0.0), centre, vsg::dvec3(0.0, 0.0, 1.0));

    auto perspective = vsg::Perspective::create(30.0, static_cast<double>(width) / static_cast<double>(height),
                                                nearFarRatio * radius, radius * 4.5);

    auto viewportstate = vsg::ViewportState::create(x, y, width, height);

    return vsg::Camera::create(perspective, lookAt, viewportstate);
}

} // namespace

VSGRenderer::VSGRenderer(const std::string &name, int moduleID, mpi::communicator comm)
: Renderer(name, moduleID, comm), m_renderManager(this), m_asyncFrames(0)
{
    // create options object that is used to guide IO operations
    auto options = vsg::Options::create();
    options->fileCache = vsg::getEnv("VSG_FILE_CACHE");
    options->paths = vsg::getEnvPaths("VSG_FILE_PATH");

#ifdef vsgXchange_all
    options->add(vsgXchange::all::create());
#endif

    m_scenegraph = vsg::Group::create();
    m_scenegraph->addChild(
        vsg::read_cast<vsg::Node>(vsg::Path("/home/hpcmdjur/git/vsgExamples/data/models/teapot.vsgt")));
    printVSGMetaData(m_scenegraph);
    m_viewer = vsg::Viewer::create();
    auto windowTraits = vsg::WindowTraits::create();
    windowTraits->windowTitle = "VsgRenderer";
    //Vsync
    windowTraits->swapchainPreferences.presentMode = VK_PRESENT_MODE_FIFO_KHR;

    auto window = vsg::Window::create(windowTraits);
    m_viewer->addWindow(window);

    // set up the camera
    auto camera = createCameraForScene(m_scenegraph, 0, 0, window->extent2D().width, window->extent2D().height);
    auto main_view = vsg::View::create(camera, m_scenegraph);

    // add close handler to respond to pressing the window close window button and pressing escape
    m_viewer->addEventHandler(vsg::CloseHandler::create(m_viewer));

    // add a trackball event handler to control the camera view use the mouse
    auto main_trackball = vsg::Trackball::create(camera);
    main_trackball->addWindow(window);
    m_viewer->addEventHandler(main_trackball);

    // create a command graph to render the scene on specified window
    auto commandGraph = vsg::createCommandGraphForView(window, camera, m_scenegraph);
    m_viewer->assignRecordAndSubmitTaskAndPresentation({commandGraph});
    m_viewer->setupThreading();

    // compile all the Vulkan objects and transfer data required to render the scene
    m_viewer->compile();
}

VSGRenderer::~VSGRenderer()
{
    CERR << "destroying" << std::endl;
    /* prepareQuit(); */
}

/* void VSGRenderer::flush() */
/* { */
/*     CERR << "flushing outstanding frames..." << std::endl; */
/*     for (size_t f = m_asyncFrames; f > 0; --f) { */
/*         composite(f - 1); */
/*     } */
/*     assert(m_numFramesToComposite == 0); */
/* } */

void VSGRenderer::prepareQuit()
{
    removeAllObjects();
    Renderer::prepareQuit();
    m_viewer->stopThreading();
    m_viewer->close();
    /* sendMessage(vistle::message::Quit()); */
}

bool VSGRenderer::composite(size_t maxQueued)
{
    return true;
}

bool VSGRenderer::render()
{
    if (!m_viewer->active()) {
        prepareQuit();
        return false;
    }
    if (!m_viewer->advanceToNextFrame())
        return false;

    m_viewer->handleEvents();
    m_viewer->recordAndSubmit();
    m_viewer->update();
    m_viewer->present();

    return true;
}

std::shared_ptr<vistle::RenderObject> VSGRenderer::addObject(int senderId, const std::string &senderPort,
                                                             vistle::Object::const_ptr container,
                                                             vistle::Object::const_ptr geometry,
                                                             vistle::Object::const_ptr normals,
                                                             vistle::Object::const_ptr texture)
{
    return std::make_shared<vistle::RenderObject>(senderId, senderPort, container, geometry, normals, texture);
}

void VSGRenderer::removeObject(std::shared_ptr<vistle::RenderObject> ro)
{
    /* auto oro = std::static_pointer_cast<OsgRenderObject>(ro); */
    /* timesteps->removeObject(oro->node, oro->timestep); */
    m_renderManager.removeObject(ro);
}

void VSGRenderer::connectionAdded(const vistle::Port *from, const vistle::Port *to)
{
    CERR << "new connection from " << *from << " to " << *to << std::endl;
    Renderer::connectionAdded(from, to);
    if (from == m_renderManager.outputPort()) {
        m_renderManager.connectionAdded(to);
    }
}

void VSGRenderer::connectionRemoved(const vistle::Port *from, const vistle::Port *to)
{
    if (from == m_renderManager.outputPort()) {
        m_renderManager.connectionRemoved(to);
    }
    Renderer::connectionRemoved(from, to);
}

bool VSGRenderer::handleMessage(const vistle::message::Message *message, const vistle::MessagePayload &payload)
{
    if (m_renderManager.handleMessage(message, payload))
        return true;

    return Renderer::handleMessage(message, payload);
}


MODULE_MAIN(VSGRenderer)
