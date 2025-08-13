#include "VSGRenderer.h"
#include "VSGRenderObject.h"
#include "VSGTimestepHandler.h"
#include "utils/visitors.h"

//vistle
#include "vistle/core/database.h"
#include "vistle/core/message.h"
#include "vistle/core/messages.h"
#include "vistle/core/object.h"
#include <vistle/core/points.h>
#include "vistle/module/module.h"
#include "vistle/renderer/renderer.h"
#include "vsg/core/ref_ptr.h"
#include "vsg/utils/ComputeBounds.h"
#include <vistle/core/placeholder.h>

//std
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

//vulkan
#include <vulkan/vulkan_core.h>

#ifdef vsgXchange_FOUND

#include <vsgXchange/all.h>

#endif

#define CERR std::cerr << "VsgRenderer: "

#define DEBUG 1

using namespace vistle;
namespace {

/**
 * @brief Print metadata of a vsg::Object add via obj->setValue("value", value).
 * @param obj the vsg::Object to print metadata from
 * @return string with metadata
 */
std::string printVSGMetaData(vsg::ref_ptr<vsg::Object> obj)
{
    std::stringstream ss;
    ss << "Show Metadata";

    if (auto auxiliary = obj->getAuxiliary()) {
        for (auto &[key, object]: auxiliary->userObjects) {
            if (auto s = dynamic_cast<vsg::stringValue *>(object.get()))
                ss << "metadata key = " << key << ", stringValue = " << s->value() << std::endl;
            else if (auto d = dynamic_cast<vsg::doubleValue *>(object.get()))
                ss << "metadata key = " << key << ", doubleValue = " << d->value() << std::endl;
            else
                ss << "metadata key = " << key << ", object = " << object << std::endl;
        }
    } else {
        ss << "No vsg::Auxiliary assigned to Object." << std::endl;
    }
    return ss.str();
}

/**
 * @brief Create a camera to view the scene graph.
 * @param scene the scenegraph to view
 * @param x the x position of the lower left corner of the view port
 * @param y the y position of the lower left corner of the view port
 * @param width the width of the view port
 * @param height the height of the view port
 * @return the camera to view the scene
 */
vsg::ref_ptr<vsg::Camera> createCameraForScene(vsg::Node *scene, int32_t x, int32_t y, uint32_t width, uint32_t height)
{
    // compute the bounds of the scene graph to help position camera
    vsg::ComputeBounds computeBounds;
    scene->accept(computeBounds);
    vsg::dvec3 centre = (computeBounds.bounds.min + computeBounds.bounds.max) * 0.5;
    /* double radius = vsg::length(computeBounds.bounds.max - computeBounds.bounds.min) * 0.6; */
    double radius = vsg::length(computeBounds.bounds.max - computeBounds.bounds.min) * 5;
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
: Renderer(name, moduleID, comm), m_renderManager(this) //, m_asyncFrames(0)
{
    // create options object that is used to guide IO operations
    auto options = vsg::Options::create();
    options->fileCache = vsg::getEnv("VSG_FILE_CACHE");
    options->paths = vsg::getEnvPaths("VSG_FILE_PATH");

#ifdef vsgXchange_all
    options->add(vsgXchange::all::create());
#endif

#if DEBUG
    std::stringstream strstream;
#endif // DEBUG

    m_viewer = vsg::Viewer::create();
    auto windowTraits = vsg::WindowTraits::create();
    windowTraits->windowTitle = "VsgRenderer";

    //Vsync => VULKAN presentation mode (Immediate, FIFO, FIFO_Relaxed, Mailbox)
    windowTraits->swapchainPreferences.presentMode = VK_PRESENT_MODE_FIFO_RELAXED_KHR;

    // number of images rendered parallel in swapchain (3 triple buffer, 2 double buffer)
    windowTraits->swapchainPreferences.imageCount = 3;

    /* // enable wireframe polygon mode */
    auto deviceFeatures = vsg::DeviceFeatures::create();
    deviceFeatures->get().fillModeNonSolid = VK_TRUE;
    deviceFeatures->get().wideLines = VK_TRUE;
    windowTraits->deviceFeatures = deviceFeatures;

    /* // disables titlebar */
    /* windowTraits->decoration = false; */

    auto window = vsg::Window::create(windowTraits);
    m_viewer->addWindow(window);

    // enable vulkan validationlayer
    bool debugLayer = true;
    bool lunargApiDumpLayer = false;
    /* bool enableGeometryShader = true; */

    // enable multisampling (anti-aliasing) => VK_SAMPLE_COUNT_1_BIT = default
    VkSampleCountFlagBits samples = VK_SAMPLE_COUNT_8_BIT;
    /* VkSampleCountFlagBits samples = VK_SAMPLE_COUNT_32_BIT; */
    windowTraits->samples = samples;
    /* VkSampleCountFlagBits samples = VK_SAMPLE_COUNT_1_BIT; */

    uint32_t vulkanVersion = VK_API_VERSION_1_0;
    if (samples != VK_SAMPLE_COUNT_1_BIT)
        vulkanVersion = VK_API_VERSION_1_2;

    // create instance
    vsg::Names instanceExtensions;
    vsg::Names requestedLayers;
    if (debugLayer || lunargApiDumpLayer) {
        instanceExtensions.push_back(VK_KHR_GET_PHYSICAL_DEVICE_PROPERTIES_2_EXTENSION_NAME);
        instanceExtensions.push_back(VK_EXT_DEBUG_REPORT_EXTENSION_NAME);
        requestedLayers.push_back("VK_LAYER_KHRONOS_validation");
        sendInfo("Vulkan validation layer enabled");
        if (lunargApiDumpLayer)
            requestedLayers.push_back("VK_LAYER_LUNARG_api_dump");
    }

    vsg::Names validatedNames = vsg::validateInstancelayerNames(requestedLayers);
    auto instance = vsg::Instance::create(instanceExtensions, validatedNames, vulkanVersion);
    /* auto [physicalDevice, queueFamily] = instance->getPhysicalDeviceAndQueueFamily(VK_QUEUE_GRAPHICS_BIT); */
    /* auto physicalDevice = instance->getPhysicalDevice(VK_QUEUE_GRAPHICS_BIT); */
    //some debug info about the physical device
#if DEBUG
    auto physicalDevice = window->getOrCreatePhysicalDevice();
    strstream << "physicalDevice = " << physicalDevice << std::endl;
    for (auto &queueFamilyProperties: physicalDevice->getQueueFamilyProperties())
        strstream << "    queueFamilyProperties.timestampValidBits = " << queueFamilyProperties.timestampValidBits
                  << std::endl;

    const auto &limits = physicalDevice->getProperties().limits;
    strstream << "    limits.timestampComputeAndGraphics = " << limits.timestampComputeAndGraphics << std::endl;
    strstream << "    limits.timestampPeriod = " << limits.timestampPeriod << " nanoseconds." << std::endl;
    if (!physicalDevice)
        strstream << "Could not create PhysicalDevice" << std::endl;
    sendInfo(strstream.str());
#endif

    /* vsg::Names deviceExtensions; */
    /* vsg::QueueSettings queueSettings{vsg::QueueSetting{queueFamily, {1.0}}}; */

    /* auto deviceFeatures = vsg::DeviceFeatures::create(); */
    /* deviceFeatures->get().samplerAnisotropy = VK_TRUE; */
    /* deviceFeatures->get().geometryShader = enableGeometryShader; */

    /* auto device = vsg::Device::create(physicalDevice, queueSettings, validatedNames, deviceExtensions, deviceFeatures); */
    /* auto context = vsg::Context::create(device); */
    /* context->commandPool = vsg::CommandPool::create(device, queueFamily); */
    /* context->graphicsQueue = device->getQueue(queueFamily); */

    m_scenegraph = vsg::StateGroup::create();
    /* auto shaderSet = vsg::createPhongShaderSet(options); */
    /* auto graphicPipelineConfig = vsg::GraphicsPipelineConfigurator::create(shaderSet); */

    /* vsg::visit<SetPipelineStates>(graphicPipelineConfig); */

    /* // instantiate dynamicstate and add the state */
    /* auto dynamicState = vsg::DynamicState::create(); */
    /* dynamicState->dynamicStates.emplace_back(VK_DYNAMIC_STATE_LINE_WIDTH); */
    /* graphicPipelineConfig->pipelineStates.push_back(dynamicState); */

    /* // set up passing of material */
    /* auto mat = vsg::PhongMaterialValue::create(); */
    /* mat->value().diffuse.set(1.0f, 1.0f, 1.0f, 1.0f); */
    /* mat->value().specular.set(1.0f, 0.0f, 0.0f, 1.0f); // red specular highlight */

    /* graphicPipelineConfig->assignUniform("material", mat); */
    /* graphicPipelineConfig->copyTo(m_scenegraph); */

    auto stateCommands = setupVulkanGraphicsPipeline();
    for (auto &cmd: stateCommands)
        m_scenegraph->addChild(cmd);
    initScene(vsg::read_cast<vsg::Node>(vsg::Path("/home/hpcmdjur/git/vsgExamples/data/models/teapot.vsgt")), window);
    /* try { */
    /*     initScene(vsg::read_cast<vsg::Node>(vsg::Path("/home/hpcmdjur/git/vsgExamples/data/models/teapot.vsgt")), */
    /*               window); */
    /* } catch (vsg::Exception &e) { */
    /*     strstream << "VSG Exception: " << e.message << ". Result: " << e.result << std::endl; */
    /*     sendInfo(strstream.str()); */
    /* } */
}

VSGRenderer::~VSGRenderer()
{
#if DEBUG
    sendInfo("destroying");
#endif
}

std::vector<vsg::ref_ptr<vsg::StateCommand>> VSGRenderer::setupVulkanGraphicsPipeline()
{
    std::vector<vsg::ref_ptr<vsg::StateCommand>> stateCommands;

    // set up graphics pipeline
    // layout determines allocation size of descriptorsets
    vsg::DescriptorSetLayoutBindings descriptorBindings{
        {0, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, 1, VK_SHADER_STAGE_FRAGMENT_BIT,
         nullptr} // { binding, descriptorType, descriptorCount, stageFlags, pImmutableSamplers}
    };

    auto descriptorSetLayout = vsg::DescriptorSetLayout::create(descriptorBindings);

    vsg::PushConstantRanges pushConstantRanges{
        {VK_SHADER_STAGE_VERTEX_BIT, 0, 128}
        // projection, view, and model matrices, actual push constant calls automatically provided by the VSG's RecordTraversal
    };

    vsg::VertexInputState::Bindings vertexBindingsDescriptions{
        VkVertexInputBindingDescription{0, sizeof(vsg::vec3), VK_VERTEX_INPUT_RATE_VERTEX}, // vertex data
        VkVertexInputBindingDescription{1, sizeof(vsg::vec3), VK_VERTEX_INPUT_RATE_VERTEX}, // colour data
        VkVertexInputBindingDescription{2, sizeof(vsg::vec2), VK_VERTEX_INPUT_RATE_VERTEX} // tex coord data
    };

    vsg::VertexInputState::Attributes vertexAttributeDescriptions{
        VkVertexInputAttributeDescription{0, 0, VK_FORMAT_R32G32B32_SFLOAT, 0}, // vertex data
        VkVertexInputAttributeDescription{1, 1, VK_FORMAT_R32G32B32_SFLOAT, 0}, // colour data
        VkVertexInputAttributeDescription{2, 2, VK_FORMAT_R32G32_SFLOAT, 0} // tex coord data
    };

    // enable wireframe
    auto rasterization = vsg::RasterizationState::create();
    rasterization->polygonMode = VK_POLYGON_MODE_LINE;
    rasterization->lineWidth = 1.0f;
    rasterization->cullMode = VK_CULL_MODE_NONE;

    vsg::GraphicsPipelineStates pipelineStates{
        vsg::VertexInputState::create(vertexBindingsDescriptions, vertexAttributeDescriptions),
        vsg::InputAssemblyState::create(), rasterization,
        // vsg::RasterizationState::create(),
        vsg::MultisampleState::create(), vsg::ColorBlendState::create(), vsg::DepthStencilState::create()};

    // set up search paths to SPIRV shaders and textures
    vsg::Paths searchPaths = vsg::getEnvPaths("VISTLE_SHADER_PATH");

    // load shaders from vsgExamples => export VSG_FILE_PATH=<path to vsgExamples>/data
    vsg::ref_ptr<vsg::ShaderStage> vertexShader = vsg::ShaderStage::read(
        VK_SHADER_STAGE_VERTEX_BIT, "main", vsg::findFile("vert_PushConstants.spv", searchPaths));
    if (!vertexShader) {
        std::cerr << "Failed to load vertex shader" << std::endl;
        return {};
    }
    vsg::ref_ptr<vsg::ShaderStage> fragmentShader = vsg::ShaderStage::read(
        VK_SHADER_STAGE_FRAGMENT_BIT, "main", vsg::findFile("frag_PushConstants.spv", searchPaths));
    if (!fragmentShader) {
        std::cerr << "Failed to load fragment shader" << std::endl;
        return {};
    }
    auto pipelineLayout =
        vsg::PipelineLayout::create(vsg::DescriptorSetLayouts{descriptorSetLayout}, pushConstantRanges);
    auto graphicsPipeline =
        vsg::GraphicsPipeline::create(pipelineLayout, vsg::ShaderStages{vertexShader, fragmentShader}, pipelineStates);
    auto bindGraphicsPipeline = vsg::BindGraphicsPipeline::create(graphicsPipeline);

    // create texture image and associated DescriptorSets and binding
    // TODO: later
    /* auto texture = vsg::DescriptorImage::create(vsg::Sampler::create(), textureData, 0, 0, */
    /*                                             VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER); */

    /* auto descriptorSet = vsg::DescriptorSet::create(descriptorSetLayout, vsg::Descriptors{texture}); */
    /* auto bindDescriptorSet = */
    /*     vsg::BindDescriptorSet::create(VK_PIPELINE_BIND_POINT_GRAPHICS, pipelineLayout, 0, descriptorSet); */

    stateCommands.push_back(bindGraphicsPipeline);
    /* stateCommands.push_back(bindDescriptorSet); */

    return stateCommands;
}

void VSGRenderer::initScene(vsg::ref_ptr<vsg::Node> node, vsg::ref_ptr<vsg::Window> window)
{
    m_timesteps = VSGTimestepHandler::create(m_viewer, 1, 0);
    m_scenegraph->addChild(m_timesteps->root());

    if (node.valid())
        m_scenegraph->addChild(node);

    // set up the camera
    auto camera = createCameraForScene(m_scenegraph, 0, 0, window->extent2D().width, window->extent2D().height);
    auto main_view = vsg::View::create(camera, m_scenegraph);

    // add close handler to respond to pressing the window close window button and pressing escape =>
    // FIXME: at the moment vistle is blocking this event
    m_viewer->addEventHandler(vsg::CloseHandler::create(m_viewer));
    m_viewer->addEventHandler(vsg::WindowResizeHandler::create());
    // m_viewer->addEventHandler(AnimateTimestepSwitch::create(m_timesteps->animated(), m_viewer->start_point()));

    // add a trackball event handler to control the camera view use the mouse
    auto main_trackball = vsg::Trackball::create(camera);
    main_trackball->addWindow(window);
    m_viewer->addEventHandler(main_trackball);

    // create a command graph to render the scene on specified window
    auto commandGraph = vsg::createCommandGraphForView(window, camera, m_scenegraph);
    m_viewer->assignRecordAndSubmitTaskAndPresentation({commandGraph});
    m_viewer->setupThreading();

#if DEBUG
    sendInfo("Init Scenegraph: " + printVSGMetaData(m_scenegraph));
#endif

    // compile all the Vulkan objects and transfer data required to render the scene
    m_viewer->compile();
}

/* void VSGRenderer::flush() */
/* { */
/*     CERR << "flushing outstanding frames..." << std::endl; */
/*     for (size_t f = m_asyncFrames; f > 0; --f) { */
/*         composite(f - 1); */
/*     } */
/*     assert(m_numFramesToComposite == 0); */
/* } */

/* bool VSGRenderer::composite(size_t maxQueued) */
/* { */
/*     return true; */
/* } */

void VSGRenderer::prepareQuit()
{
#if DEBUG
    sendInfo("Prepare for quit");
#endif
    removeAllObjects();
    Renderer::prepareQuit();
    m_viewer->stopThreading();
    m_viewer->close();
}

bool VSGRenderer::render()
{
    if (!m_viewer->active())
        return false;

    if (!m_viewer->advanceToNextFrame())
        return false;

    m_viewer->handleEvents();
    m_viewer->recordAndSubmit();
    m_viewer->update();
    m_viewer->present();
    // fps
    /* auto fs = m_viewer->getFrameStamp(); */
    /* double fps = */
    /*     static_cast<double>(fs->frameCount) / */
    /*     std::chrono::duration<double, std::chrono::seconds::period>(vsg::clock::now() - m_viewer->start_point()).count(); */
    /* std::stringstream ss; */
    /* ss << "Average frame rate = " << fps << " fps" */
    /*           << std::endl; */
    /* sendInfo(ss.str()); */

    return true;
}

std::shared_ptr<vistle::RenderObject> VSGRenderer::addObject(int senderId, const std::string &senderPort,
                                                             vistle::Object::const_ptr container,
                                                             vistle::Object::const_ptr geometry,
                                                             vistle::Object::const_ptr normals,
                                                             vistle::Object::const_ptr texture)
{
    auto vro = std::make_shared<VsgRenderObject>(senderId, senderPort, container, geometry, normals, texture);
    auto t = vro->timestep;
    m_timesteps->addNode(vro->geo(), t);
    m_renderManager.addObject(vro);
    return vro;
}

void VSGRenderer::removeObject(std::shared_ptr<vistle::RenderObject> ro)
{
    auto vro = std::static_pointer_cast<VsgRenderObject>(ro);
    auto t = vro->timestep;
    m_timesteps->removeNode(vro->geo(), t);
    m_renderManager.removeObject(vro);
}

void VSGRenderer::connectionAdded(const vistle::Port *from, const vistle::Port *to)
{
#if DEBUG
    std::stringstream ss;
    ss << "new connection from " << *from << " to " << *to << std::endl;
    sendInfo(ss.str());
#endif
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
