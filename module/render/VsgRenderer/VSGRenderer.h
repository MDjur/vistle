//TODO:
//[] in the end look at vsgHeadless example
//[] load other filetypes with vsgxchange
//[] add serialization for vistle objects (example animal in vsg tutorial) to save vistle objects (scenes as vsg natives)

#ifndef VSGRENDERER_H
#define VSGRENDERER_H

#include "VSGTimestepHandler.h"

#include <vistle/renderer/renderer.h>
#include <vistle/renderer/parrendmgr.h>

#include <vsg/all.h>

class VSGRenderer: public vistle::Renderer {
public:
    VSGRenderer(const std::string &name, int moduleID, mpi::communicator comm);
    ~VSGRenderer() override;

private:
    /* bool compute() override; */
    /* bool compute(const std::shared_ptr<vistle::BlockTask> &task) const override; */
    /** 
     * @copydoc bool render()  
     */
    bool render() override;
    /** 
     * @copydoc void prepareQuit()  
     */
    void prepareQuit() override;
    /** 
     * @copydoc void connectionAdded(const vistle::Port *from, const vistle::Port *to)
     */
    void connectionAdded(const vistle::Port *from, const vistle::Port *to) override;
    /** 
     * @copydoc void connectionRemoved(const vistle::Port *from, const vistle::Port *to)
     */
    void connectionRemoved(const vistle::Port *from, const vistle::Port *to) override;
    /** 
     * @copydoc std::shared_ptr<vistle::RenderObject> addObject(int senderId, const std::string &senderPort,
                                                    vistle::Object::const_ptr container,
                                                    vistle::Object::const_ptr geometry,
                                                    vistle::Object::const_ptr normals,
                                                    vistle::Object::const_ptr texture)
     */
    std::shared_ptr<vistle::RenderObject> addObject(int senderId, const std::string &senderPort,
                                                    vistle::Object::const_ptr container,
                                                    vistle::Object::const_ptr geometry,
                                                    vistle::Object::const_ptr normals,
                                                    vistle::Object::const_ptr texture) override;
    /** 
     * @copydoc void connectionRemoved(const vistle::Port *from, const vistle::Port *to)
     */
    bool handleMessage(const vistle::message::Message *message, const vistle::MessagePayload &payload) override;
    /** 
     * @copydoc void removeObject(std::shared_ptr<vistle::RenderObject> ro)
     */
    void removeObject(std::shared_ptr<vistle::RenderObject> ro) override;
    void initScene(vsg::ref_ptr<vsg::Node> node, vsg::ref_ptr<vsg::Window> window);

    /**
     * @brief Create descriptorsets for textures, fragement and vertex shaders. 
     * Desriptorsets are used to define how resources (Image, Sampler, Buffer, Acceleration structure) can be accessed.
     * Each command in the graphics pipeline can access one or more suitable descriptorsets (commands which can handle descriptorSet) to access resources.
     * 
     * @return vector of statecommands to be added to the graphics pipeline
    */
    std::vector<vsg::ref_ptr<vsg::StateCommand>> setupVulkanGraphicsPipeline();
    /* bool changeParameter(const vistle::Parameter *p) override; */

    /* bool composite(size_t maxQueued); */
    /* void flush(); */

    vistle::ParallelRemoteRenderManager m_renderManager;
    vsg::ref_ptr<VSGTimestepHandler> m_timesteps;
    /* int m_asyncFrames; */
    vsg::ref_ptr<vsg::Viewer> m_viewer;
    vsg::ref_ptr<vsg::StateGroup> m_scenegraph;
    vsg::ref_ptr<vsg::MatrixTransform> m_transform;

    int m_numBlocks = 0;
    int m_numAnimationsSteps = 0;
    int m_numTotalTimesteps = 0;
    bool m_firstComputeCall = true;
};

#endif
