#ifndef OPERATION_H
#define OPERATION_H

#include <iostream>
#include <vsg/all.h>

template<typename VSGGroupNodeType>
struct Merge: public vsg::Inherit<vsg::Operation, Merge<VSGGroupNodeType>> {
    Merge(vsg::observer_ptr<vsg::Viewer> in_viewer, vsg::ref_ptr<VSGGroupNodeType> in_attachmentPoint,
          vsg::ref_ptr<vsg::Node> in_node, const vsg::CompileResult &in_compileResult)
    : viewer(in_viewer), attachmentPoint(in_attachmentPoint), node(in_node), compileResult(in_compileResult)
    {}

    vsg::observer_ptr<vsg::Viewer> viewer;
    vsg::ref_ptr<VSGGroupNodeType> attachmentPoint;
    vsg::ref_ptr<vsg::Node> node;
    vsg::CompileResult compileResult;

    void run() override
    {
        std::cout << "Merge::run() attachmentpoint = " << attachmentPoint << ", " << node << std::endl;

        vsg::ref_ptr<vsg::Viewer> ref_viewer = viewer;
        if (ref_viewer)
            updateViewer(*ref_viewer, compileResult);

        attachmentPoint->addChild(node);
    }
};

/* struct LoadOperation: public vsg::Inherit<vsg::Operation, LoadOperation> { */
/*     LoadOperation(vsg::ref_ptr<vsg::Viewer> in_viewer, vsg::ref_ptr<vsg::Node> in_node, */
/*                   vsg::ref_ptr<vsg::Group> in_attachmentPoint, vsg::ref_ptr<vsg::Options> in_options) */
/*     : viewer(in_viewer), attachmentPoint(in_attachmentPoint), node(in_node), options(in_options) */
/*     {} */

/*     vsg::observer_ptr<vsg::Viewer> viewer; */
/*     vsg::ref_ptr<vsg::Group> attachmentPoint; */
/*     vsg::ref_ptr<vsg::Node> node; */
/*     vsg::ref_ptr<vsg::Options> options; */

/*     void run() override */
/*     { */
/*         vsg::ref_ptr<vsg::Viewer> ref_viewer = viewer; */

/*         // std::cout << "Loading " << filename << std::endl; */
/*         /1* if (auto node = vsg::read_cast<vsg::Node>(filename, options)) *1/ */
/*         /1* { *1/ */
/*         // std::cout << "Loaded " << filename << std::endl; */

/*         vsg::ComputeBounds computeBounds; */
/*         node->accept(computeBounds); */

/*         vsg::dvec3 centre = (computeBounds.bounds.min + computeBounds.bounds.max) * 0.5; */
/*         double radius = vsg::length(computeBounds.bounds.max - computeBounds.bounds.min) * 0.5; */
/*         auto scale = vsg::MatrixTransform::create(vsg::scale(1.0 / radius, 1.0 / radius, 1.0 / radius) * */
/*                                                   vsg::translate(-centre)); */

/*         scale->addChild(node); */

/*         auto result = ref_viewer->compileManager->compile(node); */
/*         if (result) */
/*             ref_viewer->addUpdateOperation(Merge::create(viewer, attachmentPoint, scale, result)); */
/*         /1* if (result) ref_viewer->addUpdateOperation(Merge::create(filename, viewer, attachmentPoint, scale, result)); *1/ */
/*         /1* } *1/ */
/*     } */
/* }; */

#endif
