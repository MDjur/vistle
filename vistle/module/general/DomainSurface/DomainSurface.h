#ifndef DOMAINSURFACE_H
#define DOMAINSURFACE_H

#include <module/module.h>
#include <core/unstr.h>
#include <core/polygons.h>

class DomainSurface: public vistle::Module {

 public:
   DomainSurface(const std::string &name, int moduleID, mpi::communicator comm);
   ~DomainSurface();

private:
#if 0
   vistle::UnstructuredGrid::const_ptr m_grid_in;
   vistle::Polygons::ptr m_grid_out;
#endif
   typedef std::map<vistle::Index,vistle::Index> VerticesMapping;
   bool compute(std::shared_ptr<vistle::PortTask> task) const override;
   vistle::Polygons::ptr createSurface(vistle::UnstructuredGrid::const_ptr m_grid_in, VerticesMapping &vm) const;
   //bool checkNormal(vistle::Index v1, vistle::Index v2, vistle::Index v3, vistle::Scalar x_center, vistle::Scalar y_center, vistle::Scalar z_center);
};

#endif
