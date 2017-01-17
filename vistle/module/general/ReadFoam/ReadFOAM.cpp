/**************************************************************************\
 **                                                           (C)2013 RUS  **
 **                                                                        **
 ** Description: Read FOAM data format                                     **
 **                                                                        **
 **                                                                        **
 **                                                                        **
 **                                                                        **
 **                                                                        **
 **                                                                        **
 ** History:                                                               **
 ** May   13	    C.Kopf  	    V1.0                                   **
 *\**************************************************************************/

#include "ReadFOAM.h"
#include <core/unstr.h>
#include <core/vec.h>
#include <core/message.h>

//Includes copied from covise ReadFOAM.cpp
#include <sstream>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <set>
#include <cctype>

#include <ctime>
#include <cstdio>
#include <cstdlib>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include <memory>

#include "foamtoolbox.h"
#include <util/coRestraint.h>
#include <boost/serialization/vector.hpp>
#include <boost/mpi.hpp>
#include <unordered_set>
namespace mpi = boost::mpi;

using namespace vistle;

ReadFOAM::ReadFOAM(const std::string &shmname, const std::string &name, int moduleId)
: Module("ReadFoam", shmname, name, moduleId)
, m_boundOut(nullptr)
{
   // file browser parameter
   m_casedir = addStringParameter("casedir", "OpenFOAM case directory",
      "/data/OpenFOAM", Parameter::Directory);
   //Time Parameters
   m_starttime = addFloatParameter("starttime", "start reading at the first step after this time", 0.);
   setParameterMinimum<Float>(m_starttime, 0.);
   m_stoptime = addFloatParameter("stoptime", "stop reading at the last step before this time",
         std::numeric_limits<double>::max());
   setParameterMinimum<Float>(m_stoptime, 0.);
   m_timeskip = addIntParameter("timeskip", "skip this many timesteps after reading one", 0);
   setParameterMinimum<Integer>(m_timeskip, 0);
   m_readGrid = addIntParameter("read_grid", "load the grid?", 1, Parameter::Boolean);

   //Mesh ports
   m_boundOut = createOutputPort("grid_out1");

   for (int i=0; i<NumPorts; ++i) {
      {// Data Ports
         std::stringstream s;
         s << "data_out" << i;
         m_volumeDataOut.push_back(createOutputPort(s.str()));
      }
      {// Date Choice Parameters
         std::stringstream s;
         s << "Data" << i;
         auto p =  addStringParameter(s.str(), "name of field", "(NONE)", Parameter::Choice);
         std::vector<std::string> choices;
         choices.push_back("(NONE)");
         setParameterChoices(p, choices);
         m_fieldOut.push_back(p);
      }
   }
   m_readBoundary = addIntParameter("read_boundary", "load the boundary?", 1, Parameter::Boolean);
   m_boundaryPatchesAsVariants = addIntParameter("patches_as_variants", "create sub-objects with variant attribute for boundary patches", 1, Parameter::Boolean);
   m_patchSelection = addStringParameter("patches", "select patches","all");
   for (int i=0; i<NumBoundaryPorts; ++i) {
      {// 2d Data Ports
         std::stringstream s;
         s << "data_2d_out" << i;
         m_boundaryDataOut.push_back(createOutputPort(s.str()));
      }
      {// 2d Data Choice Parameters
         std::stringstream s;
         s << "Data2d" << i;
         auto p =  addStringParameter(s.str(), "name of field", "(NONE)", Parameter::Choice);
         std::vector<std::string> choices;
         choices.push_back("(NONE)");
         setParameterChoices(p, choices);
         m_boundaryOut.push_back(p);
      }
   }
   m_buildGhostcellsParam = addIntParameter("build_ghostcells", "whether to build ghost cells", 1, Parameter::Boolean);
}


ReadFOAM::~ReadFOAM()       //Destructor
{
}

std::vector<std::string> ReadFOAM::getFieldList() const {

   std::vector<std::string> choices;
   choices.push_back("(NONE)");

   if (m_case.valid) {
      for (auto &field: m_case.constantFields)
         choices.push_back(field.first);
      for (auto &field: m_case.varyingFields)
         choices.push_back(field.first);
   }

   return choices;
}

int ReadFOAM::rankForBlock(int processor) const {

   if (m_case.numblocks == 0)
      return 0;

   if (processor == -1)
      return -1;

   return processor % size();
}

bool ReadFOAM::changeParameter(const Parameter *p)
{
   auto sp = dynamic_cast<const StringParameter *>(p);
   if (sp == m_casedir) {
      std::string casedir = sp->getValue();

      m_case = getCaseInfo(casedir);
      if (!m_case.valid) {
         std::cerr << casedir << " is not a valid OpenFOAM case" << std::endl;
         return false;
      }

      std::cerr << "# processors: " << m_case.numblocks << std::endl;
      std::cerr << "# time steps: " << m_case.timedirs.size() << std::endl;
      std::cerr << "grid topology: " << (m_case.varyingGrid?"varying":"constant") << std::endl;
      std::cerr << "grid coordinates: " << (m_case.varyingCoords?"varying":"constant") << std::endl;

      //print out a list of boundary patches to Vistle Console
      if (rank() == 0) {
         std::stringstream meshdir;
         meshdir << casedir << "/constant/polyMesh"; //<< m_case.constantdir << "/polyMesh";
         Boundaries bounds = loadBoundary(meshdir.str());
         if (bounds.valid) {
            sendInfo("boundary patches:");
            for (Index i=0;i<bounds.boundaries.size();++i) {
               std::stringstream info;
               info << bounds.boundaries[i].index<< " ## " << bounds.boundaries[i].name;
               sendInfo("%s", info.str().c_str());
            }
         } else {
            sendInfo("No global boundary file was found at:");
            sendInfo(meshdir.str());
         }
      }

      //fill choice parameters
      std::vector<std::string> choices = getFieldList();
      for (auto out: m_fieldOut) {
         setParameterChoices(out, choices);
      }
      for (auto out: m_boundaryOut) {
         setParameterChoices(out, choices);
      }
   }

   return Module::changeParameter(p);
}

bool loadCoords(const std::string &meshdir, Coords::ptr grid) {

   std::shared_ptr<std::istream> pointsIn = getStreamForFile(meshdir, "points");
   if (!pointsIn)
      return false;
   HeaderInfo pointsH = readFoamHeader(*pointsIn);
   grid->setSize(pointsH.lines);
   if (!readFloatVectorArray(pointsH, *pointsIn, grid->x().data(), grid->y().data(), grid->z().data(), pointsH.lines)) {
      std::cerr << "readFloatVectorArray for " << meshdir << "/points failed" << std::endl;
      return false;
   }
   return true;
}

GridDataContainer ReadFOAM::loadGrid(const std::string &meshdir, std::string topologyDir) {

   //std::cerr << "loadGrid(\"" << meshdir << "\", \"" << topologyDir << "\")" << std::endl;
   if (topologyDir.empty())
       topologyDir = meshdir;
   bool readGrid = m_readGrid->getValue();
   bool readBoundary = m_readBoundary->getValue();
   bool patchesAsVariants = m_boundaryPatchesAsVariants->getValue();

   std::shared_ptr<Boundaries> boundaries(new Boundaries());
   *boundaries = loadBoundary(topologyDir);
   UnstructuredGrid::ptr grid(new UnstructuredGrid(0, 0, 0));
   std::vector<Polygons::ptr> polyList;
   size_t numPatches = 0;
   for (const auto &b: boundaries->boundaries) {
       int boundaryIndex=b.index;
       if (b.numFaces>0 && m_boundaryPatches(boundaryIndex)) {
           ++numPatches;
       }
   }
   for (size_t i=0; i<(patchesAsVariants ? numPatches : 1); ++i) {
       polyList.emplace_back(new Polygons(0, 0, 0));
   }
   std::shared_ptr<std::vector<Index> > owners(new std::vector<Index>());
   GridDataContainer result(grid,polyList,owners,boundaries);
   if (!readGrid && !readBoundary) {
      return result;
   }

   //read mesh files
   std::shared_ptr<std::istream> ownersIn = getStreamForFile(topologyDir, "owner");
   if (!ownersIn)
      return result;
   HeaderInfo ownerH = readFoamHeader(*ownersIn);
   DimensionInfo dim = parseDimensions(ownerH.header);
   owners->resize(ownerH.lines);
   if (!readIndexArray(ownerH, *ownersIn, (*owners).data(), (*owners).size())) {
      std::cerr << "readIndexArray for " << topologyDir << "/owner failed" << std::endl;
      return result;
   }

   {

      std::shared_ptr<std::istream> facesIn = getStreamForFile(topologyDir, "faces");
      if (!facesIn)
         return result;
      HeaderInfo facesH = readFoamHeader(*facesIn);
      std::vector<std::vector<Index>> faces(facesH.lines);
      if (!readIndexListArray(facesH, *facesIn, faces.data(), faces.size())) {
         std::cerr << "readIndexListArray for " << topologyDir << "/faces failed" << std::endl;
         return result;
      }

      std::shared_ptr<std::istream> neighboursIn = getStreamForFile(topologyDir, "neighbour");
      if (!neighboursIn)
         return result;
      HeaderInfo neighbourH = readFoamHeader(*neighboursIn);
      if (neighbourH.lines != dim.internalFaces) {
         std::cerr << "inconsistency: #internalFaces != #neighbours (" << dim.internalFaces << " != " << neighbourH.lines << ")" << std::endl;
         return result;
      }
      std::vector<Index> neighbours(neighbourH.lines);
      if (!readIndexArray(neighbourH, *neighboursIn, neighbours.data(), neighbours.size()))
         return result;

      //Boundary Polygon
      if (readBoundary) {
          if (patchesAsVariants) {
              size_t boundIdx=0;
              for (const auto &b: boundaries->boundaries) {
                  int boundaryIndex=b.index;
                  if (b.numFaces>0 && m_boundaryPatches(boundaryIndex)) {
                      const auto &poly = polyList[boundIdx];
                      auto &polys = poly->el();
                      auto &conn = poly->cl();
                      polys.reserve(b.numFaces+1);

                      for (index_t i=b.startFace; i<b.startFace + b.numFaces; ++i) {
                          auto &face = faces[i];
                          for (Index j=0; j<face.size(); ++j) {
                              conn.push_back(face[j]);
                          }
                          polys.push_back(conn.size());
                      }
                      ++boundIdx;
                  }
              }
          } else {
              Index num_bound = 0;
              for (const auto &b: boundaries->boundaries) {
                  int boundaryIndex=b.index;
                  if (m_boundaryPatches(boundaryIndex)) {
                      num_bound+=b.numFaces;
                  }
              }
              const auto &poly = polyList[0];
              auto &polys = poly->el();
              auto &conn = poly->cl();
              polys.reserve(num_bound+1);
              for (const auto &b: boundaries->boundaries) {
                  int boundaryIndex=b.index;
                  if (m_boundaryPatches(boundaryIndex) && b.numFaces>0) {
                      for (index_t i=b.startFace; i<b.startFace + b.numFaces; ++i) {
                          auto &face = faces[i];
                          for (Index j=0; j<face.size(); ++j) {
                              conn.push_back(face[j]);
                          }
                          polys.push_back(conn.size());
                      }
                  }
              }
          }
      }

      //Grid
      if (readGrid) {
         grid->el().resize(dim.cells+1);
         grid->tl().resize(dim.cells);
         //Create CellFaceMap
         std::vector<std::vector<Index>> cellfacemap(dim.cells);
         for (Index face = 0; face < (*owners).size(); ++face) {
            cellfacemap[(*owners)[face]].push_back(face);
         }
         for (Index face = 0; face < neighbours.size(); ++face) {
            cellfacemap[neighbours[face]].push_back(face);
         }

         //Vertices lists for GhostCell creation -
         //each node creates lists of the outer vertices that are shared with other domains
         //with either its own or the neighbouring domain's face-numbering (clockwise or ccw)
         //so that two domains have the same list for a mutual border
         //therefore m_procBoundaryVertices[0][1] point to the same vertices in the same order as
         //m_procBoundaryVertices[1][0] (though each list may use different labels  for each vertice)
         if (m_buildGhost) {
            for (const auto &b: boundaries->procboundaries) {
               std::vector<Index> outerVertices;
               int myProc=b.myProc;
               int neighborProc=b.neighborProc;
               if (myProc < neighborProc) {
                  //create with own numbering
                  for (Index i=b.startFace; i<b.startFace+b.numFaces; ++i) {
                     auto face=faces[i];
                     for (Index j=0; j<face.size(); ++j) {
                        outerVertices.push_back(face[j]);
                     }
                  }
               } else {
                  //create with neighbour numbering (reverse direction)
                  for (Index i=b.startFace; i<b.startFace+b.numFaces; ++i) {
                     auto face=faces[i];
                     outerVertices.push_back(face[0]);
                     for (Index j=face.size()-1; j>0; --j) {
                        outerVertices.push_back(face[j]);
                     }
                  }
               }

               //check for ghost cells recursively
               std::unordered_set <Index> ghostCellCandidates;
               std::unordered_set <Index> notGhostCells;
               for (Index i=0;i<b.numFaces;++i) {
                  Index cell=(*owners)[b.startFace + i];
                  ghostCellCandidates.insert(cell);
               }
//               std::sort(ghostCellCandidates.begin(),ghostCellCandidates.end()); //Sort Vector by ascending Value
//               ghostCellCandidates.erase(std::unique(ghostCellCandidates.begin(), ghostCellCandidates.end()), ghostCellCandidates.end()); //Delete duplicate entries
               for (Index i=0;i<b.numFaces;++i) {
                  Index cell=(*owners)[b.startFace + i];
                  std::vector<Index> adjacentCells=getAdjacentCells(cell,dim,cellfacemap,*owners,neighbours);
                  for (Index j=0; j<adjacentCells.size(); ++j) {
                     if (!checkCell(adjacentCells[j],ghostCellCandidates,notGhostCells,dim,outerVertices,cellfacemap,faces,*owners,neighbours))
                        std::cerr << "ERROR finding GhostCellCandidates" << std::endl;
                  }
               }
               m_procGhostCellCandidates[myProc][neighborProc] = ghostCellCandidates;
               m_procBoundaryVertices[myProc][neighborProc] = outerVertices;
            }
         }

         auto types = grid->tl().data();
         Index num_conn = 0;
         //Check Shape of Cells and fill Type_List
         for (index_t i=0; i<dim.cells; i++) {
            const std::vector<Index> &cellfaces=cellfacemap[i];

            bool onlySimpleFaces = true; // only faces with 3 or 4 corners
            std::vector<Index> threeVert, fourVert;
            for (index_t j=0; j<cellfaces.size(); ++j) {
               if (faces[cellfaces[j]].size() == 4) {
                   fourVert.push_back(j);
                   if (fourVert.size() > 6)
                      break;
               } else if (faces[cellfaces[j]].size() == 3) {
                   threeVert.push_back(j);
                   if (threeVert.size() > 4)
                      break;
               } else {
                   onlySimpleFaces = false;
                   break;
               }
            }
            const Index num_faces = cellfaces.size();
            Index num_verts = 0;
            if (num_faces==6 && fourVert.size()==6 && threeVert.size() == 0 && onlySimpleFaces) {
               types[i]=UnstructuredGrid::HEXAHEDRON;
               num_verts = 8;
            } else if (num_faces==5 && fourVert.size()==3 && threeVert.size()==2 && onlySimpleFaces) {
               types[i]=UnstructuredGrid::PRISM;
               num_verts = 6;
            } else if (num_faces==5 && fourVert.size()==1 && threeVert.size()==4 && onlySimpleFaces) {
               types[i]=UnstructuredGrid::PYRAMID;
               num_verts = 5;
            } else if (num_faces==4 && fourVert.size()==0 && threeVert.size()==4 && onlySimpleFaces) {
               types[i]=UnstructuredGrid::TETRAHEDRON;
               num_verts = 4;
            } else {
               types[i]=UnstructuredGrid::POLYHEDRON;
               for (Index j=0; j<cellfaces.size(); ++j) {
                  num_verts += faces[cellfaces[j]].size() + 1;
               }
            }
            num_conn += num_verts;
         }
         //save data cell by cell to element, connectivity and type list
         auto el = grid->el().data();
         auto &connectivities = grid->cl();
         auto inserter = std::back_inserter(connectivities);
         connectivities.reserve(num_conn);
         for(index_t i=0;  i<dim.cells; i++) {
            //element list
            el[i] = connectivities.size();
            //connectivity list
            const auto &cellfaces=cellfacemap[i];//get all faces of current cell
            switch (types[i]) {
               case UnstructuredGrid::HEXAHEDRON: {
                  index_t ia=cellfaces[0];//choose the first face as starting face
                  std::vector<index_t> a=faces[ia];

                  //vistle requires that the vertices-numbering of the first face conforms with the right hand rule (pointing into the cell)
                  //so the starting face is tested if it does and the numbering is reversed if it doesn't
                  if (!isPointingInwards(ia,i,dim.internalFaces,(*owners),neighbours)) {
                     std::reverse(a.begin(), a.end());
                  }

                  // bottom face
                  std::copy(a.begin(), a.end(), inserter);
#if 0
                  int idx2Common = -1, idxDisjoint = -1;
                  int commonIndMy[2] = { -1, -1}, commonIndOther[2] = { -1, -1};
                  Index commonVerts[2];
                  for (int f=1; f<6; ++f) {
                      const auto &face = faces[cellfaces[f]];
                      int numCommon = 0;
                      for (int i=0; i<4; ++i) {
                          for (int j=0; j<4; ++j) {
                              if (a[j] == face[i] && numCommon < 2) {
                                  commonVerts[numCommon] = a[j];
                                  commonIndOther[numCommon] = i;
                                  commonIndMy[numCommon] = j;
                                  ++numCommon;
                              }
                          }
                      }
                      vassert(numCommon == 0 || numCommon == 2);
                      if (numCommon == 0) {
                          idxDisjoint = f;
                          if (idx2Common >= 0)
                              break;
                      } else if (numCommon == 2) {
                          idx2Common = f;
                          if (idxDisjoint >= 0)
                              break;
                      }
                  }
                  const auto &adjoining_face = faces[cellfaces[idx2Common]];
                  // top face - bring into correct order
                  const auto &opposite_face = faces[cellfaces[idxDisjoint]];
                  int otherVerts[2] = { -1, -1 };
                  bool reverse = false;
                  int startIdx = -1;
                  if (commonIndOther[0]+1 == commonIndOther[1]) {
                  } else if (commonIndOther[1]+3 == commonIndOther[0]) {
                  } else if (commonIndOther[0]+3 == commonIndOther[1]) {
                  } else {
                  }
#else
                  connectivities.push_back(findVertexAlongEdge(a[0],ia,cellfaces,faces));
                  connectivities.push_back(findVertexAlongEdge(a[1],ia,cellfaces,faces));
                  connectivities.push_back(findVertexAlongEdge(a[2],ia,cellfaces,faces));
                  connectivities.push_back(findVertexAlongEdge(a[3],ia,cellfaces,faces));
#endif
               }
               break;

               case UnstructuredGrid::PRISM: {
                  index_t it=1;
                  index_t ia=cellfaces[0];
                  while (faces[ia].size()>3) {//find first face with 3 vertices to use as starting face
                     ia=cellfaces[it++];
                  }

                  std::vector<index_t> a=faces[ia];

                  if(!isPointingInwards(ia,i,dim.internalFaces,(*owners),neighbours)) {
                     std::reverse(a.begin(), a.end());
                  }

                  std::copy(a.begin(), a.end(), inserter);
                  connectivities.push_back(findVertexAlongEdge(a[0],ia,cellfaces,faces));
                  connectivities.push_back(findVertexAlongEdge(a[1],ia,cellfaces,faces));
                  connectivities.push_back(findVertexAlongEdge(a[2],ia,cellfaces,faces));
               }
               break;

               case UnstructuredGrid::PYRAMID: {
                  index_t it=1;
                  index_t ia=cellfaces[0];
                  while (faces[ia].size()<4) {//find the rectangular face to use as starting face
                     ia=cellfaces[it++];
                  }

                  std::vector<index_t> a=faces[ia];

                  if(!isPointingInwards(ia,i,dim.internalFaces,(*owners),neighbours)) {
                     std::reverse(a.begin(), a.end());
                  }

                  std::copy(a.begin(), a.end(), inserter);
                  connectivities.push_back(findVertexAlongEdge(a[0],ia,cellfaces,faces));
               }
               break;

               case UnstructuredGrid::TETRAHEDRON: {
                  index_t ia=cellfaces[0];//use first face as starting face
                  std::vector<index_t> a=faces[ia];

                  if(!isPointingInwards(ia,i,dim.internalFaces,(*owners),neighbours)) {
                     std::reverse(a.begin(), a.end());
                  }

                  std::copy(a.begin(), a.end(), inserter);
                  connectivities.push_back(findVertexAlongEdge(a[0],ia,cellfaces,faces));
               }
               break;

               case UnstructuredGrid::POLYHEDRON: {
                  for (index_t j=0;j<cellfaces.size();j++) {
                     index_t ia=cellfaces[j];
                     std::vector<index_t> a=faces[ia];

                     if(!isPointingInwards(ia,i,dim.internalFaces,(*owners),neighbours)) {
                        std::reverse(a.begin(), a.end());
                     }

                     connectivities.push_back(a.size());
                     for (index_t k=0; k<a.size(); ++k) {
                        connectivities.push_back(a[k]);
                     }
                  }
               }
               break;
            }
         }
         el[dim.cells] = connectivities.size();
      }
   }

   if (readGrid) {
      loadCoords(meshdir, grid);
      grid->checkConvexity();

      if (readBoundary) {
          //if grid has been read already and boundary polygons are read also -> re-use coordinate lists for the boundary-plygon
          for (auto &poly: polyList) {
              poly->d()->x[0] = grid->d()->x[0];
              poly->d()->x[1] = grid->d()->x[1];
              poly->d()->x[2] = grid->d()->x[2];
          }
      }
   } else {
       //else read coordinate lists just for boundary polygons
       bool first = true;
       for (auto &poly: polyList) {
           if (first) {
               loadCoords(meshdir, poly);
           } else {
              poly->d()->x[0] = polyList[0]->d()->x[0];
              poly->d()->x[1] = polyList[0]->d()->x[1];
              poly->d()->x[2] = polyList[0]->d()->x[2];
           }
       }
   }
   return result;
}

DataBase::ptr ReadFOAM::loadField(const std::string &meshdir, const std::string &field) {

   std::shared_ptr<std::istream> stream = getStreamForFile(meshdir, field);
   if (!stream) {
      std::cerr << "failed to open " << meshdir << "/" << field << std::endl;
      return DataBase::ptr();
   }
   HeaderInfo header = readFoamHeader(*stream);
   if (header.fieldclass == "volScalarField") {
      Vec<Scalar>::ptr s(new Vec<Scalar>(header.lines));
      if (!readFloatArray(header, *stream, s->x().data(), s->x().size())) {
         std::cerr << "readFloatArray for " << meshdir << "/" << field << " failed" << std::endl;
         return DataBase::ptr();
      }
      return s;
   } else if (header.fieldclass == "volVectorField") {
      Vec<Scalar, 3>::ptr v(new Vec<Scalar, 3>(header.lines));
      if (!readFloatVectorArray(header, *stream, v->x().data(), v->y().data(), v->z().data(), v->x().size())) {
         std::cerr << "readFloatVectorArray for " << meshdir << "/" << field << " failed" << std::endl;
         return DataBase::ptr();
      }
      return v;
   }

   std::cerr << "cannot interpret " << meshdir << "/" << field << std::endl;
   return DataBase::ptr();
}

std::vector<DataBase::ptr> ReadFOAM::loadBoundaryField(const std::string &meshdir, const std::string &field,
                                        const int &processor) {

   bool asVariants = m_boundaryPatchesAsVariants->getValue();
   auto &boundaries = *m_boundaries[processor];
   auto owners = m_owners[processor]->data();

   //Create the dataMapping Vector
   std::vector<std::string> patchNames;
   std::vector<std::vector<index_t>> dataMapping;
   if (asVariants) {
       size_t idx=0;
       for (const auto &b: boundaries.boundaries) {
           Index boundaryIndex=b.index;
           if (b.numFaces>0 && m_boundaryPatches(boundaryIndex)) {
               dataMapping.resize(idx+1);
               for (index_t i=b.startFace; i<b.startFace + b.numFaces; ++i) {
                   dataMapping[idx].push_back(owners[i]);
               }
               patchNames.emplace_back(b.name);
               ++idx;
           }
       }
   } else {
       dataMapping.resize(1);
       for (const auto &b: boundaries.boundaries) {
           Index boundaryIndex=b.index;
           if (b.numFaces>0 && m_boundaryPatches(boundaryIndex)) {
               for (index_t i=b.startFace; i<b.startFace + b.numFaces; ++i) {
                   dataMapping[0].push_back(owners[i]);
               }
           }
       }
   }

   std::shared_ptr<std::istream> stream = getStreamForFile(meshdir, field);
   if (!stream) {
      std::cerr << "failed to open " << meshdir << "/" << field << std::endl;
   }
   HeaderInfo header = readFoamHeader(*stream);
   std::vector<scalar_t> fullX(header.lines),fullY(header.lines),fullZ(header.lines);
   if (header.fieldclass == "volScalarField") {
      fullX.resize(header.lines);
      if (!readFloatArray(header, *stream, fullX.data(), header.lines)) {
         std::cerr << "readFloatArray for " << meshdir << "/" << field << " failed" << std::endl;
      }
   } else if (header.fieldclass == "volVectorField") {
      fullX.resize(header.lines);
      fullY.resize(header.lines);
      fullZ.resize(header.lines);
      if (!readFloatVectorArray(header, *stream,fullX.data(),fullY.data(),fullZ.data(),header.lines)) {
         std::cerr << "readFloatVectorArray for " << meshdir << "/" << field << " failed" << std::endl;
      }
   } else {
       std::cerr << "cannot interpret " << meshdir << "/" << field << std::endl;
   }

   std::vector<DataBase::ptr> result;

   for (size_t idx=0; idx<dataMapping.size(); ++idx) {
       if (header.fieldclass == "volScalarField") {
           Vec<Scalar>::ptr s(new Vec<Scalar>(dataMapping[idx].size()));
           auto x = s->x().data();
           for (index_t i=0;i<dataMapping[idx].size();++i) {
               x[i] = fullX[dataMapping[idx][i]];
           }
           if (asVariants)
               s->addAttribute("_variant", patchNames[idx]);
           result.push_back(s);

       } else if (header.fieldclass == "volVectorField") {
           Vec<Scalar, 3>::ptr v(new Vec<Scalar, 3>(dataMapping[idx].size()));
           auto x = v->x().data();
           auto y = v->y().data();
           auto z = v->z().data();
           for (index_t i=0;i<dataMapping[idx].size();++i) {
               x[i] = fullX[dataMapping[idx][i]];
               y[i] = fullY[dataMapping[idx][i]];
               z[i] = fullZ[dataMapping[idx][i]];
           }
           if (asVariants)
               v->addAttribute("_variant", patchNames[idx]);
           result.push_back(v);
       }
   }
   return result;
}

void ReadFOAM::setMeta(Object::ptr obj, int processor, int timestep) const {

   if (obj) {
      Index skipfactor = m_timeskip->getValue()+1;
      obj->setTimestep(timestep);
      obj->setNumTimesteps(m_case.timedirs.size()/skipfactor);
      obj->setBlock(processor);
      obj->setNumBlocks(m_case.numblocks == 0 ? 1 : m_case.numblocks);

      if (timestep >= 0) {
         int i = 0;
         for (auto &ts: m_case.timedirs) {
            if (i == timestep) {
               obj->setRealTime(ts.first);
               break;
            }
            ++i;
         }
      }
   }
}

bool ReadFOAM::loadFields(const std::string &meshdir, const std::map<std::string, int> &fields, int processor, int timestep) {
   for (int i=0; i<NumPorts; ++i) {
      std::string field = m_fieldOut[i]->getValue();
      auto it = fields.find(field);
      if (it == fields.end())
         continue;
      DataBase::ptr obj = loadField(meshdir, field);
      setMeta(obj, processor, timestep);
      m_currentvolumedata[processor][i]= obj;
   }

   for (int i=0; i<NumBoundaryPorts; ++i) {
      std::string field = m_boundaryOut[i]->getValue();
      auto it = fields.find(field);
      if (it == fields.end())
         continue;
      auto fields = loadBoundaryField(meshdir, field, processor);
      if (fields.size() != m_currentbound[processor].size()) {
          std::cerr << "MISMATCH: trying to load field " << field << " in " << meshdir << " on proc " << processor << " for t=" << timestep << std::endl;
          std::cerr << "MISMATCH: fields.size()=" << fields.size() << ", curbound[proc].size()=" << m_currentbound[processor].size() << std::endl;
      }
      vassert(fields.size() == m_currentbound[processor].size());
      for (size_t j=0; j<fields.size(); ++j) {
          auto &obj = fields[j];
          setMeta(obj, processor, timestep);
          obj->setGrid(m_currentbound[processor][j]);
          obj->setMapping(DataBase::Element);
          addObject(m_boundaryDataOut[i], obj);
      }
   }
   return true;
}




bool ReadFOAM::readDirectory(const std::string &casedir, int processor, int timestep) {
   std::string dir = casedir;

   if (processor >= 0) {
      std::stringstream s;
      s << "/processor" << processor;
      dir += s.str();
   }

   if (timestep < 0) {
      dir += "/" + m_case.constantdir;
      if (!m_case.varyingGrid){
         auto ret = loadGrid(dir + "/polyMesh");
         UnstructuredGrid::ptr grid = ret.grid;
         setMeta(grid, processor, timestep);
         for (auto &poly: ret.polygon)
             setMeta(poly, processor, timestep);
         m_owners[processor] = ret.owners;
         m_boundaries[processor] = ret.boundaries;

         m_currentgrid[processor] = grid;
         m_currentbound[processor] = ret.polygon;
         m_basedir[processor] = dir;
      }
      loadFields(dir, m_case.constantFields, processor, timestep);
   } else {
      Index i = 0;
      Index skipfactor = m_timeskip->getValue()+1;
      std::string completeMeshDir;
      for (auto &ts: m_case.timedirs) {
         if (i == timestep*skipfactor) {
            completeMeshDir = dir;
            auto it = m_case.completeMeshDirs.find(ts.first);
            if (it != m_case.completeMeshDirs.end()) {
                completeMeshDir = dir + "/" + it->second;
            }
            dir += "/" + ts.second;
            break;
         }
         ++i;
      }
      if (i == m_case.timedirs.size()) {
         std::cerr << "no directory for timestep " << timestep << " found" << std::endl;
         return false;
      }
      if (m_case.varyingGrid || m_case.varyingCoords) {
         UnstructuredGrid::ptr grid;
         std::vector<Polygons::ptr> polygons;
         if (completeMeshDir == m_basedir[processor]) {
            {
               grid.reset(new UnstructuredGrid(0, 0, 0));
               UnstructuredGrid::Data *od = m_currentgrid[processor]->d();
               UnstructuredGrid::Data *nd = grid->d();
               nd->tl = od->tl;
               nd->el = od->el;
               nd->cl = od->cl;
            }
            loadCoords(dir + "/polyMesh", grid);
            {
                for (size_t j=0; j<m_currentbound[processor].size(); ++j) {
                    Polygons::ptr poly(new Polygons(0, 0, 0));
                    Polygons::Data *od = m_currentbound[processor][j]->d();
                    Polygons::Data *nd = poly->d();
                    nd->el = od->el;
                    nd->cl = od->cl;
                    for (int i=0; i<3; ++i)
                        poly->d()->x[i] = grid->d()->x[i];
                    polygons.push_back(poly);
                }
            }
         } else {
            auto ret = loadGrid(dir + "/polyMesh", completeMeshDir + "/polyMesh");
            grid = ret.grid;
            polygons = ret.polygon;
            m_owners[processor] = ret.owners;
            m_boundaries[processor] = ret.boundaries;
            m_basedir[processor] = completeMeshDir;
         }
         setMeta(grid, processor, timestep);
         for (auto &poly: polygons)
             setMeta(poly, processor, timestep);
         m_currentgrid[processor] = grid;
         m_currentbound[processor] = polygons;
      }
      loadFields(dir, m_case.varyingFields, processor, timestep);
   }

   return true;
}

int tag(int p, int n, int i=0) { //MPI needs a unique ID for each pair of send/receive request, this function creates unique ids for each processor pairing
   return p*10000+n*100+i;
}

std::vector<Index> ReadFOAM::getAdjacentCells(const Index &cell,
                                              const DimensionInfo &dim,
                                              const std::vector<std::vector<Index>> &cellfacemap,
                                              const std::vector<Index> &owners,
                                              const std::vector<Index> &neighbours) {
   const std::vector<Index> &cellfaces=cellfacemap[cell];
   std::vector<Index> adjacentCells;
   for (Index i=0; i<cellfaces.size(); ++i) {
      if (cellfaces[i] < dim.internalFaces) {
         Index adjCell;
         Index o=owners[cellfaces[i]];
         Index n=neighbours[cellfaces[i]];
         adjCell= (o==cell ? n : o);
         adjacentCells.push_back(adjCell);
      }
   }
   return adjacentCells;
}

bool ReadFOAM::checkCell(const Index cell,
                         std::unordered_set<Index> &ghostCellCandidates,
                         std::unordered_set<Index> &notGhostCells,
                         const DimensionInfo &dim,
                         const std::vector<Index> &outerVertices,
                         const std::vector<std::vector<Index>> &cellfacemap,
                         const std::vector<std::vector<Index>> &faces,
                         const std::vector<Index> &owners,
                         const std::vector<Index> &neighbours) {

   if (notGhostCells.count(cell) == 1) {// if cell is already known to not be a ghost-cell
      return true;
   }

   if (ghostCellCandidates.count(cell) == 1) {// if cell is not an already known ghost-cell
      return true;
   }

   bool isGhostCell=false;
   const std::vector<Index> &cellfaces=cellfacemap[cell];
   const auto cellvertices = getVerticesForCell(cellfaces, faces);
   for (auto it = cellvertices.begin(); it!=cellvertices.end(); ++it) {//check if the cell has a vertex on the outer boundary
      if (std::find(outerVertices.begin(), outerVertices.end(), *it) != outerVertices.end()) {
         isGhostCell=true;
         break;
      }
   }

   if (isGhostCell) {
      ghostCellCandidates.insert(cell);
      std::vector<Index> adjacentCells=getAdjacentCells(cell,dim,cellfacemap,owners,neighbours);
      for (Index i: adjacentCells) {
         if (!checkCell(i,ghostCellCandidates,notGhostCells,dim,outerVertices,cellfacemap,faces,owners,neighbours))
            return false;
      }
      return true;
   } else {

      notGhostCells.insert(cell);
      return true;
   }

   return false;
}

bool ReadFOAM::buildGhostCells(int processor, GhostMode mode) {
   //std::cerr << "buildGhostCells(p=" << processor << ", mode=" << mode << ")" << std::endl;
   auto &boundaries = *m_boundaries[processor];

   UnstructuredGrid::ptr grid = m_currentgrid[processor];
   auto &el = grid->el();
   auto &cl = grid->cl();
   auto &tl = grid->tl();
   auto &x = grid->x();
   auto &y = grid->y();
   auto &z = grid->z();

   for (const auto &b :boundaries.procboundaries) {
      Index neighborProc=b.neighborProc;
      std::shared_ptr<GhostCells> out(new GhostCells()); //object that will be sent to neighbor processor
      m_GhostCellsOut[processor][neighborProc] = out;
      std::vector<Index> &procBoundaryVertices = m_procBoundaryVertices[processor][neighborProc];
      std::unordered_set<Index> &procGhostCellCandidates = m_procGhostCellCandidates[processor][neighborProc];

      std::vector<Index> &elOut = out->el;
      std::vector<SIndex> &clOut = out->cl;
      std::vector<unsigned char> &tlOut = out->tl;
      std::vector<Scalar> &pointsOutX = out->x;
      std::vector<Scalar> &pointsOutY = out->y;
      std::vector<Scalar> &pointsOutZ = out->z;

      Index coordCount = 0;
      if (mode == ALL || mode == BASE) { //create ghost cell topology and send vertice-coordinates
         //build ghost cell element list and connectivity list for current boundary patch
         elOut.push_back(0);
         Index conncount=0;
         for (const Index cell: procGhostCellCandidates) {
            Index elementStart = el[cell];
            Index elementEnd = el[cell + 1];
            for (Index j=elementStart; j<elementEnd; ++j) {
               SIndex point = cl[j];
               clOut.push_back(point);
               ++conncount;
            }
            elOut.push_back(conncount);
            tlOut.push_back(tl[cell]);
         }
         //Create Vertices Mapping
         std::map<Index, SIndex> verticesMapping;
         //shared vertices (coords do not have to be sent) -> mapped to negative values
         SIndex c=-1;
         for (const Index v: procBoundaryVertices) {
            if (verticesMapping.emplace(v,c).second) {//emplace tries to create the map entry with the key/value pair (v,c) - if an entry with the key v  already exists it does not insert the new one
               --c;                                   //and returns a pair which consosts of a pointer to the already existing key/value pair (first) and a boolean that states if anything was inserted into the map (second)
            }
         }
         //vertices with coordinates that have to be sent -> mapped to positive values
         for (const SIndex v: clOut) {
            if (verticesMapping.emplace(v,coordCount).second) {
               ++coordCount;
            }
         }

         //Change connectivity list entries to the mapped values
         for (SIndex &v: clOut) {
            v = verticesMapping[v];
         }

         //save the vertices mapping for later use
         m_verticesMappings[processor][neighborProc]=verticesMapping;
      } else if (mode == COORDS) {
          const std::map<Index, SIndex> &verticesMapping = m_verticesMappings[processor][neighborProc];
          for (const auto &v: verticesMapping) {
              SIndex s=v.second;
              if (s >= 0) {
                  ++coordCount;
              }
          }
      }

      if (mode == ALL || mode == COORDS) {
          //create and fill Coordinate Vectors that have to be sent
          pointsOutX.resize(coordCount);
          pointsOutY.resize(coordCount);
          pointsOutZ.resize(coordCount);

          const std::map<Index, SIndex> &verticesMapping = m_verticesMappings[processor][neighborProc];
          for (const auto &v: verticesMapping) {
              Index f=v.first;
              SIndex s=v.second;
              if (s >= 0) {
                  pointsOutX[s] = x[f];
                  pointsOutY[s] = y[f];
                  pointsOutZ[s] = z[f];
              }
          }
      }
   }

   //build requests for Ghost Cells
   for (const auto &b :boundaries.procboundaries) {
      Index neighborProc=b.neighborProc;
      int myRank=rank();
      int neighborRank = rankForBlock(neighborProc);
      std::shared_ptr<GhostCells> out = m_GhostCellsOut[processor][neighborProc];
      if (myRank != neighborRank) {
         m_requests.push_back(comm().isend(neighborRank, tag(processor,neighborProc), *out));
         std::shared_ptr<GhostCells> in(new GhostCells());
         m_GhostCellsIn[processor][neighborProc] = in;
         m_requests.push_back(comm().irecv(neighborRank, tag(neighborProc,processor), *in));
      } else {
         m_GhostCellsIn[neighborProc][processor] = out;
      }
   }
   return true;
}

bool ReadFOAM::buildGhostCellData(int processor) {
   //std::cerr << "buildGhostCellsData(p=" << processor << ")" << std::endl;
   vassert(m_buildGhost);
   auto &boundaries = *m_boundaries[processor];
   for (const auto &b :boundaries.procboundaries) {
      Index neighborProc=b.neighborProc;
      std::unordered_set<Index> &procGhostCellCandidates = m_procGhostCellCandidates[processor][neighborProc];
      for (int i = 0; i < NumPorts; ++i) {
         auto f=m_currentvolumedata[processor].find(i);
         if (f == m_currentvolumedata[processor].end()) {
            continue;
         }
         Object::ptr obj = f->second;
         Vec<Scalar, 1>::ptr v1 = Vec<Scalar, 1>::as(obj);
         Vec<Scalar, 3>::ptr v3 = Vec<Scalar, 3>::as(obj);
         if (!v1 && !v3) {
            std::cerr << "Could not send Data - unsupported Data-Object Type" << std::endl;
            continue;
         }
         if (v1) {
            std::shared_ptr<GhostData> dataOut(new GhostData(1));
            m_GhostDataOut[processor][neighborProc][i] = dataOut;
            auto &d=v1->x(0);
            for (const Index& cell: procGhostCellCandidates) {
               (*dataOut).x[0].push_back(d[cell]);
            }
         } else if (v3) {
            std::shared_ptr<GhostData> dataOut(new GhostData(3));
            m_GhostDataOut[processor][neighborProc][i] = dataOut;
            auto &d1=v3->x(0);
            auto &d2=v3->x(1);
            auto &d3=v3->x(2);
            for (const Index& cell: procGhostCellCandidates) {
               (*dataOut).x[0].push_back(d1[cell]);
               (*dataOut).x[1].push_back(d2[cell]);
               (*dataOut).x[2].push_back(d3[cell]);
            }
         }
      }
   }

   //build requests for Ghost Data
   for (const auto &b :boundaries.procboundaries) {
      Index neighborProc=b.neighborProc;
      int myRank=rank();
      int neighborRank = rankForBlock(neighborProc);

      std::map<int, std::shared_ptr<GhostData> > &m = m_GhostDataOut[processor][neighborProc];
      for (Index i=0; i<NumPorts; ++i) {
         if (m.find(i) != m.end()) {
            std::shared_ptr<GhostData> dataOut = m[i];
            if (myRank != neighborRank) {
               m_requests.push_back(comm().isend(neighborRank, tag(processor,neighborProc,i+1), *dataOut));
               std::shared_ptr<GhostData> dataIn(new GhostData((*dataOut).dim));
               m_GhostDataIn[processor][neighborProc][i] = dataIn;
               m_requests.push_back(comm().irecv(neighborRank, tag(neighborProc,processor,i+1), *dataIn));
            } else {
               m_GhostDataIn[neighborProc][processor][i] = dataOut;
            }
         }
      }
   }
   return true;
}

void ReadFOAM::processAllRequests() {
   vassert(m_buildGhost);
   mpi::wait_all(m_requests.begin(), m_requests.end());
   m_requests.clear();
   m_GhostCellsOut.clear();
   m_GhostDataOut.clear();
}

void ReadFOAM::applyGhostCells(int processor, GhostMode mode) {
   vassert(m_buildGhost);
   auto &boundaries = *m_boundaries[processor];
   UnstructuredGrid::ptr grid = m_currentgrid[processor];
   auto &el = grid->el();
   auto &cl = grid->cl();
   auto &tl = grid->tl();
   auto &x = grid->x();
   auto &y = grid->y();
   auto &z = grid->z();
   //std::cerr << "applyGhostCells(p=" << processor << ", mode=" << mode << "), #cells=" << el.size() << ", #coords: " << x.size() << std::endl;

   for (const auto &b: boundaries.procboundaries) {
       Index neighborProc=b.neighborProc;
       std::vector<Index> &procBoundaryVertices = m_procBoundaryVertices[processor][neighborProc];
       std::vector<Index> sharedVerticesMapping;
       for (const Index &v: procBoundaryVertices) {
          // create sharedVerticesMapping vector that consists of all the shared (already known) vertices
          // without duplicates and in order of "first appearance" when going through the boundary cells
          if(std::find(sharedVerticesMapping.begin(), sharedVerticesMapping.end(), v) == sharedVerticesMapping.end()) {
              sharedVerticesMapping.push_back(v);
          }
       }
       std::shared_ptr<GhostCells> in = m_GhostCellsIn[processor][neighborProc];
       std::vector<Index> &elIn = in->el;
       std::vector<SIndex> &clIn = in->cl;
       std::vector<unsigned char> &tlIn = in->tl;
       std::vector<Scalar> &pointsInX = in->x;
       std::vector<Scalar> &pointsInY = in->y;
       std::vector<Scalar> &pointsInZ = in->z;
       Index pointsSize=x.size();
       if (pointsSize == 0) {
           std::cerr << "Warning: no coordinates loaded yet" << std::endl;
       }

       if (mode == ALL || mode == BASE) { //ghost cell topology is unnknown and has to be appended to the current topology
          for (Index cell = 0; cell < tlIn.size();++cell) {//append new topology to old grid
             Index elementStart = elIn[cell];
             Index elementEnd = elIn[cell + 1];
             for (Index i = elementStart; i < elementEnd; ++i) {
                SIndex point = clIn[i];
                if (point < 0) {//if point<0 then vertice is already known and can be looked up in sharedVerticesMapping
                   point=sharedVerticesMapping[-point-1];
                } else {//else the vertice is unknown and its coordinates will be appended (in order of first appearance) to the old coord-lists so we point to an index beyond the current size
                   point+=pointsSize;
                }
                cl.push_back(point);
             }
             el.push_back(cl.size());
             tl.push_back(tlIn[cell]|UnstructuredGrid::GHOST_BIT);
          }
       }

       if (mode == ALL || mode == COORDS) {
           //append new coordinates to old coordinate-lists
           std::copy(pointsInX.begin(), pointsInX.end(), std::back_inserter(x));
           std::copy(pointsInY.begin(), pointsInY.end(), std::back_inserter(y));
           std::copy(pointsInZ.begin(), pointsInZ.end(), std::back_inserter(z));
       }
   }
   m_GhostCellsIn[processor].clear();
   m_procBoundaryVertices[processor].clear();
   //std::cerr << "applyGhostCells(p=" << processor << ", mode=" << mode << "), final #cells=" << el.size() << ", #coords: " << x.size() << std::endl;
}

void ReadFOAM::applyGhostCellsData(int processor) {
   //std::cerr << "applyGhostCellsData(p=" << processor << ")" << std::endl;
   auto &boundaries = *m_boundaries[processor];

   for (const auto &b :boundaries.procboundaries) {
      Index neighborProc=b.neighborProc;
      for (int i = 0; i < NumPorts; ++i) {
         auto f=m_currentvolumedata[processor].find(i);
         if (f == m_currentvolumedata[processor].end()) {
            continue;
         }
         Object::ptr obj = f->second;
         Vec<Scalar, 1>::ptr v1 = Vec<Scalar, 1>::as(obj);
         Vec<Scalar, 3>::ptr v3 = Vec<Scalar, 3>::as(obj);
         if (!v1 && !v3) {
            std::cerr << "Could not apply Data - unsupported Data-Object Type" << std::endl;
            continue;
         }
         //append ghost cell data to old data objects
         if (v1) {
            std::shared_ptr<GhostData> dataIn = m_GhostDataIn[processor][neighborProc][i];
            auto &d=v1->x(0);
            std::vector<Scalar> &x=(*dataIn).x[0];
            std::copy(x.begin(), x.end(), std::back_inserter(d));
         } else if (v3) {
            std::shared_ptr<GhostData> dataIn = m_GhostDataIn[processor][neighborProc][i];
            for (Index j=0; j<3; ++j) {
               auto &d=v3->x(j);
               std::vector<Scalar> &x=(*dataIn).x[j];
               std::copy(x.begin(), x.end(), std::back_inserter(d));
            }
         }
      }
   }
   m_GhostDataIn[processor].clear();
}

bool ReadFOAM::addGridToPorts(int processor) {
   for (auto &poly: m_currentbound[processor])
       addObject(m_boundOut, poly);
   return true;
}

bool ReadFOAM::addVolumeDataToPorts(int processor) {

    for (int portnum=0; portnum<NumPorts; ++portnum) {
        auto &volumedata = m_currentvolumedata[processor];
        if (volumedata.find(portnum) != volumedata.end()) {
            volumedata[portnum]->setGrid(m_currentgrid[processor]);
            volumedata[portnum]->setMapping(DataBase::Element);
            addObject(m_volumeDataOut[portnum], volumedata[portnum]);
        } else {
            addObject(m_volumeDataOut[portnum], m_currentgrid[processor]);
        }
    }
   return true;
}

bool ReadFOAM::readConstant(const std::string &casedir)
{
   for (int i=-1; i<m_case.numblocks; ++i) {
      if (rankForBlock(i) == rank()) {
         if (!readDirectory(casedir, i, -1))
            return false;
      }
   }

   if (m_case.varyingCoords && m_case.varyingGrid)
      return true;

   bool readGrid = m_readGrid->getValue();
   if (!readGrid)
      return true;

   GhostMode buildGhostMode = m_case.varyingCoords ? BASE : ALL;

   if (m_buildGhost) {

      for (int i=0; i<m_case.numblocks; ++i) {
         if (rankForBlock(i) == rank()) {
            if (!buildGhostCells(i, buildGhostMode))
               return false;
         }
      }

      processAllRequests();
   }

   for (int i=-1; i<m_case.numblocks; ++i) {
      if (rankForBlock(i) == rank()) {
         if (m_buildGhost) {
            applyGhostCells(i,buildGhostMode);
         }
         if (!m_case.varyingCoords)
            addGridToPorts(i);
      }
   }

   m_GhostCellsIn.clear();
   m_currentvolumedata.clear();
   return true;
}

bool ReadFOAM::readTime(const std::string &casedir, int timestep) {
   for (int i=-1; i<m_case.numblocks; ++i) {
      if (rankForBlock(i) == rank()) {
         if (!readDirectory(casedir, i, timestep))
            return false;
      }
   }

   bool readGrid = m_readGrid->getValue();
   if (!readGrid)
      return true;

   if (m_case.varyingCoords || m_case.varyingGrid) {
      const GhostMode ghostMode = m_case.varyingGrid ? ALL : COORDS;
      if (m_buildGhost) {
         for (int i=0; i<m_case.numblocks; ++i) {
            if (rankForBlock(i) == rank()) {
               if (!buildGhostCells(i,ghostMode))
                  return false;
            }
         }

         processAllRequests();
      }

      for (int i=-1; i<m_case.numblocks; ++i) {
         if (rankForBlock(i) == rank()) {
            if (m_buildGhost) {
               applyGhostCells(i,ghostMode);
            }
            addGridToPorts(i);
         }
      }
      m_GhostCellsIn.clear();
   }

   if (m_buildGhost) {
      for (int i=0; i<m_case.numblocks; ++i) {
         if (rankForBlock(i) == rank()) {
            buildGhostCellData(i);
         }
      }

      processAllRequests();
   }

   for (int i=-1; i<m_case.numblocks; ++i) {
      if (rankForBlock(i) == rank()) {
         if (m_buildGhost) {
            applyGhostCellsData(i);
         }
         addVolumeDataToPorts(i);
      }
   }

   m_GhostDataIn.clear();
   m_currentvolumedata.clear();
   return true;
}

bool ReadFOAM::compute()     //Compute is called when Module is executed
{
   const std::string casedir = m_casedir->getValue();
   m_boundaryPatches.add(m_patchSelection->getValue());
   m_case = getCaseInfo(casedir);
   if (!m_case.valid) {
      std::cerr << casedir << " is not a valid OpenFOAM case" << std::endl;
      return true;
   }

   if (!checkPolyMeshDirContent(m_case)) {
      std::cerr << "failed to gather topology directories for " << casedir << std::endl;
      return true;
   }

   m_buildGhost = m_buildGhostcellsParam->getValue() && m_case.numblocks>0;

   std::cerr << "# processors: " << m_case.numblocks << std::endl;
   std::cerr << "# time steps: " << m_case.timedirs.size() << std::endl;
   std::cerr << "grid topology: " << (m_case.varyingGrid?"varying":"constant") << std::endl;
   std::cerr << "grid coordinates: " << (m_case.varyingCoords?"varying":"constant") << std::endl;

   if (!readConstant(casedir)) {
      std::cerr << "reading of constant data failed" << std::endl;
      return true;
   }
   int skipfactor = m_timeskip->getValue()+1;
   for (Index timestep=0; timestep<m_case.timedirs.size()/skipfactor; ++timestep) {
      if (!readTime(casedir, timestep)) {
         std::cerr << "reading of data for timestep " << timestep << " failed" << std::endl;
      }
   }
   vassert(m_requests.empty());
   vassert(m_GhostCellsOut.empty());
   vassert(m_GhostDataOut.empty());
   vassert(m_GhostCellsIn.empty());
   vassert(m_GhostDataIn.empty());
   m_currentbound.clear();
   m_currentgrid.clear();
   m_basedir.clear();
   m_currentvolumedata.clear();
   m_owners.clear();
   m_boundaries.clear();
   m_procBoundaryVertices.clear();
   m_procGhostCellCandidates.clear();
   m_verticesMappings.clear();

   return true;
}

MODULE_MAIN(ReadFOAM)
