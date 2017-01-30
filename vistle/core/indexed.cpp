#include "indexed.h"
#include "celltree_impl.h"
#include "assert.h"
#include "archives.h"

namespace vistle {

Indexed::Indexed(const Index numElements, const Index numCorners,
                      const Index numVertices,
                      const Meta &meta)
   : Indexed::Base(static_cast<Data *>(NULL))
{
    refreshImpl();
}

bool Indexed::isEmpty() const {

   return getNumElements()==0 || getNumCorners()==0;
}

bool Indexed::checkImpl() const {

   V_CHECK (d()->el->check());
   V_CHECK (d()->cl->check());
   V_CHECK (d()->el->size() > 0);
   V_CHECK (el()[0] == 0);
   if (getNumElements() > 0) {
      V_CHECK (el()[getNumElements()-1] < getNumCorners());
      V_CHECK (el()[getNumElements()] == getNumCorners());
   }

   if (getNumCorners() > 0) {
      V_CHECK (cl()[0] < getNumVertices());
      V_CHECK (cl()[getNumCorners()-1] < getNumVertices());
   }

   return true;
}

std::pair<Vector, Vector> Indexed::getBounds() const {
   if (hasCelltree()) {
      const auto ct = getCelltree();
      return std::make_pair(Vector(ct->min()), Vector(ct->max()));
   }

   return Base::getMinMax();
}

bool Indexed::hasCelltree() const {
   if (m_celltree)
       return true;

   return hasAttachment("celltree");
}

Indexed::Celltree::const_ptr Indexed::getCelltree() const {

   if (m_celltree)
       return m_celltree;

   boost::interprocess::scoped_lock<boost::interprocess::interprocess_recursive_mutex> lock(d()->attachment_mutex);
   if (!hasAttachment("celltree")) {
      refresh();
      createCelltree(getNumElements(), &el()[0], &cl()[0]);
   }

   m_celltree = Celltree::as(getAttachment("celltree"));
   vassert(m_celltree);
   return m_celltree;
}

void Indexed::createCelltree(Index nelem, const Index *el, const Index *cl) const {

   if (hasCelltree())
      return;

   const Scalar *coords[3] = {
      &x()[0],
      &y()[0],
      &z()[0]
   };
   const Scalar smax = std::numeric_limits<Scalar>::max();
   Vector vmin, vmax;
   vmin.fill(-smax);
   vmax.fill(smax);

   std::vector<Vector> min(nelem, vmax);
   std::vector<Vector> max(nelem, vmin);

   Vector gmin=vmax, gmax=vmin;
   for (Index i=0; i<nelem; ++i) {
      const Index start = el[i], end = el[i+1];
      for (Index c = start; c<end; ++c) {
         const Index v = cl[c];
         for (int d=0; d<3; ++d) {
            min[i][d] = std::min(min[i][d], coords[d][v]);
#if 0
            if (min[i][d] > coords[d][v]) {
               min[i][d] = coords[d][v];
               if (gmin[d] > min[i][d])
                  gmin[d] = min[i][d];
            }
#endif
            max[i][d] = std::max(max[i][d], coords[d][v]);
#if 0
            if (max[i][d] < coords[d][v]) {
               max[i][d] = coords[d][v];
               if (gmax[d] < max[i][d])
                  gmax[d] = max[i][d];
            }
#endif
         }
      }
      for (int d=0; d<3; ++d) {
         gmin[d] = std::min(gmin[d], min[i][d]);
         gmax[d] = std::max(gmax[d], max[i][d]);
      }
   }

   typename Celltree::ptr ct(new Celltree(nelem));
   ct->init(min.data(), max.data(), gmin, gmax);
   addAttachment("celltree", ct);
   if (!validateCelltree()) {
       std::cerr << "ERROR: Celltree validation failed." << std::endl;
   }
}


bool Indexed::hasVertexOwnerList() const {

   if (m_vertexOwnerList)
       return true;

   return hasAttachment("vertexownerlist");
}

Indexed::VertexOwnerList::const_ptr Indexed::getVertexOwnerList() const {

   if (m_vertexOwnerList)
       return m_vertexOwnerList;

   boost::interprocess::scoped_lock<boost::interprocess::interprocess_recursive_mutex> lock(d()->attachment_mutex);
   if (!hasAttachment("vertexownerlist")) {
      refresh();
      createVertexOwnerList();
   }

   m_vertexOwnerList = VertexOwnerList::as(getAttachment("vertexownerlist"));
   vassert(m_vertexOwnerList);
   return m_vertexOwnerList;
}

void Indexed::createVertexOwnerList() const {

   if (hasVertexOwnerList())
      return;

   Index numelem=getNumElements();
   Index numcoord=getNumVertices();

   VertexOwnerList::ptr vol(new VertexOwnerList(numcoord));
   const auto cl = &this->cl()[0];
   const auto el = &this->el()[0];
   auto vertexList=vol->vertexList().data();

   std::fill(vol->vertexList().begin(), vol->vertexList().end(), 0);
   std::vector<Index> used_vertex_list(numcoord, InvalidIndex);
   auto uvl = used_vertex_list.data();

   // Calculation of the number of cells that contain a certain vertex:
   // temporarily stored in vertexList
   for (Index i = 0; i < numelem; i++) {
      const Index begin = el[i], end = el[i+1];
      for (Index j = begin; j<end; ++j) {
         const Index v = cl[j];
         if (uvl[v] != j) {
             vertexList[v]++;
             uvl[v] = j;
         }
      }
   }

   // create the vertexList: prefix sum
   // vertexList will index into cellList
   std::vector<Index> outputIndex(numcoord);
   auto outIdx = outputIndex.data();
   {
      Index numEnt = 0;
      for (Index i = 0; i < numcoord; i++)
      {
         Index n = numEnt;
         numEnt += vertexList[i];
         vertexList[i] = n;
         outIdx[i] = n;
      }
      vertexList[numcoord] = numEnt;
      vol->cellList().resize(numEnt);
   }

   //fill the cellList
   std::fill(used_vertex_list.begin(), used_vertex_list.end(), InvalidIndex);
   auto cellList=vol->cellList().data();
   for (Index i = 0; i < numelem; i++) {
      const Index begin = el[i], end = el[i+1];
      for (Index j = begin; j<end; ++j) {
         const Index v = cl[j];
         if (uvl[v] != j) {
            uvl[v] = j;
            cellList[outIdx[v]] = i;
            outIdx[v]++;
         }
      }
   }

   addAttachment("vertexownerlist",vol);
}

void Indexed::removeVertexOwnerList() const {
    removeAttachment("vertexownerlist");
}

Indexed::NeighborFinder::NeighborFinder(const Indexed *indexed)
    : indexed(indexed)
{
    auto ol = indexed->getVertexOwnerList();
    numElem = indexed->getNumElements();
    numVert = indexed->getNumCorners();
    el = indexed->el();
    cl = indexed->cl();
    vl = ol->vertexList();
    vol = ol->cellList();
}

Index Indexed::NeighborFinder::getNeighborElement(Index elem, Index v1, Index v2, Index v3) const {

   if (v1 == v2 || v1 == v3 || v2 == v3) {
      std::cerr << "WARNING: getNeighborElement was not called with 3 unique vertices." << std::endl;
      return InvalidIndex;
   }

   const Index *elems = &vol[vl[v1]];
   Index nelem = vl[v1+1] - vl[v1];
   for (Index j=0; j<nelem; ++j) {
       Index e = elems[j];
       if (e == elem)
           continue;
       bool f2=false, f3=false;
       for (Index i=el[e]; i<el[e+1]; ++i) {
           const Index v = cl[i];
           if (v == v2) {
               f2 = true;
               if (f3)
                   return e;
           } else if (v == v3) {
               f3 = true;
               if (f2)
                   return e;
           }
       }
   }
   return InvalidIndex;
}

std::vector<Index> Indexed::NeighborFinder::getContainingElements(Index vert) const {

   const Index begin = vl[vert], end = vl[vert+1];
   return std::vector<Index>(&vol[begin], &vol[end]);
}

std::vector<Index> Indexed::NeighborFinder::getNeighborElements(Index elem) const {

   std::vector<Index> elems;
   if (elem == InvalidIndex)
       return elems;

   const auto vert = indexed->cellVertices(elem);
   for (auto v: vert) {
       const auto e = getContainingElements(v);
       std::copy_if(e.begin(), e.end(), std::back_inserter(elems), [elem](Index el) -> bool { return el!=elem; });
   }

   std::sort(elems.begin(), elems.end());
   auto last = std::unique(elems.begin(), elems.end());

   elems.resize(last - elems.begin());
   return elems;
}

const Indexed::NeighborFinder &Indexed::getNeighborFinder() const {

    if (!m_neighborfinder) {
        refresh();
        m_neighborfinder.reset(new NeighborFinder(this));
    }

    return *m_neighborfinder;
}

struct CellBoundsFunctor: public Indexed::Celltree::CellBoundsFunctor {

   CellBoundsFunctor(const Indexed *indexed)
      : m_indexed(indexed)
      {}

   bool operator()(Index elem, Vector *min, Vector *max) const {

      auto b = m_indexed->elementBounds(elem);
      *min = b.first;
      *max = b.second;
      return true;
   }

   const Indexed *m_indexed;
};

bool Indexed::validateCelltree() const {

   if (!hasCelltree())
      return false;

   CellBoundsFunctor boundFunc(this);
   auto ct = getCelltree();
   if (!ct->validateTree(boundFunc)) {
       std::cerr << "Indexed: Celltree validation failed with " << getNumElements() << " elements total, bounds: " << getBounds().first << "-" << getBounds().second << std::endl;
       return false;
   }
   return true;
}

std::pair<Vector, Vector> Indexed::elementBounds(Index elem) const {
   const Index *el = &this->el()[0];
   const Index *cl = &this->cl()[0];
   const Index begin = el[elem], end = el[elem+1];
   const Scalar *x[3] = { &this->x()[0], &this->y()[0], &this->z()[0] };

   const Scalar smax = std::numeric_limits<Scalar>::max();
   Vector min(smax, smax, smax), max(-smax, -smax, -smax);
   for (Index i=begin; i<end; ++i) {
       Index v=cl[i];
       for (int c=0; c<3; ++c) {
           min[c] = std::min(min[c], x[c][v]);
           max[c] = std::max(max[c], x[c][v]);
       }
   }
   return std::make_pair(min, max);
}

std::vector<Index> Indexed::cellVertices(Index elem) const {
    const Index *el = &this->el()[0];
    const Index *cl = &this->cl()[0];
    const Index begin = el[elem], end = el[elem+1];
    return std::vector<Index>(&cl[begin], &cl[end]);
}

void Indexed::refreshImpl() const {

    const Data *d = static_cast<Data *>(m_data);
    m_el = (d && d->el.valid()) ? d->el->data() : nullptr;
    m_cl = (d && d->cl.valid()) ? d->cl->data() : nullptr;
}

Indexed::Data::Data(const Index numElements, const Index numCorners,
             const Index numVertices,
             Type id, const std::string & name,
             const Meta &meta)
   : Indexed::Base::Data(numVertices, id, name,
         meta)
{
   el.construct(numElements+1);
   cl.construct(numCorners);
   (*el)[0] = 0;
}

Indexed::Data::Data(const Indexed::Data &o, const std::string &name)
: Indexed::Base::Data(o, name)
, el(o.el)
, cl(o.cl)
{
}

Indexed::Data *Indexed::Data::create(const std::string &objId, Type id,
            const Index numElements, const Index numCorners,
            const Index numVertices,
            const Meta &meta) {

   // required for boost::serialization
   assert("should never be called" == NULL);

   return NULL;
}


Index Indexed::getNumElements() const {

   return d()->el->size()-1;
}

Index Indexed::getNumCorners() const {

   return d()->cl->size();
}

Index Indexed::getNumVertices() const {

    return getSize();
}

Object::Type Indexed::type() {
   vassert("should not be called" == 0);
   return Object::UNKNOWN;
}

V_SERIALIZERS(Indexed);
V_OBJECT_CTOR(Indexed);

} // namespace vistle
