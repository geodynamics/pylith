#ifndef included_ALE_Mesh_hh
#define included_ALE_Mesh_hh

#ifndef  included_ALE_CoSifter_hh
#include <CoSifter.hh>
#endif
#ifndef  included_ALE_ParDelta_hh
#include <ParDelta.hh>
#endif
#ifndef  included_ALE_Partitioner_hh
#include <Partitioner.hh>
#endif

#include <petscvec.h>

namespace ALE {
  class PetscMesh;
} // namespace ALE

class ALE::PetscMesh {
public:
  typedef ALE::Point point_type;
  typedef std::vector<point_type> PointArray;
  typedef ALE::Sieve<point_type,int,int> sieve_type;
  typedef point_type patch_type;
  typedef CoSifter<sieve_type, patch_type, point_type, int> bundle_type;
  typedef CoSifter<sieve_type, patch_type, point_type, double> field_type;
  typedef CoSifter<sieve_type, ALE::pair<patch_type,int>, point_type, double> foliation_type;
  typedef std::map<std::string, Obj<field_type> > FieldContainer;
  typedef std::map<int, Obj<bundle_type> > BundleContainer;
  int debug;

private:
  Obj<sieve_type> topology;
  Obj<field_type> coordinates;
  Obj<field_type> boundary;
  Obj<foliation_type> boundaries;
  FieldContainer  fields;
  BundleContainer bundles;
  MPI_Comm        _comm;
  int             _commRank;
  int             _commSize;
  int             dim;

  //FIX:
public:
  bool            distributed;

public:
  
  PetscMesh(MPI_Comm comm, 
       int dimension, 
       int debug = 0);

  MPI_Comm comm() const {return this->_comm;};
  void setComm(MPI_Comm comm) {
    this->_comm = comm;
    MPI_Comm_rank(comm, &this->_commRank);
    MPI_Comm_size(comm, &this->_commSize);
  };
  int             commRank() const {return this->_commRank;};
  int             commSize() const {return this->_commSize;};
  Obj<sieve_type> getTopology() const {return this->topology;};
  void            setTopology(const Obj<sieve_type>& topology) {this->topology = topology;};
  int             getDimension() const {return this->dim;};
  void            setDimension(int dim) {this->dim = dim;};
  Obj<field_type> getCoordinates() const {return this->coordinates;};
  void            setCoordinates(const Obj<field_type>& coordinates) {this->coordinates = coordinates;};
  Obj<field_type> getBoundary() const {return this->boundary;};
  void            setBoundary(const Obj<field_type>& boundary) {this->boundary = boundary;};
  Obj<foliation_type> getBoundaries() const {return this->boundaries;};

  Obj<bundle_type> getBundle(const int dim);

  Obj<field_type> getField(const std::string& name);

  bool hasField(const std::string& name);
  
  Obj<std::set<std::string> > getFields();

  void buildHexFaces(int dim, 
		     std::map<int, int*> *curSimplex, 
		     PointArray *boundary, 
		     point_type& simplex);

  void buildFaces(int dim, 
		  std::map<int, int*> *curSimplex, 
		  PointArray *boundary, 
		  point_type& simplex);

  void buildTopology(int numSimplices, 
		     int simplices[], 
		     int numVertices, 
		     bool interpolate = true, 
		     int corners = -1);

  void createVertexBundle(int numSimplices, 
			  int simplices[], 
			  int elementOffset = 0, 
			  int corners = -1);

  void createSerialCoordinates(int embedDim, 
			       int numSimplices, 
			       double coords[]);

  void createParallelVertexReorder(Obj<bundle_type> serialVertexBundle);

  void populate(int numSimplices, 
		int simplices[], 
		int numVertices, 
		double coords[], 
		bool interpolate = true, 
		int corners = -1);

  void populateBd(int numSimplices, 
		  int simplices[], 
		  int numVertices, 
		  double coords[], 
		  bool interpolate = true, 
		  int corners = -1);


}; // PetscMesh

#endif
