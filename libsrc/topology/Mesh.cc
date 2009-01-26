// -*- C++ -*-
//
// ======================================================================
//
//                           Brad T. Aagaard
//                        U.S. Geological Survey
//
// {LicenseText}
//
// ======================================================================
//

#include <portinfo>

#include "Mesh.hh" // implementation of class methods

#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys

// ----------------------------------------------------------------------
// Default constructor
pylith::topology::Mesh::Mesh(void) :
  _coordsys(0),
  _comm(PETSC_COMM_WORLD),
  _debug(false)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Default constructor
pylith::topology::Mesh::Mesh(const int dim,
			     const MPI_Comm& comm) :
  _mesh(new SieveMesh(comm, dim)),
  _coordsys(0),
  _comm(comm),
  _debug(false)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Default destructor
pylith::topology::Mesh::~Mesh(void)
{ // destructor
  delete _coordsys; _coordsys = 0;
} // destructor

// ----------------------------------------------------------------------
// Create Sieve mesh.
void
pylith::topology::Mesh::createSieveMesh(const int dim)
{ // createSieveMesh
  _mesh.destroy();
  _mesh = new SieveMesh(_comm, dim);
  _mesh->setDebug(_debug);
} // createSieveMesh

// ----------------------------------------------------------------------
// Set coordinate system.
void
pylith::topology::Mesh::coordsys(const spatialdata::geocoords::CoordSys* cs)
{ // coordsys
  delete _coordsys; _coordsys = (0 != cs) ? cs->clone() : 0;
} // coordsys

// ----------------------------------------------------------------------
// Initialize the finite-element mesh.
void 
pylith::topology::Mesh::initialize(void)
{ // initialize
  if (0 != _coordsys)
    _coordsys->initialize();
} // initialize

// ----------------------------------------------------------------------
// Create boundary mesh from mesh and label.
void
pylith::topology::Mesh::boundaryMeshFromMesh(const char* label)
{ // boundaryMeshFromMesh
#if 0
  assert(!_mesh.isNull());

  const ALE::Obj<SieveMesh::int_section_type>& groupField = 
    _mesh->getIntSection(_label);
  if (groupField.isNull()) {
    std::ostringstream msg;
    msg << "Could not find group of points '" << label << "' in mesh.";
    throw std::runtime_error(msg.str());
  } // if
  ALE::Obj<SubSieveMesh> boundarySieveMesh = 
    ALE::Selection<SieveMesh>::submeshV<SieveSubMesh>(_mesh, groupField);
  if (boundarySieveMesh.isNull()) {
    std::ostringstream msg;
    msg << "Could not construct boundary mesh for boundary '"
	<< label << "'.";
    throw std::runtime_error(msg.str());
  } // if
  boundarySieveMesh->setRealSection("coordinates", 
				    _mesh->getRealSection("coordinates"));

  // Create the parallel overlap
  ALE::Obj<SieveSubMesh::send_overlap_type> sendParallelMeshOverlap =
    boundarySieveMesh->getSendOverlap();
  ALE::Obj<SieveSubMesh::recv_overlap_type> recvParallelMeshOverlap =
    boundarySieveMesh->getRecvOverlap();
  SieveMesh::renumbering_type& renumbering = _mesh->getRenumbering();
  //   Can I figure this out in a nicer way?
  ALE::SetFromMap<std::map<SieveMesh::point_type,SieveMesh::point_type> > globalPoints(renumbering);

  ALE::OverlapBuilder<>::constructOverlap(globalPoints, renumbering,
					  sendParallelMeshOverlap,
					  recvParallelMeshOverlap);
  boundarySieveMesh->setCalculatedOverlap(true);

  Mesh boundaryMesh;

  return boundaryMesh
#endif
} // createBoundaryMesh


// End of file 
