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

#include "SubMesh.hh" // implementation of class methods

#include "Mesh.hh" // USES Mesh
#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys

#include <Selection.hh> // USES ALE::Selection

// ----------------------------------------------------------------------
// Default constructor
pylith::topology::SubMesh::SubMesh(void) :
  _coordsys(0),
  _debug(false)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Constructor with mesh and label for vertices marking boundary.
pylith::topology::SubMesh::SubMesh(const Mesh& mesh,
				   const char* label) :
  _coordsys(0),
  _debug(false)
{ // constructor
  createSubMesh(mesh, label);
} // constructor

// ----------------------------------------------------------------------
// Default destructor
pylith::topology::SubMesh::~SubMesh(void)
{ // destructor
  delete _coordsys; _coordsys = 0;
} // destructor

// ----------------------------------------------------------------------
// Create Sieve mesh.
void
pylith::topology::SubMesh::createSubMesh(const Mesh& mesh,
					 const char* label)
{ // createSieveMesh
  _mesh.destroy();

  const ALE::Obj<SieveMesh>& meshSieveMesh = mesh.sieveMesh();
  assert(!meshSieveMesh.isNull());

  const ALE::Obj<SieveMesh::int_section_type>& groupField = 
    meshSieveMesh->getIntSection(label);
  if (groupField.isNull()) {
    std::ostringstream msg;
    msg << "Could not find group of points '" << label << "' in mesh.";
    throw std::runtime_error(msg.str());
  } // if
  _mesh = 
    ALE::Selection<SieveMesh>::submeshV<SieveSubMesh>(meshSieveMesh,
						      groupField);
  if (_mesh.isNull()) {
    std::ostringstream msg;
    msg << "Could not construct boundary mesh for boundary '"
	<< label << "'.";
    throw std::runtime_error(msg.str());
  } // if
  _mesh->setRealSection("coordinates", 
			meshSieveMesh->getRealSection("coordinates"));

  // Create the parallel overlap
  ALE::Obj<SieveSubMesh::send_overlap_type> sendParallelMeshOverlap =
    _mesh->getSendOverlap();
  ALE::Obj<SieveSubMesh::recv_overlap_type> recvParallelMeshOverlap =
    _mesh->getRecvOverlap();
  SieveMesh::renumbering_type& renumbering = meshSieveMesh->getRenumbering();
  //   Can I figure this out in a nicer way?
  ALE::SetFromMap<std::map<SieveMesh::point_type,SieveMesh::point_type> > globalPoints(renumbering);

  ALE::OverlapBuilder<>::constructOverlap(globalPoints, renumbering,
					  sendParallelMeshOverlap,
					  recvParallelMeshOverlap);
  _mesh->setCalculatedOverlap(true);

  // Set data from mesh.
  _mesh->setDebug(mesh.debug());
  delete _coordsys; _coordsys = 0;
  const spatialdata::geocoords::CoordSys* cs = mesh.coordsys();
  if (0 != cs)
    _coordsys = cs->clone();
} // createSubMesh

// ----------------------------------------------------------------------
// Initialize the finite-element mesh.
void 
pylith::topology::SubMesh::initialize(void)
{ // initialize
  if (0 != _coordsys)
    _coordsys->initialize();
} // initialize


// End of file 
