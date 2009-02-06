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

#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys

#include <Selection.hh> // USES ALE::Selection

// ----------------------------------------------------------------------
// Default constructor
template<typename mesh_type>
pylith::topology::SubMesh<mesh_type>::SubMesh(void) :
  _coordsys(0),
  _debug(false)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Constructor with mesh and label for vertices marking boundary.
template<typename mesh_type>
pylith::topology::SubMesh<mesh_type>::SubMesh(const mesh_type& mesh,
					      const char* label) :
  _coordsys(0),
  _debug(false)
{ // constructor
  createSubMesh(mesh, label);
} // constructor

// ----------------------------------------------------------------------
// Default destructor
template<typename mesh_type>
pylith::topology::SubMesh<mesh_type>::~SubMesh(void)
{ // destructor
  delete _coordsys; _coordsys = 0;
} // destructor

// ----------------------------------------------------------------------
// Create Sieve mesh.
template<typename mesh_type>
void
pylith::topology::SubMesh<mesh_type>::createSubMesh(const mesh_type& mesh,
						    const char* label)
{ // createSieveMesh
  _mesh.destroy();

  const ALE::Obj<DomainSieveMesh>& meshSieveMesh = mesh.sieveMesh();
  assert(!meshSieveMesh.isNull());

  const ALE::Obj<IntSection>& groupField = meshSieveMesh->getIntSection(label);
  if (groupField.isNull()) {
    std::ostringstream msg;
    msg << "Could not find group of points '" << label << "' in mesh.";
    throw std::runtime_error(msg.str());
  } // if
#if 0
  // QUESTION FOR MATT
  // Why doesn't this work?
  // SieveMesh and SieveSubMesh are typedefs in SubMesh.hh
  _mesh = 
    ALE::Selection<DomainSieveMesh>::submeshV<SieveMesh>(meshSieveMesh,
							 groupField);
#else
  _mesh = 
    ALE::Selection<ALE::IMesh<> >::submeshV<SieveMesh>(meshSieveMesh,
						       groupField);
#endif
  if (_mesh.isNull()) {
    std::ostringstream msg;
    msg << "Could not construct boundary mesh for boundary '"
	<< label << "'.";
    throw std::runtime_error(msg.str());
  } // if
  _mesh->setRealSection("coordinates", 
			meshSieveMesh->getRealSection("coordinates"));

  // Create the parallel overlap
  ALE::Obj<typename SieveMesh::send_overlap_type> sendParallelMeshOverlap =
    _mesh->getSendOverlap();
  ALE::Obj<typename SieveMesh::recv_overlap_type> recvParallelMeshOverlap =
    _mesh->getRecvOverlap();
  typename DomainSieveMesh::renumbering_type& renumbering = 
    meshSieveMesh->getRenumbering();
  //   Can I figure this out in a nicer way?
  ALE::SetFromMap<std::map<typename DomainSieveMesh::point_type,
    typename DomainSieveMesh::point_type> > globalPoints(renumbering);

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
template<typename mesh_type>
void 
pylith::topology::SubMesh<mesh_type>::initialize(void)
{ // initialize
  if (0 != _coordsys)
    _coordsys->initialize();
} // initialize


// End of file 
