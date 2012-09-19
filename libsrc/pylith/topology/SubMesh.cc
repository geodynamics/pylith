// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2012 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

#include <portinfo>

#include "SubMesh.hh" // implementation of class methods

#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys

#include <Selection.hh> // USES ALE::Selection

#include "pylith/utils/petscerror.h" // USES CHECK_PETSC_ERROR

#include <stdexcept> // USES std::runtime_error
#include <sstream> // USES std::ostringstream
#include <cassert> // USES assert()

// ----------------------------------------------------------------------
// Default constructor
pylith::topology::SubMesh::SubMesh(void) :
  _newMesh(PETSC_NULL),
  _coordsys(0),
  _debug(false)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Constructor with mesh and label for vertices marking boundary.
pylith::topology::SubMesh::SubMesh(const Mesh& mesh,
				   const char* label) :
  _newMesh(PETSC_NULL),
  _coordsys(0),
  _debug(false)
{ // constructor
  createSubMesh(mesh, label);
} // constructor

// ----------------------------------------------------------------------
// Default destructor
pylith::topology::SubMesh::~SubMesh(void)
{ // destructor
  deallocate();
} // destructor

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::topology::SubMesh::deallocate(void)
{ // deallocate
  delete _coordsys; _coordsys = 0;
  _mesh.destroy();
} // deallocate
  
// ----------------------------------------------------------------------
// Create Sieve mesh.
void
pylith::topology::SubMesh::createSubMesh(const Mesh& mesh,
					 const char* label)
{ // createSieveMesh
  _mesh.destroy();

  const ALE::Obj<DomainSieveMesh>& meshSieveMesh = mesh.sieveMesh();
  assert(!meshSieveMesh.isNull());
  DM dmMesh = mesh.dmMesh();
  assert(dmMesh);

  if (!meshSieveMesh->hasIntSection(label)) {
    std::ostringstream msg;
    msg << "Could not find group of points '" << label << "' in mesh.";
    throw std::runtime_error(msg.str());
  } // if
  PetscBool      hasLabel;
  PetscErrorCode err;
  err = DMComplexHasLabel(dmMesh, label, &hasLabel);CHECK_PETSC_ERROR(err);
  if (!hasLabel) {
    std::ostringstream msg;
    msg << "Could not find group of points '" << label << "' in DM mesh.";
    throw std::runtime_error(msg.str());
  } // if

  const ALE::Obj<IntSection>& groupField = meshSieveMesh->getIntSection(label);
  assert(!groupField.isNull());
  _mesh = 
    ALE::Selection<DomainSieveMesh>::submeshV<SieveMesh>(meshSieveMesh,
							 groupField);
  if (_mesh.isNull()) {
    std::ostringstream msg;
    msg << "Could not construct boundary mesh for boundary '"
	<< label << "'.";
    throw std::runtime_error(msg.str());
  } // if
  _mesh->setRealSection("coordinates", 
			meshSieveMesh->getRealSection("coordinates"));
  if (meshSieveMesh->hasRealSection("coordinates_dimensioned"))
    _mesh->setRealSection("coordinates_dimensioned", 
			  meshSieveMesh->getRealSection("coordinates_dimensioned"));

  /* TODO: Implement subMesh() for DMComplex */
  _newMesh = PETSC_NULL;

  // Create the parallel overlap
  const ALE::Obj<SieveMesh::sieve_type>& sieve = _mesh->getSieve();
  assert(!sieve.isNull());
  ALE::Obj<SieveMesh::send_overlap_type> sendParallelMeshOverlap =
    _mesh->getSendOverlap();
  assert(!sendParallelMeshOverlap.isNull());
  ALE::Obj<SieveMesh::recv_overlap_type> recvParallelMeshOverlap =
    _mesh->getRecvOverlap();
  assert(!recvParallelMeshOverlap.isNull());

  SieveMesh::renumbering_type& renumbering = _mesh->getRenumbering();

  DomainSieveMesh::renumbering_type& oldRenumbering = 
    meshSieveMesh->getRenumbering();
  const SieveMesh::renumbering_type::const_iterator oldBegin = 
    oldRenumbering.begin();
  const SieveMesh::renumbering_type::const_iterator oldEnd = 
    oldRenumbering.end();
  for (SieveMesh::renumbering_type::const_iterator r_iter = oldBegin;
       r_iter != oldEnd;
       ++r_iter)
    if (sieve->getChart().hasPoint(r_iter->second) && 
	(sieve->getConeSize(r_iter->second) || 
	 sieve->getSupportSize(r_iter->second)))
      renumbering[r_iter->first] = r_iter->second;
  
  //   Can I figure this out in a nicer way?
  ALE::SetFromMap<std::map<DomainSieveMesh::point_type,
    DomainSieveMesh::point_type> > globalPoints(renumbering);

  ALE::OverlapBuilder<>::constructOverlap(globalPoints, renumbering,
					  sendParallelMeshOverlap,
					  recvParallelMeshOverlap);
  _mesh->setCalculatedOverlap(true);

  // Set data from mesh.
  _mesh->setDebug(mesh.debug());
  coordsys(mesh);

  // Set name
  std::string meshLabel = "subdomain_" + std::string(label);
  _mesh->setName(meshLabel);

  int maxConeSizeLocal = sieve->getMaxConeSize();
  int maxConeSize = 0;
  err = MPI_Allreduce(&maxConeSizeLocal, &maxConeSize, 1, MPI_INT, MPI_MAX,
                      sieve->comm()); CHECK_PETSC_ERROR(err);

  if (maxConeSize <= 0) {
    std::ostringstream msg;
    msg << "Error while creating submesh. Submesh '" 
	<< label << "' does not contain any cells.\n"
	<< "Submeshes must be one dimension lower than the domain mesh.";
    throw std::runtime_error(msg.str());
  } // if
} // createSubMesh

// ----------------------------------------------------------------------
// Set coordinate system using mesh.
void
pylith::topology::SubMesh::coordsys(const Mesh& mesh)
{ // coordsys
  delete _coordsys; _coordsys = 0;
  const spatialdata::geocoords::CoordSys* cs = mesh.coordsys();
  if (0 != cs) {
    _coordsys = cs->clone();
    assert(0 != _coordsys);
    _coordsys->initialize();
  } // if
} // coordsys

// ----------------------------------------------------------------------
// Initialize the finite-element mesh.
void 
pylith::topology::SubMesh::initialize(void)
{ // initialize
  if (0 != _coordsys)
    _coordsys->initialize();
} // initialize


// End of file 
