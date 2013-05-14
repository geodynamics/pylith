// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2013 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "BCIntegratorSubMesh.hh" // implementation of object methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Fields.hh" // USES Fields
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/Stratum.hh" // USES Stratum

#include "pylith/feassemble/Quadrature.hh" // USES Quadrature

// ----------------------------------------------------------------------
// Default constructor.
pylith::bc::BCIntegratorSubMesh::BCIntegratorSubMesh(void) :
  _boundaryMesh(0),
  _parameters(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::bc::BCIntegratorSubMesh::~BCIntegratorSubMesh(void)
{ // destructor
  deallocate();
} // destructor

// ----------------------------------------------------------------------
// Deallocate data structures.
void
pylith::bc::BCIntegratorSubMesh::deallocate(void)
{ // deallocate
  PYLITH_METHOD_BEGIN;

  BoundaryCondition::deallocate();
  feassemble::Integrator::deallocate();

  delete _boundaryMesh; _boundaryMesh = 0;
  delete _parameters; _parameters = 0;

  PYLITH_METHOD_END;
} // deallocate
  
// ----------------------------------------------------------------------
// Get submesh associated with boundary condition.
void
pylith::bc::BCIntegratorSubMesh::createSubMesh(const topology::Mesh& mesh)
{ // createSubMesh
  PYLITH_METHOD_BEGIN;

  delete _boundaryMesh; _boundaryMesh = 0;
  delete _parameters; _parameters = 0;

  _boundaryMesh = new topology::Mesh(mesh, _label.c_str());
  assert(_boundaryMesh);

  PYLITH_METHOD_END;
} // createSubMesh

// ----------------------------------------------------------------------
// Verify configuration is acceptable.
void
pylith::bc::BCIntegratorSubMesh::verifyConfiguration(const topology::Mesh& mesh) const 
{ // verifyConfiguration
  PYLITH_METHOD_BEGIN;

  assert(_quadrature);
  assert(_boundaryMesh);

  BoundaryCondition::verifyConfiguration(mesh);

  // check compatibility of quadrature and boundary mesh
  if (_quadrature->cellDim() != _boundaryMesh->dimension()) {
    std::ostringstream msg;
    msg << "Quadrature is incompatible with cells for boundary condition '"
	<< _label << "'.\n"
	<< "Dimension of boundary mesh: " << _boundaryMesh->dimension()
	<< ", dimension of quadrature: " << _quadrature->cellDim()
	<< ".";
    throw std::runtime_error(msg.str());
  } // if
  const int numCorners = _quadrature->numBasis();

  // Get 'surface' cells (1 dimension lower than top-level cells)
  const PetscDM dmSubMesh = _boundaryMesh->dmMesh();assert(dmSubMesh);
  topology::Stratum cellsStratum(dmSubMesh, topology::Stratum::HEIGHT, 1);
  const PetscInt cStart = cellsStratum.begin();
  const PetscInt cEnd = cellsStratum.end();

  // Make sure surface cells are compatible with quadrature.
  PetscInt depth = 0;
  PetscErrorCode err = DMPlexGetDepth(dmSubMesh, &depth);PYLITH_CHECK_ERROR(err);
  for (PetscInt c = cStart; c < cEnd; ++c) {
    PetscInt cellNumCorners = 0;
    err = DMPlexGetConeSize(dmSubMesh, c, &cellNumCorners);PYLITH_CHECK_ERROR(err);
    if (numCorners != cellNumCorners) {
      std::ostringstream msg;
      msg << "Quadrature is incompatible with cell for boundary condition '"
	  << _label << "'.\n"
	  << "Cell " << c << " has " << cellNumCorners
	  << " vertices but quadrature reference cell has "
	  << numCorners << " vertices.";
      throw std::runtime_error(msg.str());
    } // if
  } // for

  PYLITH_METHOD_END;
} // verifyConfiguration


// End of file 
