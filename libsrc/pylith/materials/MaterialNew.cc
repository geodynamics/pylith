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
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "MaterialNew.hh" // implementation of object methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // HOLDSA Field
#include "pylith/topology/FieldQuery.hh" // HOLDSA FieldQuery
#include "pylith/topology/CoordsVisitor.hh" // USES CoordsVisitor
#include "pylith/topology/VisitorMesh.hh" // USES VecVisitorMesh
#include "pylith/topology/Stratum.hh" // USES StratumIS
#include "pylith/utils/array.hh" // USES scalar_array, std::vector

#include "spatialdata/spatialdb/SpatialDB.hh" // USES SpatialDB
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include <petscds.h> // USES PetscDS

#include <strings.h> // USES strcasecmp()
#include <cassert> // USES assert()
#include <stdexcept> // USES std::runtime_error
#include <sstream> // USES std::ostringstream

// ----------------------------------------------------------------------
// Default constructor.
pylith::materials::MaterialNew::MaterialNew(const int dimension) :
  _materialIS(NULL),
  _auxFieldsDB(NULL),
  _auxFieldsQuery(NULL),
  _dimension(dimension),
  _id(0),
  _label("")
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::materials::MaterialNew::~MaterialNew(void)
{ // destructor
  deallocate();
} // destructor

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::materials::MaterialNew::deallocate(void)
{ // deallocate
  PYLITH_METHOD_BEGIN;

  delete _materialIS; _materialIS = NULL;
  delete _auxFieldsQuery; _auxFieldsQuery = NULL;

  _auxFieldsDB = 0; // :TODO: Use shared pointer.

  PYLITH_METHOD_END;
} // deallocate
  
// ----------------------------------------------------------------------
// Get physical property parameters and initial state (if used) from database.
void
pylith::materials::MaterialNew::initialize(const topology::Mesh& mesh)
{ // initialize
  PYLITH_METHOD_BEGIN;

  // Get cells associated with material
  PetscDM dmMesh = mesh.dmMesh();assert(dmMesh);
  const bool includeOnlyCells = true;
  delete _materialIS; _materialIS = new topology::StratumIS(dmMesh, "material-id", _id, includeOnlyCells);assert(_materialIS);
  //const PetscInt numCells = _materialIS->size();
  //const PetscInt* cells = _materialIS->points();

  delete _auxFields; _auxFields = new topology::Field(mesh);assert(_auxFields);
  _auxFields->label("auxiliary fields");

  // :TODO: Update setup of auxiliary field
#if 0
  int_array cellsTmp(cells, numCells);
  _auxFields->newSection(cellsTmp, propsFiberDim);
  _auxFields->allocate();
  _auxFields->zeroAll();
#endif

  if (_auxFieldsDB) {
    assert(_normalizer);
    _auxFieldsQuery->queryDB(_auxFields, _auxFieldsDB, _normalizer->lengthScale());
  } else { // else
    assert(0);
    throw std::logic_error("Unknown case for setting up auxiliary fields.");
  } // if/else

  PYLITH_METHOD_END;
} // initialize


// End of file 
