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
  const topology::FieldBase::DiscretizeInfo defaultInfo = {-1, -1, true};
  _auxFieldsFEInfo["default"] = defaultInfo;
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
pylith::materials::MaterialNew::initialize(const pylith::topology::Field& solution)
{ // initialize
  PYLITH_METHOD_BEGIN;

  // Get cells associated with material
  const pylith::topology::Mesh& mesh = solution.mesh();
  PetscDM dmMesh = mesh.dmMesh();assert(dmMesh);
  const bool includeOnlyCells = true;
  delete _materialIS; _materialIS = new topology::StratumIS(dmMesh, "material-id", _id, includeOnlyCells);assert(_materialIS);

  delete _auxFields; _auxFields = new topology::Field(mesh);assert(_auxFields);
  delete _auxFieldsQuery; _auxFieldsQuery = new topology::FieldQuery(*_auxFields);assert(_auxFieldsQuery);
  _auxFields->label("auxiliary fields");
  _auxFieldsSetup();
  _auxFields->subfieldsSetup();
  _auxFields->allocate();
  _auxFields->zeroAll();

  if (_auxFieldsDB) {
    assert(_normalizer);
    _auxFieldsQuery->openDB(_auxFieldsDB, _normalizer->lengthScale());
    _auxFieldsQuery->queryDB();
    _auxFieldsQuery->closeDB(_auxFieldsDB);
  } else { // else
    assert(0);
    throw std::logic_error("Unknown case for setting up auxiliary fields.");
  } // if/else

  // Set finite-element kernels
  _setFEKernels(solution);

  PYLITH_METHOD_END;
} // initialize


// ----------------------------------------------------------------------
// Set discretization information for subfield.
void
pylith::materials::MaterialNew::discretization(const char* name,
					       const pylith::topology::FieldBase::DiscretizeInfo& feInfo)
{ // discretization
  _auxFieldsFEInfo[name] = feInfo;
} // discretization


// ----------------------------------------------------------------------
// Get discretization information for subfield.
const pylith::topology::FieldBase::DiscretizeInfo& 
pylith::materials::MaterialNew::discretization(const char* name) const
{ // discretization
  PYLITH_METHOD_BEGIN;

  discretizations_type::const_iterator iter = _auxFieldsFEInfo.find(name);
  if (iter != _auxFieldsFEInfo.end()) {
    PYLITH_METHOD_RETURN(iter->second);
  } else { // not found so try default
    iter = _auxFieldsFEInfo.find("default");
    if (iter == _auxFieldsFEInfo.end()) {
      throw std::logic_error("Default discretization not set for auxiliary fields.");
    } // if
  } // if/else

  PYLITH_METHOD_RETURN(iter->second); // default
} // discretization


// End of file 
