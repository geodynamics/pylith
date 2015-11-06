// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, Rice University
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

#include <portinfo>

#include "IntegratorPointwise.hh" // implementation of class methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/Jacobian.hh" // USES Jacobian
#include "pylith/topology/Stratum.hh" // USES Stratum
#include "pylith/topology/VisitorMesh.hh" // USES VecVisitorMesh
#include "pylith/topology/CoordsVisitor.hh" // USES CoordsVisitor

#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include "pylith/utils/array.hh" // USES scalar_array
#include "pylith/utils/EventLogger.hh" // USES EventLogger

#include <strings.h> // USES strcasecmp()
#include <cassert> // USES assert()
#include <stdexcept> // USES std::runtime_error
#include <iostream> // USES std::cerr

// ----------------------------------------------------------------------
// Constructor
pylith::feassemble::IntegratorPointwise::IntegratorPointwise(void) :
  _auxFields(NULL)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::feassemble::IntegratorPointwise::~IntegratorPointwise(void)
{ // destructor
  deallocate();
} // destructor
  
// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::feassemble::IntegratorPointwise::deallocate(void)
{ // deallocate
  PYLITH_METHOD_BEGIN;

  delete _auxFields; _auxFields = 0;
  delete _auxFieldsQuery; _auxFieldsQuery = NULL;

  _auxFieldsDB = 0; // :TODO: Use shared pointer.

  PYLITH_METHOD_END;
} // deallocate

// ----------------------------------------------------------------------
// Return auxiliary fields for this problem
const pylith::topology::Field&
pylith::feassemble::IntegratorPointwise::auxFields(void) const
{ // auxFields
  PYLITH_METHOD_BEGIN;

  assert(_auxFields);

  PYLITH_METHOD_RETURN(*_auxFields);
} // auxFields

// ----------------------------------------------------------------------
// Check whether material has a given auxilirary field.
bool
pylith::feassemble::IntegratorPointwise::hasAuxField(const char* name)
{ // hasAuxField
  PYLITH_METHOD_BEGIN;

  assert(_auxFields);

  PYLITH_METHOD_RETURN(_auxFields->hasSubfield(name));
} // hasAuxField


// ----------------------------------------------------------------------
// Get auxiliary field.
void
pylith::feassemble::IntegratorPointwise::getAuxField(pylith::topology::Field *field,
						     const char* name) const
{ // getAuxField
  PYLITH_METHOD_BEGIN;

  assert(field);
  assert(_auxFields);

  field->copySubfield(*_auxFields, name);

  PYLITH_METHOD_END;
} // getAuxField

  
// ----------------------------------------------------------------------
// Set discretization information for auxiliary subfield.
void
pylith::feassemble::IntegratorPointwise::discretization(const char* name,
							const pylith::topology::FieldBase::DiscretizeInfo& feInfo)
{ // discretization
  _auxFieldsFEInfo[name] = feInfo;
} // discretization


// ----------------------------------------------------------------------
// Get discretization information for auxiliary subfield.
const pylith::topology::FieldBase::DiscretizeInfo& 
pylith::feassemble::IntegratorPointwise::discretization(const char* name) const
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


// ----------------------------------------------------------------------
// Verify configuration is acceptable.
void
pylith::feassemble::IntegratorPointwise::verifyConfiguration(const pylith::topology::Mesh& mesh) const
{ // verifyConfiguration
  PYLITH_METHOD_BEGIN;

  PYLITH_METHOD_END;
} // verifyConfiguration


// ----------------------------------------------------------------------
// Update state variables as needed.
void
pylith::feassemble::IntegratorPointwise::updateStateVars(const pylith::topology::Field& solution)
{ // updateState
  PYLITH_METHOD_BEGIN;


  PYLITH_METHOD_END;
} // updateStateVars



// End of file 
