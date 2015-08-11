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
#include "pylith/topology/Fields.hh" // USES Fields
#include "pylith/topology/Jacobian.hh" // USES Jacobian
#include "pylith/topology/SolutionFields.hh" // USES SolutionFields
#include "pylith/topology/Stratum.hh" // USES Stratum
#include "pylith/topology/VisitorMesh.hh" // USES VecVisitorMesh
#include "pylith/topology/CoordsVisitor.hh" // USES CoordsVisitor

#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional
#include "spatialdata/spatialdb/GravityField.hh" // USES GravityField

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

  Integrator::deallocate();
  delete _auxFields; _auxFields = 0;

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
// Update state variables as needed.
void
pylith::feassemble::IntegratorPointwise::updateStateVars(const PylithScalar t,
							  topology::SolutionFields* const fields)
{ // updateState
  PYLITH_METHOD_BEGIN;

  assert(_quadrature);
  assert(fields);

#if 0
  // No need to update state vars if material doesn't have any.
  if (!_material->hasStateVars())
    PYLITH_METHOD_END;

  // Get cell information that doesn't depend on particular cell
  const int cellDim = _quadrature->cellDim();
  const int numQuadPts = _quadrature->numQuadPts();
  const int numBasis = _quadrature->numBasis();
  const int spaceDim = _quadrature->spaceDim();
  const int numCorners = _quadrature->refGeometry().numCorners();
  const int tensorSize = _material->tensorSize();
  totalStrain_fn_type calcTotalStrainFn;
  if (2 == cellDim) {
    calcTotalStrainFn = &pylith::feassemble::IntegratorPointwise::_calcTotalStrain2D;
  } else if (3 == cellDim) {
    calcTotalStrainFn = &pylith::feassemble::IntegratorPointwise::_calcTotalStrain3D;
  } else {
      std::cerr << "Bad cell dimension '" << cellDim << "'." << std::endl;
      assert(0);
      throw std::logic_error("Bad cell dimension in IntegratorPointwise::updateStateVars().");
  } // else

  // Allocate arrays for cell data.
  scalar_array strainCell(numQuadPts*tensorSize);
  strainCell = 0.0;

  // Get cell information
  PetscDM dmMesh = fields->mesh().dmMesh();assert(dmMesh);
  assert(_materialIS);
  const PetscInt* cells = _materialIS->points();
  const PetscInt numCells = _materialIS->size();

  // Setup visitors.
  scalar_array dispCell(numBasis*spaceDim);
  topology::VecVisitorMesh dispVisitor(fields->get("disp(t)"), "displacement");
  dispVisitor.optimizeClosure();

  scalar_array coordsCell(numCorners*spaceDim);
  topology::CoordsVisitor coordsVisitor(dmMesh);

  // Loop over cells
  for(PetscInt c = 0; c < numCells; ++c) {
    const PetscInt cell = cells[c];
    // Retrieve geometry information for current cell
    coordsVisitor.getClosure(&coordsCell, cell);
    _quadrature->computeGeometry(&coordsCell[0], coordsCell.size(), cell);
    const scalar_array& basisDeriv = _quadrature->basisDeriv();

    // Get physical properties and state variables for cell.
    _material->retrievePropsAndVars(cell);

    // Restrict input fields to cell
    dispVisitor.getClosure(&dispCell, cell);

    // Compute strains
    calcTotalStrainFn(&strainCell, basisDeriv, &dispCell[0], numBasis, spaceDim, numQuadPts);

    // Update material state
    _material->updateStateVars(strainCell, cell);
  } // for
#endif

  PYLITH_METHOD_END;
} // updateStateVars

// ----------------------------------------------------------------------
// Verify configuration is acceptable.
void
pylith::feassemble::IntegratorPointwise::verifyConfiguration(const topology::Mesh& mesh) const
{ // verifyConfiguration
  PYLITH_METHOD_BEGIN;

  PYLITH_METHOD_END;
} // verifyConfiguration

// End of file 
