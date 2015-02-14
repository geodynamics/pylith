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

#include "PointForce.hh" // implementation of object methods

#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/Fields.hh" // USES Fields
#include "pylith/topology/SolutionFields.hh" // USES SolutionFields
#include "pylith/topology/VisitorMesh.hh" // USES VecVisitorMesh
#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include <cstring> // USES strcpy()
#include <cassert> // USES assert()
#include <stdexcept> // USES std::runtime_error
#include <sstream> // USES std::ostringstream

// ----------------------------------------------------------------------
// Default constructor.
pylith::bc::PointForce::PointForce(void)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::bc::PointForce::~PointForce(void)
{ // destructor
  deallocate();
} // destructor

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::bc::PointForce::deallocate(void)
{ // deallocate
  PYLITH_METHOD_BEGIN;

  TimeDependentPoints::deallocate();
  feassemble::Integrator::deallocate();

  PYLITH_METHOD_END;
} // deallocate
  
// ----------------------------------------------------------------------
// Initialize boundary condition.
void
pylith::bc::PointForce::initialize(const topology::Mesh& mesh,
				    const PylithScalar upDir[3])
{ // initialize
  PYLITH_METHOD_BEGIN;

  if (0 == _bcDOF.size())
    PYLITH_METHOD_END;

  _getPoints(mesh);

  assert(_normalizer);
  const PylithScalar lengthScale = _normalizer->lengthScale();
  const PylithScalar pressureScale = _normalizer->pressureScale();
  const PylithScalar forceScale = pressureScale * lengthScale * lengthScale;

  _queryDatabases(mesh, forceScale, "force");

  PYLITH_METHOD_END;
} // initialize

// ----------------------------------------------------------------------
// Integrate contributions to residual term (r) for operator.
void
pylith::bc::PointForce::integrateResidual(const topology::Field& residual,
					  const PylithScalar t,
					  topology::SolutionFields* const fields)
{ // integrateResidualAssembled
  PYLITH_METHOD_BEGIN;

  assert(_parameters);
  assert(_normalizer);

  // Calculate spatial and temporal variation of value for BC.
  _calculateValue(t);

  const topology::Mesh& mesh = residual.mesh();
  const spatialdata::geocoords::CoordSys* cs = mesh.coordsys();
  assert(cs);
  const int spaceDim = cs->spaceDim();

  topology::VecVisitorMesh residualVisitor(residual);
  PetscScalar* residualArray = residualVisitor.localArray();

  topology::Field& valueField = _parameters->get("value");
  topology::VecVisitorMesh valueVisitor(valueField);
  PetscScalar* valueArray = valueVisitor.localArray();

  // Get global order
  PetscSection globalSection = residual.globalSection();assert(globalSection);
  PetscErrorCode err = 0;

  const int numPoints = _points.size();
  const int numBCDOF = _bcDOF.size();
  for (int iPoint=0; iPoint < numPoints; ++iPoint) {
    const PetscInt p_bc = _points[iPoint]; // Get point label.

    // Contribute to residual only if point is local.
    PetscInt goff = -1;
    err = PetscSectionGetOffset(globalSection, p_bc, &goff);PYLITH_CHECK_ERROR(err);
    if (goff < 0) continue;

    const PetscInt roff = residualVisitor.sectionOffset(p_bc);
    const PetscInt voff = valueVisitor.sectionOffset(p_bc);
    assert(numBCDOF == valueVisitor.sectionDof(p_bc));
    assert(spaceDim == residualVisitor.sectionDof(p_bc));

    for (int iDOF=0; iDOF < numBCDOF; ++iDOF) {
      residualArray[roff+_bcDOF[iDOF]] += valueArray[voff+iDOF];
    } // for
  } // for

  PYLITH_METHOD_END;
} // integrateResidualAssembled

// ----------------------------------------------------------------------
// Verify configuration is acceptable.
void
pylith::bc::PointForce::verifyConfiguration(const topology::Mesh& mesh) const
{ // verifyConfiguration
  PYLITH_METHOD_BEGIN;

  BoundaryCondition::verifyConfiguration(mesh);
  TimeDependent::verifyConfiguration(mesh);

  PYLITH_METHOD_END;
} // verifyConfiguration


// End of file 
