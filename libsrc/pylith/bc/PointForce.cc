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
// Copyright (c) 2010-2012 University of California, Davis
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
  TimeDependentPoints::deallocate();
  feassemble::Integrator<feassemble::Quadrature<topology::Mesh> >::deallocate();
} // deallocate
  
// ----------------------------------------------------------------------
// Initialize boundary condition.
void
pylith::bc::PointForce::initialize(const topology::Mesh& mesh,
				    const PylithScalar upDir[3])
{ // initialize
  if (0 == _bcDOF.size())
    return;

  _getPoints(mesh);

  assert(_normalizer);
  const PylithScalar lengthScale = _normalizer->lengthScale();
  const PylithScalar pressureScale = _normalizer->pressureScale();
  const PylithScalar forceScale = pressureScale * lengthScale * lengthScale;

  _queryDatabases(mesh, forceScale, "force");
} // initialize

// ----------------------------------------------------------------------
// Integrate contributions to residual term (r) for operator.
void
pylith::bc::PointForce::integrateResidual(const topology::Field<topology::Mesh>& residual,
					  const PylithScalar t,
					  topology::SolutionFields* const fields)
{ // integrateResidualAssembled
  assert(_parameters);
  assert(_normalizer);

  // Calculate spatial and temporal variation of value for BC.
  _calculateValue(t);

  const int numPoints = _points.size();
  const int numBCDOF = _bcDOF.size();

  const topology::Mesh& mesh = residual.mesh();
  const spatialdata::geocoords::CoordSys* cs = mesh.coordsys();
  assert(cs);
  const int spaceDim = cs->spaceDim();

  PetscSection residualSection = residual.petscSection();assert(residualSection);
  PetscVec residualVec = residual.localVector();assert(residualVec);
  PetscScalar *residualArray;
  
  // Get global order
  PetscDM dmMesh = residual.dmMesh();assert(dmMesh);
  PetscSection globalSection = NULL;
  PetscErrorCode err;
  err = DMGetDefaultGlobalSection(dmMesh, &globalSection);CHECK_PETSC_ERROR(err);

  PetscSection valueSection = _parameters->get("value").petscSection();assert(valueSection);
  PetscVec valueVec = _parameters->get("value").localVector();assert(valueVec);
  PetscScalar* valueArray = NULL;
  err = VecGetArray(valueVec, &valueArray);CHECK_PETSC_ERROR(err);
  err = VecGetArray(residualVec, &residualArray);CHECK_PETSC_ERROR(err);

  for (int iPoint=0; iPoint < numPoints; ++iPoint) {
    const PetscInt p_bc = _points[iPoint]; // Get point label.
    PetscInt       goff;

    // Contribute to residual only if point is local.
    err = PetscSectionGetOffset(globalSection, p_bc, &goff);CHECK_PETSC_ERROR(err);
    if (goff < 0) continue;

    PetscInt vdof, voff, rdof, roff;

    err = PetscSectionGetDof(valueSection, p_bc, &vdof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(valueSection, p_bc, &voff);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetDof(residualSection, p_bc, &rdof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(residualSection, p_bc, &roff);CHECK_PETSC_ERROR(err);
    assert(vdof == numBCDOF);
    assert(rdof == spaceDim);

    for (int iDOF=0; iDOF < numBCDOF; ++iDOF)
      residualArray[roff+_bcDOF[iDOF]] += valueArray[voff+iDOF];
  } // for
  err = VecRestoreArray(valueVec, &valueArray);CHECK_PETSC_ERROR(err);
  err = VecRestoreArray(residualVec, &residualArray);CHECK_PETSC_ERROR(err);
} // integrateResidualAssembled

// ----------------------------------------------------------------------
// Verify configuration is acceptable.
void
pylith::bc::PointForce::verifyConfiguration(const topology::Mesh& mesh) const
{ // verifyConfiguration
  BoundaryCondition::verifyConfiguration(mesh);
  TimeDependent::verifyConfiguration(mesh);
} // verifyConfiguration


// End of file 
