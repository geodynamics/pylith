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

#include "Implicit.hh" // implementation of class methods

#include "pylith/topology/SolutionFields.hh" // USES SolutionFields

// ----------------------------------------------------------------------
typedef pylith::topology::Mesh::SieveMesh SieveMesh;
typedef pylith::topology::Mesh::RealSection RealSection;

// ----------------------------------------------------------------------
// Constructor
pylith::problems::Implicit::Implicit(void)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::problems::Implicit::~Implicit(void)
{ // destructor
  deallocate();
} // destructor

// ----------------------------------------------------------------------
// Compute velocity at time t.
void
pylith::problems::Implicit::calcRateFields(void)
{ // calcRateFields
  assert(0 != _fields);

  // vel(t) = (disp(t+dt) - disp(t)) / dt
  //        = dispIncr(t+dt) / dt
  const PylithScalar dt = _dt;

  topology::Field<topology::Mesh>& dispIncr = _fields->get("dispIncr(t->t+dt)");
  const spatialdata::geocoords::CoordSys* cs = dispIncr.mesh().coordsys();
  assert(0 != cs);
  const int spaceDim = cs->spaceDim();
  
  // Get sections.
  PetscSection dispIncrSection = dispIncr.petscSection();
  Vec          dispIncrVec     = dispIncr.localVector();
  PetscScalar *dispIncrArray;
  assert(dispIncrSection);assert(dispIncrVec);
	 
  PetscSection velSection = _fields->get("velocity(t)").petscSection();
  Vec          velVec     = _fields->get("velocity(t)").localVector();
  PetscScalar *velArray;
  assert(velSection);assert(velVec);

  // Get mesh vertices.
  DM             dmMesh = dispIncr.dmMesh();
  PetscInt       vStart, vEnd;
  PetscErrorCode err;

  assert(dmMesh);
  err = DMComplexGetDepthStratum(dmMesh, 0, &vStart, &vEnd);CHECK_PETSC_ERROR(err);
  err = VecGetArray(dispIncrVec, &dispIncrArray);CHECK_PETSC_ERROR(err);
  err = VecGetArray(velVec, &velArray);CHECK_PETSC_ERROR(err);
  for(PetscInt v = vStart; v < vEnd; ++v) {
    PetscInt dof, off, vdof, voff;

    err = PetscSectionGetDof(dispIncrSection, v, &dof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(dispIncrSection, v, &off);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetDof(velSection, v, &vdof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(velSection, v, &voff);CHECK_PETSC_ERROR(err);
    assert(dof == spaceDim);assert(vdof == spaceDim);
    for(PetscInt d = 0; d < dof; ++d) {
      velArray[voff+d] = dispIncrArray[off+d] / dt;
    }
  } // for
  err = VecRestoreArray(dispIncrVec, &dispIncrArray);CHECK_PETSC_ERROR(err);
  err = VecRestoreArray(velVec, &velArray);CHECK_PETSC_ERROR(err);
  PetscLogFlops((vEnd - vStart) * spaceDim);
} // calcRateFields


// End of file
