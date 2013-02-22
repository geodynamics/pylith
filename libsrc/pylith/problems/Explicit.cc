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

#include "Explicit.hh" // implementation of class methods

#include "pylith/topology/SolutionFields.hh" // USES SolutionFields

// ----------------------------------------------------------------------
typedef pylith::topology::Mesh::SieveMesh SieveMesh;
typedef pylith::topology::Mesh::RealSection RealSection;

// ----------------------------------------------------------------------
// Constructor
pylith::problems::Explicit::Explicit(void)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::problems::Explicit::~Explicit(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Compute velocity and acceleration at time t.
void
pylith::problems::Explicit::calcRateFields(void)
{ // calcRateFields
  assert(0 != _fields);

  // vel(t) = (disp(t+dt) - disp(t-dt)) / (2*dt)
  //        = (dispIncr(t+dt) + disp(t) - disp(t-dt)) / (2*dt)
  //
  // acc(t) = (disp(t+dt) - 2*disp(t) + disp(t-dt)) / (dt*dt)
  //        = (dispIncr(t+dt) - disp(t) + disp(t-dt)) / (dt*dt)

  const PylithScalar dt = _dt;
  const PylithScalar dt2 = dt*dt;
  const PylithScalar twodt = 2.0*dt;

  topology::Field<topology::Mesh>& dispIncr = _fields->get("dispIncr(t->t+dt)");
  const spatialdata::geocoords::CoordSys* cs = dispIncr.mesh().coordsys();
  assert(0 != cs);
  const int spaceDim = cs->spaceDim();
  
  // Get sections.
  PetscSection dispIncrSection = dispIncr.petscSection();
  Vec          dispIncrVec     = dispIncr.localVector();
  PetscScalar *dispIncrArray;
  assert(dispIncrSection);assert(dispIncrVec);
	 
  PetscSection dispTSection = _fields->get("disp(t)").petscSection();
  Vec          dispTVec     = _fields->get("disp(t)").localVector();
  PetscScalar *dispTArray;
  assert(dispTSection);assert(dispTVec);

  PetscSection dispTmdtSection = _fields->get("disp(t-dt)").petscSection();
  Vec          dispTmdtVec     = _fields->get("disp(t-dt)").localVector();
  PetscScalar *dispTmdtArray;
  assert(dispTmdtSection);assert(dispTmdtVec);

  PetscSection velSection = _fields->get("velocity(t)").petscSection();
  Vec          velVec     = _fields->get("velocity(t)").localVector();
  PetscScalar *velArray;
  assert(velSection);assert(velVec);

  PetscSection accSection = _fields->get("acceleration(t)").petscSection();
  Vec          accVec     = _fields->get("acceleration(t)").localVector();
  PetscScalar *accArray;
  assert(accSection);assert(accVec);

  // Get mesh vertices.
  DM             dmMesh = dispIncr.mesh().dmMesh();
  PetscInt       vStart, vEnd;
  PetscErrorCode err;

  assert(dmMesh);
  err = DMPlexGetDepthStratum(dmMesh, 0, &vStart, &vEnd);CHECK_PETSC_ERROR(err);

  err = VecGetArray(dispIncrVec, &dispIncrArray);CHECK_PETSC_ERROR(err);
  err = VecGetArray(dispTVec, &dispTArray);CHECK_PETSC_ERROR(err);
  err = VecGetArray(dispTmdtVec, &dispTmdtArray);CHECK_PETSC_ERROR(err);
  err = VecGetArray(velVec, &velArray);CHECK_PETSC_ERROR(err);
  err = VecGetArray(accVec, &accArray);CHECK_PETSC_ERROR(err);
  for(PetscInt v = vStart; v < vEnd; ++v) {
    PetscInt didof, dioff;

    err = PetscSectionGetDof(dispIncrSection, v, &didof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(dispIncrSection, v, &dioff);CHECK_PETSC_ERROR(err);
    assert(spaceDim == didof);
    PetscInt dtdof, dtoff;

    err = PetscSectionGetDof(dispTSection, v, &dtdof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(dispTSection, v, &dtoff);CHECK_PETSC_ERROR(err);
    assert(spaceDim == dtdof);
    PetscInt dmdof, dmoff;

    err = PetscSectionGetDof(dispTmdtSection, v, &dmdof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(dispTmdtSection, v, &dmoff);CHECK_PETSC_ERROR(err);
    assert(spaceDim == dmdof);
    PetscInt vdof, voff;

    err = PetscSectionGetDof(velSection, v, &vdof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(velSection, v, &voff);CHECK_PETSC_ERROR(err);
    assert(spaceDim == vdof);
    PetscInt adof, aoff;

    err = PetscSectionGetDof(accSection, v, &adof);CHECK_PETSC_ERROR(err);
    err = PetscSectionGetOffset(accSection, v, &aoff);CHECK_PETSC_ERROR(err);
    assert(spaceDim == adof);

    // TODO: I am not sure why these were updateAll() before, but if BCs need to be changed, then
    // the global update will probably need to be modified
    for (PetscInt i = 0; i < spaceDim; ++i) {
      velArray[voff+i] = (dispIncrArray[dioff+i] + dispTArray[dtoff+i] - dispTmdtArray[dmoff+i]) / twodt;
      accArray[aoff+i] = (dispIncrArray[dioff+i] - dispTArray[dtoff+i] + dispTmdtArray[dmoff+i]) / dt2;
    } // for
  } // for
  err = VecRestoreArray(dispIncrVec, &dispIncrArray);CHECK_PETSC_ERROR(err);
  err = VecRestoreArray(dispTVec, &dispTArray);CHECK_PETSC_ERROR(err);
  err = VecRestoreArray(dispTmdtVec, &dispTmdtArray);CHECK_PETSC_ERROR(err);
  err = VecRestoreArray(velVec, &velArray);CHECK_PETSC_ERROR(err);
  err = VecRestoreArray(accVec, &accArray);CHECK_PETSC_ERROR(err);

  PetscLogFlops((vEnd - vStart) * 6*spaceDim);
} // calcRateFields


// End of file
