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
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

#include <portinfo>

#include "Explicit.hh" // implementation of class methods

#include "pylith/topology/SolutionFields.hh" // USES SolutionFields
#include "pylith/topology/Stratum.hh" // USES Stratum
#include "pylith/topology/VisitorMesh.hh" // USES VecVisitorMesh

#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys

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
  PYLITH_METHOD_BEGIN;

  assert(_fields);

  // vel(t) = (disp(t+dt) - disp(t-dt)) / (2*dt)
  //        = (dispIncr(t+dt) + disp(t) - disp(t-dt)) / (2*dt)
  //
  // acc(t) = (disp(t+dt) - 2*disp(t) + disp(t-dt)) / (dt*dt)
  //        = (dispIncr(t+dt) - disp(t) + disp(t-dt)) / (dt*dt)

  const PylithScalar dt = _dt;
  const PylithScalar dt2 = dt*dt;
  const PylithScalar twodt = 2.0*dt;

  topology::Field& dispIncr = _fields->get("dispIncr(t->t+dt)");
  const spatialdata::geocoords::CoordSys* cs = dispIncr.mesh().coordsys();assert(cs);
  const int spaceDim = cs->spaceDim();

  // Get sections.
  topology::VecVisitorMesh dispIncrVisitor(dispIncr);
  PetscScalar* dispIncrArray = dispIncrVisitor.localArray();

  topology::Field& dispT = _fields->get("disp(t)");
  topology::VecVisitorMesh dispTVisitor(dispT);
  PetscScalar* dispTArray = dispTVisitor.localArray();

  topology::Field& dispTmdt = _fields->get("disp(t-dt)");
  topology::VecVisitorMesh dispTmdtVisitor(dispTmdt);
  PetscScalar* dispTmdtArray = dispTmdtVisitor.localArray();

  topology::Field& velocity = _fields->get("velocity(t)");
  topology::VecVisitorMesh velVisitor(velocity);
  PetscScalar* velArray = velVisitor.localArray();

  topology::Field& acceleration = _fields->get("acceleration(t)");
  topology::VecVisitorMesh accVisitor(acceleration);
  PetscScalar* accArray = accVisitor.localArray();

  // Get mesh vertices.
  PetscDM dmMesh = dispIncr.mesh().dmMesh();assert(dmMesh);
  topology::Stratum verticesStratum(dmMesh, topology::Stratum::DEPTH, 0);
  const PetscInt vStart = verticesStratum.begin();
  const PetscInt vEnd = verticesStratum.end();

  for(PetscInt v = vStart; v < vEnd; ++v) {
    const PetscInt dioff = dispIncrVisitor.sectionOffset(v);
    assert(spaceDim == dispIncrVisitor.sectionDof(v));

    const PetscInt dtoff = dispTVisitor.sectionOffset(v);
    assert(spaceDim == dispTVisitor.sectionDof(v));

    const PetscInt dmoff = dispTmdtVisitor.sectionOffset(v);
    assert(spaceDim == dispTmdtVisitor.sectionDof(v));

    const PetscInt voff = velVisitor.sectionOffset(v);
    assert(spaceDim == velVisitor.sectionDof(v));

    const PetscInt aoff = accVisitor.sectionOffset(v);
    assert(spaceDim == accVisitor.sectionDof(v));

    // TODO: I am not sure why these were updateAll() before, but if BCs need to be changed, then
    // the global update will probably need to be modified
    for (PetscInt i = 0; i < spaceDim; ++i) {
      velArray[voff+i] = (dispIncrArray[dioff+i] + dispTArray[dtoff+i] - dispTmdtArray[dmoff+i]) / twodt;
      accArray[aoff+i] = (dispIncrArray[dioff+i] - dispTArray[dtoff+i] + dispTmdtArray[dmoff+i]) / dt2;
    } // for
  } // for

  PetscLogFlops((vEnd - vStart) * 6*spaceDim);

  PYLITH_METHOD_END;
} // calcRateFields


// End of file
