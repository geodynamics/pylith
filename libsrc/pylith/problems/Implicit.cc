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

#include "Implicit.hh" // implementation of class methods

#include "pylith/topology/SolutionFields.hh" // USES SolutionFields
#include "pylith/topology/Stratum.hh" // USES Stratum
#include "pylith/topology/VisitorMesh.hh" // USES VecVisitorMesh

#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys

// ----------------------------------------------------------------------
// Constructor
pylith::problems::Implicit::Implicit(void)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::problems::Implicit::~Implicit(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Compute velocity at time t.
void
pylith::problems::Implicit::calcRateFields(void)
{ // calcRateFields
  PYLITH_METHOD_BEGIN;

  assert(_fields);

  // vel(t) = (disp(t+dt) - disp(t)) / dt
  //        = dispIncr(t+dt) / dt
  const PylithScalar dt = _dt;

  topology::Field& dispIncr = _fields->get("dispIncr(t->t+dt)");
  const spatialdata::geocoords::CoordSys* cs = dispIncr.mesh().coordsys();assert(cs);
  const int spaceDim = cs->spaceDim();

  // Get sections.
  topology::VecVisitorMesh dispIncrVisitor(dispIncr);
  PetscScalar* dispIncrArray = dispIncrVisitor.localArray();

  topology::Field& velocity = _fields->get("velocity(t)");
  topology::VecVisitorMesh velVisitor(velocity);
  PetscScalar* velArray = velVisitor.localArray();

  // Get mesh vertices.
  PetscDM dmMesh = dispIncr.mesh().dmMesh();assert(dmMesh);
  topology::Stratum verticesStratum(dmMesh, topology::Stratum::DEPTH, 0);
  const PetscInt vStart = verticesStratum.begin();
  const PetscInt vEnd = verticesStratum.end();

  for(PetscInt v = vStart; v < vEnd; ++v) {
    const PetscInt dioff = dispIncrVisitor.sectionOffset(v);
    assert(spaceDim == dispIncrVisitor.sectionDof(v));

    const PetscInt voff = velVisitor.sectionOffset(v);
    assert(spaceDim == velVisitor.sectionDof(v));

    for(PetscInt d = 0; d < spaceDim; ++d) {
      velArray[voff+d] = dispIncrArray[dioff+d] / dt;
    }
  } // for

  PetscLogFlops((vEnd - vStart) * spaceDim);

  PYLITH_METHOD_END;
} // calcRateFields


// End of file
