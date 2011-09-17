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
// Copyright (c) 2010-2011 University of California, Davis
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
  deallocate();
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
  const ALE::Obj<RealSection>& dispIncrSection = dispIncr.section();
  assert(!dispIncrSection.isNull());
	 
  const ALE::Obj<RealSection>& dispTSection = _fields->get("disp(t)").section();
  assert(!dispTSection.isNull());

  const ALE::Obj<RealSection>& dispTmdtSection =
    _fields->get("disp(t-dt)").section();
  assert(!dispTmdtSection.isNull());

  scalar_array velVertex(spaceDim);
  const ALE::Obj<RealSection>& velSection = 
    _fields->get("velocity(t)").section();
  assert(!velSection.isNull());

  scalar_array accVertex(spaceDim);
  const ALE::Obj<RealSection>&  accSection = 
    _fields->get("acceleration(t)").section();
  assert(!accSection.isNull());

  // Get mesh vertices.
  const ALE::Obj<SieveMesh>& sieveMesh = dispIncr.mesh().sieveMesh();
  assert(!sieveMesh.isNull());
  const ALE::Obj<SieveMesh::label_sequence>& vertices = 
    sieveMesh->depthStratum(0);
  assert(!vertices.isNull());
  const SieveMesh::label_sequence::iterator verticesBegin = vertices->begin();
  const SieveMesh::label_sequence::iterator verticesEnd = vertices->end();
  
  for (SieveMesh::label_sequence::iterator v_iter=verticesBegin; 
       v_iter != verticesEnd;
       ++v_iter) {
    assert(spaceDim == dispIncrSection->getFiberDimension(*v_iter));
    const PylithScalar* dispIncrVertex = dispIncrSection->restrictPoint(*v_iter);

    assert(spaceDim == dispTSection->getFiberDimension(*v_iter));
    const PylithScalar* dispTVertex = dispTSection->restrictPoint(*v_iter);

    assert(spaceDim == dispTmdtSection->getFiberDimension(*v_iter));
    const PylithScalar* dispTmdtVertex = dispTmdtSection->restrictPoint(*v_iter);

    for (int i=0; i < spaceDim; ++i) {
      velVertex[i] = 
	(dispIncrVertex[i] + dispTVertex[i] - dispTmdtVertex[i]) / twodt;
      accVertex[i] = 
	(dispIncrVertex[i] - dispTVertex[i] + dispTmdtVertex[i]) / dt2;
    } // for
    
    assert(velSection->getFiberDimension(*v_iter) == spaceDim);
    velSection->updatePointAll(*v_iter, &velVertex[0]);

    assert(accSection->getFiberDimension(*v_iter) == spaceDim);
    accSection->updatePointAll(*v_iter, &accVertex[0]);
  } // for

  PetscLogFlops(vertices->size() * 6*spaceDim);
} // calcRateFields


// End of file
