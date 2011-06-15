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
  const double dt = _dt;

  topology::Field<topology::Mesh>& dispIncr = _fields->get("dispIncr(t->t+dt)");
  const spatialdata::geocoords::CoordSys* cs = dispIncr.mesh().coordsys();
  assert(0 != cs);
  const int spaceDim = cs->spaceDim();
  
  // Get sections.
  double_array dispIncrVertex(spaceDim);
  const ALE::Obj<RealSection>& dispIncrSection = dispIncr.section();
  assert(!dispIncrSection.isNull());
	 
  double_array velVertex(spaceDim);
  const ALE::Obj<RealSection>& velSection = 
    _fields->get("velocity(t)").section();
  assert(!velSection.isNull());

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
    dispIncrSection->restrictPoint(*v_iter, &dispIncrVertex[0],
				   dispIncrVertex.size());
    velVertex = dispIncrVertex / dt;
    
    assert(velSection->getFiberDimension(*v_iter) == spaceDim);
    velSection->updatePointAll(*v_iter, &velVertex[0]);
  } // for
  PetscLogFlops(vertices->size() * spaceDim);

} // calcRateFields


// End of file
