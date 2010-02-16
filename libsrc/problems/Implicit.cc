// -*- C++ -*-
//
// ======================================================================
//
//                           Brad T. Aagaard
//                        U.S. Geological Survey
//
// {LicenseText}
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
pylith::problems::Implicit::_calcVelocity(void)
{ // _calcVelocity
  assert(0 != _fields);

  // vel(t) = (disp(t+dt) - disp(t)) / dt
  //        = dispIncr(t+dt) / dt
  const double dt = _dt;

  topology::Field<topology::Mesh>& dispIncr = _fields->get("dispIncr(t->t+dt)");
  const spatialdata::geocoords::CoordSys* cs = dispIncr.mesh().coordsys();
  assert(0 != cs);
  const int spaceDim = cs->spaceDim();
  
  if (!_fields->hasField("velocity(t)")) {
    _fields->add("velocity(t)", "velocity");
    topology::Field<topology::Mesh>& velocity = _fields->get("velocity(t)");
    velocity.cloneSection(dispIncr);
  } // if

  // Get sections.
  double_array dispIncrVertex(spaceDim);
  const ALE::Obj<RealSection>& dispIncrSection = dispIncr.section();
  assert(!dispIncrSection.isNull());
	 
  double_array velocityVertex(spaceDim);
  topology::Field<topology::Mesh>& velocity = _fields->get("velocity(t)");
  const ALE::Obj<RealSection>& velocitySection = velocity.section();
  assert(!velocitySection.isNull());

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
    velocityVertex = dispIncrVertex / dt;
    
    assert(velocitySection->getFiberDimension(*v_iter) == spaceDim);
    velocitySection->updatePoint(*v_iter, &velocityVertex[0]);
  } // for
  PetscLogFlops(vertices->size() * spaceDim);

} // _calcVelocity


// End of file
