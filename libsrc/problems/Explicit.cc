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
// Compute velocity at time t.
void
pylith::problems::Explicit::_calcVelocity(void)
{ // _calcVelocity
  assert(0 != _fields);

  // vel(t) = (disp(t+dt) - disp(t-dt)) / (2*dt)
  //        = (dispIncr(t+dt) + disp(t) - disp(t-dt)) / (2*dt)
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
	 
  double_array dispTVertex(spaceDim);
  const ALE::Obj<RealSection>& dispTSection = _fields->get("disp(t)").section();
  assert(!dispTSection.isNull());

  double_array dispTmdtVertex(spaceDim);
  const ALE::Obj<RealSection>& dispTmdtSection =
    _fields->get("disp(t-dt)").section();
  assert(!dispTmdtSection.isNull());

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
    dispTSection->restrictPoint(*v_iter, &dispTVertex[0],
				dispTVertex.size());
    dispTmdtSection->restrictPoint(*v_iter, &dispTmdtVertex[0],
				   dispTmdtVertex.size());
    velocityVertex = (dispIncrVertex + dispTVertex - dispTmdtVertex) / (2.0 * dt);
    
    assert(velocitySection->getFiberDimension(*v_iter) == spaceDim);
    velocitySection->updatePoint(*v_iter, &velocityVertex[0]);
  } // for
  PetscLogFlops(vertices->size() * spaceDim);

} // _calcVelocity


// End of file
