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
// Compute velocity and acceleration at time t.
void
pylith::problems::Explicit::_calcRateFields(void)
{ // _calcRateFields
  assert(0 != _fields);

  // vel(t) = (disp(t+dt) - disp(t-dt)) / (2*dt)
  //        = (dispIncr(t+dt) + disp(t) - disp(t-dt)) / (2*dt)
  //
  // acc(t) = (disp(t+dt) - 2*disp(t) + disp(t-dt)) / (dt*dt)
  //        = (dispIncr(t+dt) - disp(t) + disp(t-dt)) / (dt*dt)

  const double dt = _dt;

  topology::Field<topology::Mesh>& dispIncr = _fields->get("dispIncr(t->t+dt)");
  const spatialdata::geocoords::CoordSys* cs = dispIncr.mesh().coordsys();
  assert(0 != cs);
  const int spaceDim = cs->spaceDim();
  
  if (!_fields->hasField("velocity(t)")) {
    _fields->add("velocity(t)", "velocity");
    topology::Field<topology::Mesh>& velocity = _fields->get("velocity(t)");
    velocity.cloneSection(dispIncr);

    _fields->add("acceleration(t)", "acceleration");
    topology::Field<topology::Mesh>& acceleration = 
      _fields->get("acceleration(t)");
    acceleration.cloneSection(dispIncr);
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

  double_array velVertex(spaceDim);
  const ALE::Obj<RealSection>& velSection = 
    _fields->get("velocity(t)").section();
  assert(!velSection.isNull());

  double_array accVertex(spaceDim);
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
    dispIncrSection->restrictPoint(*v_iter, &dispIncrVertex[0],
				   dispIncrVertex.size());
    dispTSection->restrictPoint(*v_iter, &dispTVertex[0],
				dispTVertex.size());
    dispTmdtSection->restrictPoint(*v_iter, &dispTmdtVertex[0],
				   dispTmdtVertex.size());

    velVertex = (dispIncrVertex + dispTVertex - dispTmdtVertex) / (2.0 * dt);
    accVertex = (dispIncrVertex - dispTVertex + dispTmdtVertex) / (dt * dt);
    
    assert(velSection->getFiberDimension(*v_iter) == spaceDim);
    velSection->updatePoint(*v_iter, &velVertex[0]);

    assert(accSection->getFiberDimension(*v_iter) == spaceDim);
    accSection->updatePoint(*v_iter, &accVertex[0]);
  } // for

  PetscLogFlops(vertices->size() * 6*spaceDim);
} // _calcRateFields



// End of file
