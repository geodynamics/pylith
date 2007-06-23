// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
//                           Brad T. Aagaard
//                        U.S. Geological Survey
//
// {LicenseText}
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "FaultCohesive.hh" // implementation of object methods

#include "CohesiveTopology.hh" // USES CohesiveTopology::create()

#include "pylith/utils/sievetypes.hh" // USES PETSc Mesh
#include "pylith/utils/array.hh" // USES double_array

#include <assert.h> // USES assert()

// ----------------------------------------------------------------------
// Default constructor.
pylith::faults::FaultCohesive::FaultCohesive(void) :
  _faultMesh(new ALE::Obj<ALE::Mesh>)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::faults::FaultCohesive::~FaultCohesive(void)
{ // destructor
  delete _faultMesh; _faultMesh = 0;
} // destructor

// ----------------------------------------------------------------------
// Adjust mesh topology for fault implementation.
void
pylith::faults::FaultCohesive::adjustTopology(const ALE::Obj<ALE::Mesh>& mesh)
{ // adjustTopology
  assert(std::string("") != label());

  // Get group of vertices associated with fault
  const ALE::Obj<int_section_type>& groupField = 
    mesh->getIntSection(label());
  assert(!groupField.isNull());

  CohesiveTopology::create(_faultMesh, mesh, groupField, id(),
                           _useLagrangeConstraints());
} // adjustTopology

// ----------------------------------------------------------------------
// Compute weighted orientation of fault for cohesive cell between
// 1-D elements.
void
pylith::faults::FaultCohesive::_orient1D(double_array* orientation,
					 const double_array& jacobian,
					 const double jacobianDet,
					 const double_array& upDir)
{ // _orient1D
  assert(0 != orientation);
  assert(1 == orientation->size());
  (*orientation) = 1.0;
} // _orient1D
		
// ----------------------------------------------------------------------
// Compute weighted orientation of fault for cohesive cell between
// 2-D elements.
void
pylith::faults::FaultCohesive::_orient2D(double_array* orientation,
					 const double_array& jacobian,
					 const double jacobianDet,
					 const double_array& upDir)
{ // _orient2D
  const int orientSize = 4;
  assert(0 != orientation);
  assert(orientSize == orientation->size());
  const int jacobianSize = 2;
  assert(jacobianSize == jacobian.size());

  // cellDim is 1
  const int spaceDim = 2;

  const double j1 = jacobian[0];
  const double j2 = jacobian[1];
  (*orientation)[0] =  j1;
  (*orientation)[1] =  j2;
  (*orientation)[2] =  j2;
  (*orientation)[3] = -j1;
} // _orient2D
		
// ----------------------------------------------------------------------
// Compute weighted orientation of fault for cohesive cell between
// 3-D elements.
void
pylith::faults::FaultCohesive::_orient3D(double_array* orientation,
					 const double_array& jacobian,
					 const double jacobianDet,
					 const double_array& upDir)
{ // _orient3D
  const int orientSize = 9;
  assert(0 != orientation);
  assert(orientSize == orientation->size());
  const int jacobianSize = 6;
  assert(jacobianSize == jacobian.size());
  assert(3 == upDir.size());

  const int cellDim = 2;
  const int spaceDim = 3;

  const double j00 = jacobian[0];
  const double j01 = jacobian[1];
  const double j10 = jacobian[2];
  const double j11 = jacobian[3];
  const double j20 = jacobian[4];
  const double j21 = jacobian[5];

  // Compute normal using Jacobian
  double r0 =  j10*j21 - j20*j11;
  double r1 = -j00*j21 + j20*j01;
  double r2 =  j00*j11 - j10*j01;
  // Make unit vector
  double mag = sqrt(r0*r0 + r1*r1 + r2*r2);
  assert(mag > 0.0);
  r0 /= mag;
  r1 /= mag;
  r2 /= mag;
  
  // Compute along-strike direction; cross product of "up" and normal
  double p0 =  upDir[1]*r2 - upDir[2]*r1;
  double p1 = -upDir[0]*r2 + upDir[2]*r0;
  double p2 =  upDir[0]*r1 - upDir[1]*r0;
  // Make unit vector
  mag = sqrt(p0*p0 + p1*p1 + p2*p2);
  assert(mag > 0.0);
  p0 /= mag;
  p1 /= mag;
  p2 /= mag;
  
  // Compute up-dip direction; cross product of normal and along-strike
  const double q0 =  r1*p2 - r2*p1;
  const double q1 = -r0*p2 + r2*p0;
  const double q2 =  r0*p1 - r1*p0;
  mag = sqrt(q0*q0 + q1*q1 + q2*q2);
  assert(mag > 0.0);
  
  const double wt = jacobianDet;
  (*orientation)[0] =  p0*wt;
  (*orientation)[1] =  p1*wt;
  (*orientation)[2] =  p2*wt;
  (*orientation)[3] =  q0*wt;
  (*orientation)[4] =  q1*wt;
  (*orientation)[5] =  q2*wt;
  (*orientation)[6] =  r0*wt;
  (*orientation)[7] =  r1*wt;
  (*orientation)[8] =  r2*wt;
} // _orient3D


// End of file 
