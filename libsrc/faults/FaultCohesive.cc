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

#include "pylith/feassemble/Quadrature.hh" // USES Quadrature

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
  assert("" != label());

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
					 const double_array& jacobianDet,
					 const double_array& upDir,
					 const int numLocs)
{ // _orient1D
  assert(0 != orientation);
  assert(numLocs == orientation->size());
  (*orientation) = 1.0;
} // _orient1D
		
// ----------------------------------------------------------------------
// Compute weighted orientation of fault for cohesive cell between
// 2-D elements.
void
pylith::faults::FaultCohesive::_orient2D(double_array* orientation,
					 const double_array& jacobian,
					 const double_array& jacobianDet,
					 const double_array& upDir,
					 const int numLocs)
{ // _orient2D
  const int orientSize = 4;
  assert(0 != orientation);
  assert(orientSize*numLocs == orientation->size());
  const int jacobianSize = 2;
  assert(numLocs*jacobianSize == jacobian.size());

  // cellDim is 1
  const int spaceDim = 2;
  for (int iLoc=0; iLoc < numLocs; ++iLoc) {
    const double j1 = jacobian[iLoc*spaceDim  ];
    const double j2 = jacobian[iLoc*spaceDim+1];
    (*orientation)[iLoc*orientSize  ] =  j1;
    (*orientation)[iLoc*orientSize+1] =  j2;
    (*orientation)[iLoc*orientSize+2] =  j2;
    (*orientation)[iLoc*orientSize+3] = -j1;
  } // for
} // _orient2D
		
// ----------------------------------------------------------------------
// Compute weighted orientation of fault for cohesive cell between
// 3-D elements.
void
pylith::faults::FaultCohesive::_orient3D(double_array* orientation,
					 const double_array& jacobian,
					 const double_array& jacobianDet,
					 const double_array& upDir,
					 const int numLocs)
{ // _orient3D
  const int orientSize = 9;
  assert(0 != orientation);
  assert(orientSize*numLocs == orientation->size());
  const int jacobianSize = 6;
  assert(numLocs*jacobianSize == jacobian.size());
  assert(numLocs == jacobianDet.size());
  assert(3 == upDir.size());

  const int cellDim = 2;
  const int spaceDim = 3;
  for (int iLoc=0; iLoc < numLocs; ++iLoc) {
    const double j00 = jacobian[iLoc*spaceDim*cellDim  ];
    const double j10 = jacobian[iLoc*spaceDim*cellDim+1];
    const double j20 = jacobian[iLoc*spaceDim*cellDim+2];
    const double j01 = jacobian[iLoc*spaceDim*cellDim+3];
    const double j11 = jacobian[iLoc*spaceDim*cellDim+4];
    const double j21 = jacobian[iLoc*spaceDim*cellDim+5];

    // Compute normal using Jacobian
    double r0 =  j10*j21 - j20*j11;
    double r1 = -j00*j21 + j20*j01;
    double r2 =  j00*j11 - j10*j01;
    // Make unit vector
    double mag = sqrt(r0*r0 + r1*r1 + r2*r2);
    r0 /= mag;
    r1 /= mag;
    r2 /= mag;

    // Compute along-strike direction; cross product of "up" and normal
    double p0 =  upDir[1]*r2 - upDir[2]*r1;
    double p1 = -upDir[0]*r2 + upDir[2]*r0;
    double p2 =  upDir[0]*r1 - upDir[1]*r0;
    // Make unit vector
    mag = sqrt(p0*p0 + p1*p1 + p2*p2);
    p0 /= mag;
    p1 /= mag;
    p2 /= mag;

    // Compute up-dip direction; cross product of normal and along-strike
    const double q0 =  r1*p2 - r2*p1;
    const double q1 = -r0*p2 + r2*p0;
    const double q2 =  r0*p1 - r1*p0;

    const double wt = jacobianDet[iLoc];
    (*orientation)[iLoc*orientSize  ] =  p0*wt;
    (*orientation)[iLoc*orientSize+1] =  q0*wt;
    (*orientation)[iLoc*orientSize+2] =  r0*wt;
    (*orientation)[iLoc*orientSize+3] =  p1*wt;
    (*orientation)[iLoc*orientSize+4] =  q1*wt;
    (*orientation)[iLoc*orientSize+5] =  r1*wt;
    (*orientation)[iLoc*orientSize+6] =  p2*wt;
    (*orientation)[iLoc*orientSize+7] =  q2*wt;
    (*orientation)[iLoc*orientSize+8] =  r2*wt;
  } // for
} // _orient3D


// End of file 
