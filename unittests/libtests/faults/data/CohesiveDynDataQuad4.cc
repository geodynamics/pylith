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


/* Original mesh
 *
 * Cells are 0-1, vertices are 2-7.
 *
 *       3 -------- 5 -------- 7
 *       |          |          |
 *       |          |          |
 *       |          |          |
 *       |          |          |
 *       |          |          |
 *       |          |          |
 *       |          |          |
 *       |          |          |
 *       2 -------- 4 -------- 6
 *
 * After adding cohesive elements
 *
 * Cells are 0-1,2 vertices are 3-10,11-12.
 *
 *       4 -------- 6 -12--10 -------- 8
 *       |          |       |          |
 *       |          |       |          |
 *       |          |       |          |
 *       |          |       |          |
 *       |          |       |          |
 *       |          |       |          |
 *       |          |       |          |
 *       |          |       |          |
 *       3 -------- 5 -11-- 9 -------- 7
 */

#include "CohesiveDynDataQuad4.hh"

const char* pylith::faults::CohesiveDynDataQuad4::_meshFilename =
  "data/quad4.mesh";

const int pylith::faults::CohesiveDynDataQuad4::_spaceDim = 2;

const int pylith::faults::CohesiveDynDataQuad4::_cellDim = 1;

const int pylith::faults::CohesiveDynDataQuad4::_numBasis = 2;

const int pylith::faults::CohesiveDynDataQuad4::_numQuadPts = 2;

const PylithScalar pylith::faults::CohesiveDynDataQuad4::_quadPts[] = {
  -1.0, 1.0,
};

const PylithScalar pylith::faults::CohesiveDynDataQuad4::_quadWts[] = {
  1.0, 1.0
};

const PylithScalar pylith::faults::CohesiveDynDataQuad4::_basis[] = {
  1.0, 0.0,
  0.0, 1.0,
};

const PylithScalar pylith::faults::CohesiveDynDataQuad4::_basisDeriv[] = {
  -0.5, 0.5,
  -0.5, 0.5,
};

const PylithScalar pylith::faults::CohesiveDynDataQuad4::_verticesRef[] = {
  -1.0, 1.0
};

const int pylith::faults::CohesiveDynDataQuad4::_id = 10;

const char* pylith::faults::CohesiveDynDataQuad4::_label = "fault";

const char* pylith::faults::CohesiveDynDataQuad4::_initialTractFilename = 
  "data/quad4_initialtract.spatialdb";

const PylithScalar pylith::faults::CohesiveDynDataQuad4::_fieldT[] = {
  8.1, 9.1,
  8.2, 9.2,
  8.3, 9.3, // 4
  8.4, 9.4, // 5
  8.5, 9.5,
  8.6, 9.6,
  8.3, 9.3, // 8
  8.4, 9.4, // 9
  8.8, 9.8, // 10
  8.0, 9.0, // 11
};

const PylithScalar pylith::faults::CohesiveDynDataQuad4::_jacobian[] = {
 0.1, 0.1, // 2x
 0.2, 0.3,
 0.4, 0.5,
 0.6, 0.7,
 0.8, 0.9,
 1.0, 1.1,
 1.2, 1.3,
 1.4, 1.5,
 1.6, 1.7,
 1.8, 1.9,
 2.0, 2.1, // 2y
 2.1, 2.2,
 2.3, 2.4,
 2.5, 2.6,
 2.7, 2.8,
 2.9, 3.0,
 3.1, 3.2,
 3.3, 3.4,
 3.5, 3.6,
 3.7, 3.8,
 3.9, 4.0, // 3x
 4.0, 4.1,
 4.2, 4.3,
 4.4, 4.5,
 4.6, 4.7,
 4.8, 4.9,
 5.0, 5.1,
 5.2, 5.3,
 5.4, 5.5,
 5.6, 5.7,
 5.8, 5.9, // 3y
 6.0, 6.1,
 6.1, 6.2,
 6.3, 6.4,
 6.5, 6.6,
 6.7, 6.8,
 6.9, 7.0,
 7.1, 7.2,
 7.3, 7.4,
 7.5, 7.6,
 7.7, 7.8, // 4x
 7.9, 8.0,
+4.0,-1.2, // 4
-2.2,-2.3, // 5
 8.4, 8.5,
 8.6, 8.7,
 8.8, 8.9,
 9.0, 9.1,
 9.2, 9.3,
 9.3, 9.4,
 3.7, 4.8, // 4y
 3.9, 4.0,
-1.2,+5.0, // 4
-1.3,-3.2, // 5
 4.4, 5.5,
 4.6, 5.7,
 4.8, 5.9,
 4.0, 5.1,
 4.2, 5.3,
 4.3, 5.4,
 7.7, 7.8, // 5x
 7.9, 8.0,
-2.2,-1.3, // 4
+4.1,-4.3, // 5
 8.4, 8.5,
 8.6, 8.7,
 8.8, 8.9,
 9.0, 9.1,
 9.2, 9.3,
 9.3, 9.4,
 3.7, 4.8, // 5y
 3.9, 4.0,
-2.3,-3.2, // 4
-4.3,+5.1, // 5
 4.4, 5.5,
 4.6, 5.7,
 4.8, 5.9,
 4.0, 5.1,
 4.2, 5.3,
 4.3, 5.4,
 0.1, 0.1, // 6x
 0.2, 0.3,
 0.4, 0.5,
 0.6, 0.7,
 0.8, 0.9,
 1.0, 1.1,
 1.2, 1.3,
 1.4, 1.5,
 1.6, 1.7,
 1.8, 1.9,
 2.0, 2.1, // 6y
 2.1, 2.2,
 2.3, 2.4,
 2.5, 2.6,
 2.7, 2.8,
 2.9, 3.0,
 3.1, 3.2,
 3.3, 3.4,
 3.5, 3.6,
 3.7, 3.8,
 3.9, 4.0, // 7x
 4.0, 4.1,
 4.2, 4.3,
 4.4, 4.5,
 4.6, 4.7,
 4.8, 4.9,
 5.0, 5.1,
 5.2, 5.3,
 5.4, 5.5,
 5.6, 5.7,
 5.8, 5.9, // 7y
 6.0, 6.1,
 6.1, 6.2,
 6.3, 6.4,
 6.5, 6.6,
 6.7, 6.8,
 6.9, 7.0,
 7.1, 7.2,
 7.3, 7.4,
 7.5, 7.6,
 7.7, 7.8, // 8x
 7.9, 8.0,
 8.4, 8.5,
 8.6, 8.7,
 8.8, 8.9,
 9.0, 9.1,
+5.0,-1.2, // 8
-2.2,-2.3, // 9
 9.2, 9.3,
 9.3, 9.4,
 3.7, 4.8, // 8y
 3.9, 4.0,
 4.4, 5.5,
 4.6, 5.7,
 4.8, 5.9,
 4.0, 5.1,
-1.2,+4.0, // 8
-1.3,-3.2, // 9
 4.2, 5.3,
 4.3, 5.4,
 7.7, 7.8, // 9x
 7.9, 8.0,
 8.4, 8.5,
 8.6, 8.7,
 8.8, 8.9,
 9.0, 9.1,
-2.2,-1.3, // 8
+5.1,-4.3, // 9
 9.2, 9.3,
 9.3, 9.4,
 3.7, 4.8, // 9y
 3.9, 4.0,
 4.4, 5.5,
 4.6, 5.7,
 4.8, 5.9,
 4.0, 5.1,
-2.3,-3.2, // 8
-4.3,+4.1, // 9
 4.2, 5.3,
 4.3, 5.4,
};

// ----------------------------------------------------------------------
// Computed values
// ----------------------------------------------------------------------

const PylithScalar pylith::faults::CohesiveDynDataQuad4::_orientation[] = {
  0.0,  -1.0,  -1.0, 0.0,
  0.0,  -1.0,  -1.0, 0.0
};

const PylithScalar pylith::faults::CohesiveDynDataQuad4::_area[] = {
  1.0,
  1.0,
};

const PylithScalar pylith::faults::CohesiveDynDataQuad4::_initialTractions[] = {
  // Fault coordinate frame
  1.0, -2.0,
  1.1, -2.1,
};


const int pylith::faults::CohesiveDynDataQuad4::_numConstraintVert = 2;
const int pylith::faults::CohesiveDynDataQuad4::_constraintVertices[] = {
  11, 12
};
const int pylith::faults::CohesiveDynDataQuad4::_negativeVertices[] = {
   5,  6
};

// ----------------------------------------------------------------------
// Stick case
// ----------------------------------------------------------------------
// Input
const PylithScalar pylith::faults::CohesiveDynDataQuad4::_fieldIncrStick[] = {
  1.1, 2.1,
  1.2, 2.2,
  1.3, 2.3, // 4
  1.4, 2.4, // 5
  1.5, 2.5,
  1.6, 2.6,
  1.3, 2.7, // 8
  1.4, 2.9, // 9
  11.2, 2.20, // 10
  16.0, 5.40, // 11
};

// No change in slip
const PylithScalar pylith::faults::CohesiveDynDataQuad4::_slipStickE[] = {
 -0.4,  0.0,
 -0.5,  0.0,
};

// ----------------------------------------------------------------------
// Slip case
// ----------------------------------------------------------------------
// Input
const PylithScalar pylith::faults::CohesiveDynDataQuad4::_fieldIncrSlip[] = {
  1.1, 2.1,
  1.2, 2.2,
  1.3, 2.3, // 4
  1.4, 2.4, // 5
  1.5, 2.5,
  1.6, 2.6,
  1.3, 2.3, // 8
  1.4, 2.4, // 9
  1.8, 2.8, // 10
  1.0, 2.0, // 11
};

// Output
const PylithScalar pylith::faults::CohesiveDynDataQuad4::_fieldIncrSlipE[] = {
   1.100000000000,   2.100000000000,
   1.200000000000,   2.200000000000,
   1.300000000000,   2.300106970298,
   1.400000000000,   2.400969389712,
   1.500000000000,   2.500000000000,
   1.600000000000,   2.600000000000,
   1.300000000000,   2.299893029702,
   1.400000000000,   2.399030610288,
   1.800000000000,  -3.440000000000,
   1.000000000000,  -3.600000000000,
};

const PylithScalar pylith::faults::CohesiveDynDataQuad4::_slipSlipE[] = {
   0.000213940595,   0.000000000000,
   0.001938779425,   0.000000000000,
};

// ----------------------------------------------------------------------
// Open case
// ----------------------------------------------------------------------
// Input
const PylithScalar pylith::faults::CohesiveDynDataQuad4::_fieldIncrOpen[] = {
  1.1, 2.1,
  1.2, 2.2,
  1.3, 2.3, // 4
  1.4, 2.4, // 5
  1.5, 2.5,
  1.6, 2.6,
  1.3, 2.3, // 8
  1.4, 2.4, // 9
  -10.8, 2.8, // 10
  -10.0, 2.0, // 11
};

// Output
const PylithScalar pylith::faults::CohesiveDynDataQuad4::_fieldIncrOpenE[] = {
   1.100000000000,   2.100000000000,
   1.200000000000,   2.200000000000,
   1.302515834136,   2.299739404157,
   1.402887400343,   2.401420863468,
   1.500000000000,   2.500000000000,
   1.600000000000,   2.600000000000,
   1.297484165864,   2.300260595843,
   1.397112599657,   2.398579136532,
  -8.800000000000,  -9.800000000000,
  -8.000000000000,  -9.000000000000,
};

const PylithScalar pylith::faults::CohesiveDynDataQuad4::_slipOpenE[] = {
  -0.000521191686,   0.005031668272,
   0.002841726937,   0.005774800686,
};

// ----------------------------------------------------------------------
pylith::faults::CohesiveDynDataQuad4::CohesiveDynDataQuad4(void)
{ // constructor
  meshFilename = const_cast<char*>(_meshFilename);
  spaceDim = _spaceDim;
  cellDim = _cellDim;
  numBasis = _numBasis;
  numQuadPts = _numQuadPts;
  quadPts = const_cast<PylithScalar*>(_quadPts);
  quadWts = const_cast<PylithScalar*>(_quadWts);
  basis = const_cast<PylithScalar*>(_basis);
  basisDeriv = const_cast<PylithScalar*>(_basisDeriv);
  verticesRef = const_cast<PylithScalar*>(_verticesRef);
  id = _id;
  label = const_cast<char*>(_label);
  initialTractFilename = const_cast<char*>(_initialTractFilename);

  fieldT = const_cast<PylithScalar*>(_fieldT);
  jacobian = const_cast<PylithScalar*>(_jacobian);
  orientation = const_cast<PylithScalar*>(_orientation);
  area = const_cast<PylithScalar*>(_area);
  initialTractions = const_cast<PylithScalar*>(_initialTractions);

  constraintVertices = const_cast<int*>(_constraintVertices);
  negativeVertices = const_cast<int*>(_negativeVertices);
  numConstraintVert = _numConstraintVert;  

  // Stick
  fieldIncrStick = const_cast<PylithScalar*>(_fieldIncrStick);
  slipStickE = const_cast<PylithScalar*>(_slipStickE);

  // Slip
  fieldIncrSlip = const_cast<PylithScalar*>(_fieldIncrSlip);
  fieldIncrSlipE = const_cast<PylithScalar*>(_fieldIncrSlipE);
  slipSlipE = const_cast<PylithScalar*>(_slipSlipE);

  // Open
  fieldIncrOpen = const_cast<PylithScalar*>(_fieldIncrOpen);
  fieldIncrOpenE = const_cast<PylithScalar*>(_fieldIncrOpenE);
  slipOpenE = const_cast<PylithScalar*>(_slipOpenE);
} // constructor

pylith::faults::CohesiveDynDataQuad4::~CohesiveDynDataQuad4(void)
{}


// End of file
