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
 * Cells are 0-1,10 vertices are 2-9.
 *
 *       3 -------- 5 -11-- 9 -------- 7
 *       |          |       |          |
 *       |          |       |          |
 *       |          |       |          |
 *       |          |       |          |
 *       |          |       |          |
 *       |          |       |          |
 *       |          |       |          |
 *       |          |       |          |
 *       2 -------- 4 -10-- 8 -------- 6
 */

#include "CohesiveDynDataQuad4.hh"

const char* pylith::faults::CohesiveDynDataQuad4::_meshFilename =
  "data/quad4.mesh";

const int pylith::faults::CohesiveDynDataQuad4::_spaceDim = 2;

const int pylith::faults::CohesiveDynDataQuad4::_cellDim = 1;

const int pylith::faults::CohesiveDynDataQuad4::_numBasis = 2;

const int pylith::faults::CohesiveDynDataQuad4::_numQuadPts = 1;

const double pylith::faults::CohesiveDynDataQuad4::_quadPts[] = {
  0.0,
};

const double pylith::faults::CohesiveDynDataQuad4::_quadWts[] = {
  2.0,
};

const double pylith::faults::CohesiveDynDataQuad4::_basis[] = {
  0.5,
  0.5
};

const double pylith::faults::CohesiveDynDataQuad4::_basisDeriv[] = {
  -0.5,
   0.5
};

const double pylith::faults::CohesiveDynDataQuad4::_verticesRef[] = {
  -1.0, 1.0
};

const int pylith::faults::CohesiveDynDataQuad4::_id = 10;

const char* pylith::faults::CohesiveDynDataQuad4::_label = "fault";

const char* pylith::faults::CohesiveDynDataQuad4::_initialTractFilename = 
  "data/quad4_initialtract.spatialdb";

const double pylith::faults::CohesiveDynDataQuad4::_fieldT[] = {
  8.1, 9.1,
  8.2, 9.2,
  8.3, 9.3, // 4
  8.4, 9.4, // 5
  8.5, 9.5,
  8.6, 9.6,
  8.7, 9.7, // 8
  8.9, 9.9, // 9
  8.8, 9.8, // 10
  8.0, 9.0, // 11
};

// :TODO: Make sensible values for Jacobian for DOF on positive and
// negative sides of the fault. Add semi-random values for other DOF.
const double pylith::faults::CohesiveDynDataQuad4::_jacobian[] = {
   1, 0.1,
 0.2, 0.3,
 0.4, 0.5,
 0.6, 0.7,
 0.8, 0.9,
   1, 1.1,
 1.2, 1.3,
 1.4, 1.5,
 1.6, 1.7,
 1.8, 1.9,
   2,   1,
 2.1, 2.2,
 2.3, 2.4,
 2.5, 2.6,
 2.7, 2.8,
 2.9,   3,
 3.1, 3.2,
 3.3, 3.4,
 3.5, 3.6,
 3.7, 3.8,
 3.9,   4,
   1, 4.1,
 4.2, 4.3,
 4.4, 4.5,
 4.6, 4.7,
 4.8, 4.9,
   5, 5.1,
 5.2, 5.3,
 5.4, 5.5,
 5.6, 5.7,
 5.8, 5.9,
   6,   1,
 6.1, 6.2,
 6.3, 6.4,
 6.5, 6.6,
 6.7, 6.8,
 6.9,   7,
 7.1, 7.2,
 7.3, 7.4,
 7.5, 7.6,
 7.7, 7.8,
 7.9,   8,
   1, 8.1,
 8.2, 8.3,
 8.4, 8.5,
 8.6, 8.7,
 8.8, 8.9,
   9, 9.1,
 9.2,  -1,
 9.3, 9.4,
 9.5, 9.6,
 9.7, 9.8,
 9.9,   1,
  10,10.1,
10.2,10.3,
10.4,10.5,
10.6,10.7,
10.8,10.9,
  -1,  11,
11.1,11.2,
11.3,11.4,
11.5,11.6,
11.7,11.8,
   1,11.9,
  12,12.1,
12.2,12.3,
12.4,12.5,
12.6,12.7,
12.8,12.9,
  13,  -1,
13.1,13.2,
13.3,13.4,
13.5,13.6,
13.7,   1,
13.8,13.9,
  14,14.1,
14.2,14.3,
14.4,14.5,
14.6,14.7,
  -1,14.8,
14.9,  15,
15.1,15.2,
15.3,15.4,
15.5,15.6,
   1,15.7,
15.8,15.9,
  16,16.1,
16.2,16.3,
16.4,16.5,
16.6,16.7,
16.8,16.9,
  17,17.1,
17.2,17.3,
17.4,17.5,
17.6,   1,
17.7,17.8,
17.9,  18,
18.1,18.2,
18.3,18.4,
18.5,18.6,
18.7,18.8,
18.9,  19,
19.1,19.2,
19.3,19.4,
19.5,19.6,
   1,19.7,
19.8,19.9,
  20,20.1,
20.2,20.3,
20.4,20.5,
20.6,20.7,
20.8,20.9,
  21,21.1,
21.2,21.3,
21.4,21.5,
21.6,   1,
21.7,21.8,
21.9,  22,
22.1,22.2,
22.3,22.4,
22.5,22.6,
22.7,22.8,
22.9,  23,
23.1,23.2,
23.3,23.4,
23.5,23.6,
   1,23.7,
23.8,23.9,
  24,   1,
24.1,24.2,
24.3,24.4,
24.5,24.6,
24.7,24.8,
24.9,  25,
25.1,25.2,
25.3,25.4,
25.5,   1,
25.6,25.7,
   1,25.8,
25.9,  26,
26.1,26.2,
26.3,26.4,
26.5,26.6,
26.7,26.8,
26.9,  27,
27.1,27.2,
27.3,27.4,
   1,27.5,
27.6,27.7,
27.8,   1,
27.9,  28,
28.1,28.2,
28.3,28.4,
28.5,28.6,
28.7,28.8,
28.9,  29,
29.1,29.2,
29.3,   1,
29.4,29.5,
   1,29.6,
29.7,29.8,
29.9,  30,
30.1,  -1,
30.2,30.3,
30.4,30.5,
30.6,30.7,
30.8,   1,
30.9,  31,
31.1,31.2,
31.3,31.4,
31.5,31.6,
31.7,31.8,
  -1,31.9,
  32,32.1,
32.2,32.3,
32.4,32.5,
   1,32.6,
32.7,32.8,
32.9,  33,
33.1,33.2,
33.3,33.4,
33.5,33.6,
33.7,33.8,
33.9,  -1,
  34,34.1,
34.2,34.3,
34.4,34.5,
34.6,   1,
34.7,34.8,
34.9,  35,
35.1,35.2,
35.3,35.4,
35.5,35.6,
  -1,35.7,
35.8,35.9,
  36,36.1,
36.2,36.3,
   1,36.4,
36.5,36.6,
36.7,36.8,
};

// ----------------------------------------------------------------------
// Computed values
// ----------------------------------------------------------------------

const double pylith::faults::CohesiveDynDataQuad4::_orientation[] = {
  0.0,  -1.0,  -1.0, 0.0,
  0.0,  -1.0,  -1.0, 0.0
};

const double pylith::faults::CohesiveDynDataQuad4::_area[] = {
  1.0,
  1.0,
};

const double pylith::faults::CohesiveDynDataQuad4::_initialTractions[] = {
  2.0, -1.0,
  2.1, -1.1,
};


const int pylith::faults::CohesiveDynDataQuad4::_numConstraintVert = 2;
const int pylith::faults::CohesiveDynDataQuad4::_constraintVertices[] = {
  10, 11
};

// ----------------------------------------------------------------------
// Stick case
// ----------------------------------------------------------------------
// Input
const double pylith::faults::CohesiveDynDataQuad4::_fieldIncrStick[] = {
  1.1, 2.1,
  1.2, 2.2,
  1.3, 2.3, // 4
  1.4, 2.4, // 5
  1.5, 2.5,
  1.6, 2.6,
  1.7, 2.7, // 8
  1.9, 2.9, // 9
  21.8, 2.8, // 10
  21.0, 2.0, // 11
};

// No change in fieldIncr
// Zero slip

// ----------------------------------------------------------------------
// Slip case
// ----------------------------------------------------------------------
// Input
const double pylith::faults::CohesiveDynDataQuad4::_fieldIncrSlip[] = {
  9.1, 10.1,
  9.2, 10.2,
  9.3, 10.3, // 4
  9.4, 10.4, // 5
  9.5, 10.5,
  9.6, 10.6,
  9.7, 10.7, // 8
  9.9, 10.9, // 9
  -9.8, -10.8, // 10
  -9.0, -10.0, // 11
};

// Output
// :TODO: Update Lagrange multiplier values
const double pylith::faults::CohesiveDynDataQuad4::_fieldIncrSlipE[] = {
  9.1,  10.100000000000,
  9.2,  10.200000000000,
  9.3,  10.313079001869,
  9.4,  10.405865253910,
  9.5,  10.500000000000,
  9.6,  10.600000000000,
  9.7,  10.686920998131,
  9.9,  10.894134746090,
 -9.4, -10.800000000000,
 -8.6, -10.000000000000,
};

// Update slip values based on changes in Lagrange multiplier values
const double pylith::faults::CohesiveDynDataQuad4::_slipSlipE[] = {
  0.026158003737426,                   0.0,
  0.011730507820273,                   0.0,
};

// ----------------------------------------------------------------------
// Open case
// ----------------------------------------------------------------------
// Input
const double pylith::faults::CohesiveDynDataQuad4::_fieldIncrOpen[] = {
  9.1, 10.1,
  9.2, 10.2,
  9.3, 10.3, // 4
  9.4, 10.4, // 5
  9.5, 10.5,
  9.6, 10.6,
  9.7, 10.7, // 8
  9.9, 10.9, // 9
  9.8, 10.8, // 10
  9.0, 10.0, // 11
};

// Output
const double pylith::faults::CohesiveDynDataQuad4::_fieldIncrOpenE[] = {
   9.100000000000,  10.100000000000,
   9.200000000000,  10.200000000000,
   9.300000000000,  10.685047196672,
   9.932709332305,  11.176949275526,
   9.500000000000,  10.500000000000,
   9.600000000000,  10.600000000000,
   9.700000000000,  10.314952803328,
   9.367290667695,  10.123050724474,
  -8.800000000000,  -9.800000000000,
  -8.000000000000,  -9.000000000000,
};

const double pylith::faults::CohesiveDynDataQuad4::_slipOpenE[] = {
  0.770094393343286,  0.0,
  1.553898551051246,   1.065418664609568,
};

// ----------------------------------------------------------------------
pylith::faults::CohesiveDynDataQuad4::CohesiveDynDataQuad4(void)
{ // constructor
  meshFilename = const_cast<char*>(_meshFilename);
  spaceDim = _spaceDim;
  cellDim = _cellDim;
  numBasis = _numBasis;
  numQuadPts = _numQuadPts;
  quadPts = const_cast<double*>(_quadPts);
  quadWts = const_cast<double*>(_quadWts);
  basis = const_cast<double*>(_basis);
  basisDeriv = const_cast<double*>(_basisDeriv);
  verticesRef = const_cast<double*>(_verticesRef);
  id = _id;
  label = const_cast<char*>(_label);
  initialTractFilename = const_cast<char*>(_initialTractFilename);

  fieldT = const_cast<double*>(_fieldT);
  jacobian = const_cast<double*>(_jacobian);
  orientation = const_cast<double*>(_orientation);
  area = const_cast<double*>(_area);
  initialTractions = const_cast<double*>(_initialTractions);

  constraintVertices = const_cast<int*>(_constraintVertices);
  numConstraintVert = _numConstraintVert;  

  // Stick
  fieldIncrStick = const_cast<double*>(_fieldIncrStick);

  // Slip
  fieldIncrSlip = const_cast<double*>(_fieldIncrSlip);
  fieldIncrSlipE = const_cast<double*>(_fieldIncrSlipE);
  slipSlipE = const_cast<double*>(_slipSlipE);

  // Open
  fieldIncrOpen = const_cast<double*>(_fieldIncrOpen);
  fieldIncrOpenE = const_cast<double*>(_fieldIncrOpenE);
  slipOpenE = const_cast<double*>(_slipOpenE);
} // constructor

pylith::faults::CohesiveDynDataQuad4::~CohesiveDynDataQuad4(void)
{}


// End of file
