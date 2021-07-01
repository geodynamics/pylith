// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University at Buffalo
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2021 University of California, Davis
//
// See LICENSE.md for license information.
//
// ======================================================================
//

/* Original mesh
 *
 * Cells are 0-3, vertices are 4-9.
 *
 *
 *         9
 *        / \
 *       /   \
 *      /     \
 *     /       \
 *    8---------5
 *     \       /|\
 *      \     / | \
 *       \   /  |  \
 *        \ /   |   \
 *         4    |    7
 *          \   |   /
 *           \  |  /
 *            \ | /
 *             \|/
 *              6
 *
 *
 * After adding cohesive elements
 *
 * Cells are 0-3, 4-5, vertices are 6-14,15-17.
 *
 *        11
 *        / \
 *       /   \
 *      /     \
 *     /       \
 *   10---------  7
 * 28 |        26/|
 *   14--------12 |
 *     \       /| |\
 *      \     / | | \
 *       \   /  | |  \
 *        \ /   | |   \
 *         6    | |    9
 *          \   | |   /
 *           \  | |  /
 *            \ | | /
 *             \| |/
 *             13-8
 *               27
 */


#include "CohesiveDynDataTri3d.hh"

const char* pylith::faults::CohesiveDynDataTri3d::_meshFilename =
  "data/tri3d.mesh";

const int pylith::faults::CohesiveDynDataTri3d::_spaceDim = 2;

const int pylith::faults::CohesiveDynDataTri3d::_cellDim = 1;

const int pylith::faults::CohesiveDynDataTri3d::_numBasis = 2;

const int pylith::faults::CohesiveDynDataTri3d::_numQuadPts = 2;

const PylithScalar pylith::faults::CohesiveDynDataTri3d::_quadPts[] = {
  -1.0, 1.0,
};

const PylithScalar pylith::faults::CohesiveDynDataTri3d::_quadWts[] = {
  1.0, 1.0
};

const PylithScalar pylith::faults::CohesiveDynDataTri3d::_basis[] = {
  1.0, 0.0,
  0.0, 1.0,
};

const PylithScalar pylith::faults::CohesiveDynDataTri3d::_basisDeriv[] = {
  -0.5, 0.5,
  -0.5, 0.5,
};

const PylithScalar pylith::faults::CohesiveDynDataTri3d::_verticesRef[] = {
  -1.0, 1.0
};

const int pylith::faults::CohesiveDynDataTri3d::_id = 10;

const char* pylith::faults::CohesiveDynDataTri3d::_label = "fault";

const char* pylith::faults::CohesiveDynDataTri3d::_initialTractFilename = 
  "data/tri3d_initialtract.spatialdb";

const PylithScalar pylith::faults::CohesiveDynDataTri3d::_fieldT[] = {
  6.1, 8.1,
  6.2, 8.2, // 7
  6.3, 8.3, // 8
  6.4, 8.4,
  6.5, 8.5, // 10
  6.6, 8.6,
  6.2, 8.2, // 12
  6.3, 8.3, // 13
  6.5, 8.5, // 14
 -3.8,-4.8, // 26
 -3.0, 4.0, // 27
  3.2,-4.2, // 28
};

const PylithScalar pylith::faults::CohesiveDynDataTri3d::_jacobian[] = {
  1.0, 1.1, // 6x
  1.1, 2.1,
  1.2, 2.2,
  1.3, 2.3,
  1.4, 2.4,
  1.5, 2.5,
  1.6, 2.6,
  1.7, 2.7,
  1.8, 2.8,
  1.9, 2.9,
  2.0, 3.0,
  2.1, 3.1,
  3.1, 1.0, // 6y
  3.2, 4.1,
  3.3, 4.2,
  3.4, 4.3,
  3.5, 4.4,
  3.6, 4.5,
  3.7, 4.6,
  3.8, 4.7,
  3.9, 4.8,
  3.0, 4.9,
  3.1, 4.0,
  3.2, 4.1,
  6.1, 7.2, // 7x
 +6.0,-1.0, // 7
 -1.1,-1.2, // 8
  6.8, 7.5,
 -1.3,-1.4, // 10
  6.6, 7.7,
  6.5, 7.8,
  6.4, 7.9,
  6.3, 7.0,
  6.5, 6.6,
  6.2, 7.1,
  6.1, 7.2,
  3.0, 5.1, // 7y
 -1.0,+6.1, // 7
 -0.9,-0.8, // 8
  3.3, 5.1,
 -0.7,-0.6, // 10
  3.5, 5.0,
  3.6, 5.4,
  3.7, 5.9,
  3.8, 5.5,
  3.9, 5.5,
  4.2, 5.6,
  4.1, 5.8,
  3.0, 2.7, // 8x
 -1.1,-0.9, // 7
 +6.2,-2.1, // 8
  3.7, 2.6,
 -2.2,-2.3, // 10
  3.5, 2.5,
  3.4, 2.1,
  3.3, 2.4,
  3.2, 2.2,
  3.1, 2.3,
  2.0, 2.0,
  3.0, 2.4,
  4.1, 6.1, // 8y
 -1.2,-0.8, // 7
 -2.1,+6.3, // 8
  4.3, 6.3,
 -1.4,-1.5, // 10
  4.3, 6.4,
  4.4, 6.4,
  4.4, 6.3,
  4.5, 6.5,
  4.5, 6.2,
  4.6, 7.0,
  4.6, 6.1,
  7.2, 8.7, // 9x
  7.9, 8.1,
  7.3, 8.6,
  7.1, 8.2,
  7.4, 8.5,
  7.2, 8.3,
  7.5, 8.4,
  7.3, 8.4,
  7.6, 8.3,
  7.4, 8.5,
  7.7, 8.6,
  7.8, 8.7,
  6.5, 3.3, // 9y
  6.4, 3.4,
  6.3, 3.8,
  6.2, 3.9,
  6.1, 3.8,
  6.0, 3.7,
  6.9, 3.6,
  5.7, 3.5,
  5.6, 3.4,
  5.5, 3.3,
  5.4, 3.2,
  5.5, 3.3,
  4.6, 6.6, // 10x
 -1.3,-0.7, // 7
 -2.2,-1.4, // 8
  5.8, 9.8,
 +6.4,-1.1, // 10
  5.9, 9.9,
  5.1, 9.7,
  5.2, 9.1,
  5.3, 9.6,
  5.1, 9.2,
  5.4, 9.5,
  5.2, 9.3,
  5.5, 9.4, // 10y
 -1.4,-0.6, // 7
 -2.3,-1.5, // 8
  6.4, 9.5,
 -1.1,+6.5, // 10
  6.5, 9.6,
  6.8, 9.1,
  6.6, 9.7,
  7.9, 8.9,
  7.7, 8.8,
  7.1, 8.8,
  7.8, 8.9,
  4.4, 2.0, // 11x
  4.3, 2.1,
  4.2, 2.2,
  4.1, 2.3,
  4.0, 2.3,
  4.9, 2.4,
  4.8, 2.5,
  3.7, 2.6,
  3.6, 2.8,
  3.5, 2.9,
  3.4, 2.1,
  3.3, 2.9,
  2.2, 2.8, // 11y
  2.1, 3.7,
  2.0, 3.6,
  2.1, 3.5,
  2.2, 3.4,
  2.3, 3.3,
  2.4, 3.1,
  2.5, 3.0,
  2.6, 3.9,
  2.7, 3.8,
  2.8, 3.7,
  2.9, 3.6,
  6.1, 7.2, // 12x
  6.8, 7.5,
  6.6, 7.7,
  6.5, 7.8,
  6.4, 7.9,
  6.3, 7.0,
 +5.0,-1.0, // 12
 -1.1,-1.2, // 13
 -1.3,-1.4, // 14
  6.5, 6.6,
  6.2, 7.1,
  6.1, 7.2,
  3.0, 5.1, // 12y
  3.3, 5.1,
  3.5, 5.0,
  3.6, 5.4,
  3.7, 5.9,
  3.8, 5.5,
 -1.0,+5.1, // 12
 -0.9,-0.8, // 13
 -0.7,-0.6, // 14
  3.9, 5.5,
  4.2, 5.6,
  4.1, 5.8,
  3.0, 2.7, // 13x
  3.7, 2.6,
  3.5, 2.5,
  3.4, 2.1,
  3.3, 2.4,
  3.2, 2.2,
 -1.1,-0.9, // 12
 +5.2,-2.1, // 13
 -2.2,-2.3, // 14
  3.1, 2.3,
  2.0, 2.0,
  3.0, 2.4,
  4.1, 6.1, // 13y
  4.3, 6.3,
  4.3, 6.4,
  4.4, 6.4,
  4.4, 6.3,
  4.5, 6.5,
 -1.2,-0.8, // 12
 -2.1,+5.3, // 13
 -1.4,-1.5, // 14
  4.5, 6.2,
  4.6, 7.0,
  4.6, 6.1,
  4.6, 6.6, // 14x
  5.8, 9.8,
  5.9, 9.9,
  5.1, 9.7,
  5.2, 9.1,
  5.3, 9.6,
 -1.3,-0.7, // 12
 -2.2,-1.4, // 13
 +5.4,-1.1, // 14
  5.1, 9.2,
  5.4, 9.5,
  5.2, 9.3,
  5.5, 9.4, // 14y
  6.4, 9.5,
  6.5, 9.6,
  6.8, 9.1,
  6.6, 9.7,
  7.9, 8.9,
 -1.4,-0.6, // 12
 -2.3,-1.5, // 13
 -1.1,+5.5, // 14
  7.7, 8.8,
  7.1, 8.8,
  7.8, 8.9,
  3.2, 8.3, // 26x
  3.3, 8.4,
  5.4, 9.3,
  5.6, 9.7,
  3.7, 9.0,
  5.9, 9.9,
  6.0, 9.8,
  4.4, 4.8,
  4.6, 4.7,
  4.8, 4.6,
  4.9, 4.4,
  4.0, 4.2,
  4.2, 4.3, // 26y
  4.3, 4.4,
  7.5, 3.4,
  6.7, 3.5,
  6.4, 3.6,
  4.6, 3.9,
  4.7, 4.0,
  8.9, 2.8,
  7.6, 2.7,
  6.4, 2.6,
  5.3, 2.5,
  3.8, 2.3,
  4.5, 2.2, // 27x
  8.5, 2.4,
  0.0, 1.0,
  7.4, 3.6,
  6.6, 3.5,
  4.7, 3.4,
  3.8, 3.5,
  0.0, 1.0,
  5.9, 3.7,
  8.7, 4.6,
  7.6, 4.5,
  6.5, 4.4,
  5.5, 4.3, // 27y
  4.3, 4.8,
  1.0, 0.0,
  4.3, 4.7,
  6.5, 4.6,
  9.6, 4.5,
  8.7, 4.3,
  1.0, 0.0,
  7.9, 4.5,
  6.7, 5.3,
  5.6, 5.8,
  3.5, 5.7,
  2.4, 5.6, // 28x
  3.3, 5.5,
  4.2, 5.3,
  5.4, 5.6,
  1.0, 6.0,
  0.8, 6.6,
  9.8, 6.5,
  8.5, 6.5,
  1.0, 0.0,
  7.5, 7.3,
  6.4, 7.6,
  5.2, 7.8,
  4.5, 7.7, // 28y
  3.7, 8.6,
  3.8, 8.5,
  2.9, 8.3,
  0.0, 1.0,
  2.0, 9.9,
  2.7, 9.8,
  1.6, 9.7,
  0.0, 1.0,
  1.5, 5.5,
  1.2, 5.4,
  1.1, 5.3,
};


// ----------------------------------------------------------------------
// Computed values
// ----------------------------------------------------------------------

const PylithScalar pylith::faults::CohesiveDynDataTri3d::_orientation[] = {
  -0.70710678118654757, +0.70710678118654757,  
  +0.70710678118654757, +0.70710678118654757,
  0.0, +1.0,  +1.0,  0.0,
 -1.0,  0.0,   0.0, +1.0
};

const PylithScalar pylith::faults::CohesiveDynDataTri3d::_area[] = {
  2.0,
  1.0,
  1.0,
};

const int pylith::faults::CohesiveDynDataTri3d::_numConstraintEdges = 3;
const int pylith::faults::CohesiveDynDataTri3d::_constraintEdges[] = {
  26, 27, 28
};
const int pylith::faults::CohesiveDynDataTri3d::_negativeVertices[] = {
   7,  8, 10
};

const PylithScalar pylith::faults::CohesiveDynDataTri3d::_initialTractions[] = {
  // Fault coordinate frame
  1.0, -2.0,
  1.1, -2.1,
  1.2, -2.2,
};


// ----------------------------------------------------------------------
// Stick case
// ----------------------------------------------------------------------
// Input
const PylithScalar pylith::faults::CohesiveDynDataTri3d::_fieldIncrStick[] = {
  1.1, 2.1,
  1.2, 2.2, // 5
  1.3, 2.3, // 6
  1.4, 2.4,
  1.5, 2.5, // 8
  1.6, 2.6,
  1.2, 2.2, // 10
  1.3, 2.3, // 11
  1.5, 2.5, // 12
 -21.8,-22.8, // 26
 -21.0, 2.0, // 27
  2.2,-22.2, // 28
};

// No slip
const PylithScalar pylith::faults::CohesiveDynDataTri3d::_slipStickE[] = {
 0.0,  0.0,
 0.0,  0.0,
 0.0,  0.0,
};

// ----------------------------------------------------------------------
// Slip case
// ----------------------------------------------------------------------
// Input
const PylithScalar pylith::faults::CohesiveDynDataTri3d::_fieldIncrSlip[] = {
  1.1, 2.1,
  1.2, 2.2, // 5
  1.3, 2.3, // 6
  1.4, 2.4,
  1.5, 2.5, // 8
  1.6, 2.6,
  1.2, 2.2, // 10
  1.3, 2.3, // 11
  1.5, 2.5, // 12
 -1.8,+3.6, // 26
 -1.0, 1.1, // 27
  1.7,-1.2, // 28
};

// Output
const PylithScalar pylith::faults::CohesiveDynDataTri3d::_fieldIncrSlipE[] = {
   1.100000000000,   2.100000000000,
   1.199970441179,   2.200029558821,
   1.300000000000,   2.299168310929,
   1.400000000000,   2.400000000000,
   1.499489949674,   2.500000000000,
   1.600000000000,   2.600000000000,
   1.200029558821,   2.199970441179,
   1.300000000000,   2.300831689071,
   1.500510050326,   2.500000000000,
  -1.640000000000,   3.440000000000,
  -1.000000000000,  -1.600000000000,
   0.040000000000,  -1.200000000000,
};

const PylithScalar pylith::faults::CohesiveDynDataTri3d::_slipSlipE[] = {
  -0.000083604971,   0.000000000000,
   0.001663378142,   0.000000000000,
  -0.001020100652,   0.000000000000,
};

// ----------------------------------------------------------------------
// Open case
// ----------------------------------------------------------------------
// Input
const PylithScalar pylith::faults::CohesiveDynDataTri3d::_fieldIncrOpen[] = {
  1.1, 2.1,
  1.2, 2.2, // 5
  1.3, 2.3, // 6
  1.4, 2.4,
  1.5, 2.5, // 8
  1.6, 2.6,
  1.2, 2.2, // 10
  1.3, 2.3, // 11
  1.5, 2.5, // 12
+11.8, 11.8, // 26
+10.0, 0.1, // 27
  1.2,+10.2, // 28
};

// Output
const PylithScalar pylith::faults::CohesiveDynDataTri3d::_fieldIncrOpenE[] = {
   1.100000000000,   2.100000000000,
   1.190062443638,   2.192265294477,
   1.292981353233,   2.293412522606,
   1.400000000000,   2.400000000000,
   1.495037703623,   2.494851226153,
   1.600000000000,   2.600000000000,
   1.209937556362,   2.207734705523,
   1.307018646767,   2.306587477394,
   1.504962296377,   2.505148773847,
   3.800000000000,   4.800000000000,
   3.000000000000,  -4.000000000000,
  -3.200000000000,   4.200000000000,
};

const PylithScalar pylith::faults::CohesiveDynDataTri3d::_slipOpenE[] = {
  -0.003115301533,   0.024992352435,
   0.013174954788,   0.014037293534,
  -0.009924592753,   0.010297547694,
};

// ----------------------------------------------------------------------
pylith::faults::CohesiveDynDataTri3d::CohesiveDynDataTri3d(void)
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

  constraintEdges = const_cast<int*>(_constraintEdges);
  negativeVertices = const_cast<int*>(_negativeVertices);
  numConstraintEdges = _numConstraintEdges;  

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

pylith::faults::CohesiveDynDataTri3d::~CohesiveDynDataTri3d(void)
{}


// End of file
