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
 * 17 |        15/|
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
 *               16
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
  6.2, 8.2, // 5
  6.3, 8.3, // 6
  6.4, 8.4,
  6.5, 8.5, // 8
  6.6, 8.6,
  6.2, 8.2, // 10
  6.3, 8.3, // 11
  6.5, 8.5, // 12
 -3.8, 4.8, // 13
  3.0, 4.0, // 14
  3.2, 4.2, // 15
};

const PylithScalar pylith::faults::CohesiveDynDataTri3d::_jacobian[] = {
  1.0, 1.1, // 4x
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
  3.1, 1.0, // 4y
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
  6.1, 7.2, // 5x
 +6.0,-1.0, // 5
 -1.1,-1.2, // 6
  6.8, 7.5,
 -1.3,-1.4, // 8
  6.6, 7.7,
  6.5, 7.8,
  6.4, 7.9,
  6.3, 7.0,
  6.5, 6.6,
  6.2, 7.1,
  6.1, 7.2,
  3.0, 5.1, // 5y
 -1.0,+6.1, // 5
 -0.9,-0.8, // 6
  3.3, 5.1,
 -0.7,-0.6, // 8
  3.5, 5.0,
  3.6, 5.4,
  3.7, 5.9,
  3.8, 5.5,
  3.9, 5.5,
  4.2, 5.6,
  4.1, 5.8,
  3.0, 2.7, // 6x
 -1.1,-0.9, // 5
 +6.2,-2.1, // 6
  3.7, 2.6,
 -2.2,-2.3, // 8
  3.5, 2.5,
  3.4, 2.1,
  3.3, 2.4,
  3.2, 2.2,
  3.1, 2.3,
  2.0, 2.0,
  3.0, 2.4,
  4.1, 6.1, // 6y
 -1.2,-0.8, // 5
 -2.1,+6.3, // 6
  4.3, 6.3,
 -1.4,-1.5, // 8
  4.3, 6.4,
  4.4, 6.4,
  4.4, 6.3,
  4.5, 6.5,
  4.5, 6.2,
  4.6, 7.0,
  4.6, 6.1,
  7.2, 8.7, // 7x
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
  6.5, 3.3, // 7y
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
  4.6, 6.6, // 8x
 -1.3,-0.7, // 5
 -2.2,-1.4, // 6
  5.8, 9.8,
 +6.4,-1.1, // 8
  5.9, 9.9,
  5.1, 9.7,
  5.2, 9.1,
  5.3, 9.6,
  5.1, 9.2,
  5.4, 9.5,
  5.2, 9.3,
  5.5, 9.4, // 8y
 -1.4,-0.6, // 5
 -2.3,-1.5, // 6
  6.4, 9.5,
 -1.1,+6.5, // 8
  6.5, 9.6,
  6.8, 9.1,
  6.6, 9.7,
  7.9, 8.9,
  7.7, 8.8,
  7.1, 8.8,
  7.8, 8.9,
  4.4, 2.0, // 9x
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
  2.2, 2.8, // 9y
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
  6.1, 7.2, // 10x
  6.8, 7.5,
  6.6, 7.7,
  6.5, 7.8,
  6.4, 7.9,
  6.3, 7.0,
 +5.0,-1.0, // 10
 -1.1,-1.2, // 11
 -1.3,-1.4, // 12
  6.5, 6.6,
  6.2, 7.1,
  6.1, 7.2,
  3.0, 5.1, // 10y
  3.3, 5.1,
  3.5, 5.0,
  3.6, 5.4,
  3.7, 5.9,
  3.8, 5.5,
 -1.0,+5.1, // 10
 -0.9,-0.8, // 11
 -0.7,-0.6, // 12
  3.9, 5.5,
  4.2, 5.6,
  4.1, 5.8,
  3.0, 2.7, // 11x
  3.7, 2.6,
  3.5, 2.5,
  3.4, 2.1,
  3.3, 2.4,
  3.2, 2.2,
 -1.1,-0.9, // 10
 +5.2,-2.1, // 11
 -2.2,-2.3, // 12
  3.1, 2.3,
  2.0, 2.0,
  3.0, 2.4,
  4.1, 6.1, // 11y
  4.3, 6.3,
  4.3, 6.4,
  4.4, 6.4,
  4.4, 6.3,
  4.5, 6.5,
 -1.2,-0.8, // 10
 -2.1,+5.3, // 11
 -1.4,-1.5, // 12
  4.5, 6.2,
  4.6, 7.0,
  4.6, 6.1,
  4.6, 6.6, // 12x
  5.8, 9.8,
  5.9, 9.9,
  5.1, 9.7,
  5.2, 9.1,
  5.3, 9.6,
 -1.3,-0.7, // 10
 -2.2,-1.4, // 11
 +5.4,-1.1, // 12
  5.1, 9.2,
  5.4, 9.5,
  5.2, 9.3,
  5.5, 9.4, // 12y
  6.4, 9.5,
  6.5, 9.6,
  6.8, 9.1,
  6.6, 9.7,
  7.9, 8.9,
 -1.4,-0.6, // 10
 -2.3,-1.5, // 11
 -1.1,+5.5, // 12
  7.7, 8.8,
  7.1, 8.8,
  7.8, 8.9,
  3.2, 8.3, // 13x
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
  4.2, 4.3, // 13y
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
  4.5, 2.2, // 14x
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
  5.5, 4.3, // 14y
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
  2.4, 5.6, // 15x
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
  4.5, 7.7, // 15y
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
  +0.70710678118654757, -0.70710678118654757,  
  -0.70710678118654757, -0.70710678118654757,
  0.0, -1.0,  -1.0,  0.0,
 +1.0,  0.0,   0.0, -1.0
};

const PylithScalar pylith::faults::CohesiveDynDataTri3d::_area[] = {
  2.0,
  1.0,
  1.0,
};

const int pylith::faults::CohesiveDynDataTri3d::_numConstraintVert = 3;
const int pylith::faults::CohesiveDynDataTri3d::_constraintVertices[] = {
  15, 16, 17
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
  21.8, 22.8, // 13
  21.0, 2.0, // 14
  2.2, 22.2, // 15
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
  1.8, 0.8, // 13
  1.0, 0.1, // 14
  1.2, 0.2, // 15
};

// Output
const PylithScalar pylith::faults::CohesiveDynDataTri3d::_fieldIncrSlipE[] = {
   1.100000000000,   2.100000000000,
   1.200805761009,   2.199194238991,
   1.300000000000,   2.299628356338,
   1.400000000000,   2.400000000000,
   1.499716364403,   2.500000000000,
   1.600000000000,   2.600000000000,
   1.199194238991,   2.200805761009,
   1.300000000000,   2.300371643662,
   1.500283635597,   2.500000000000,
   4.520000000000,  -1.920000000000,
   1.000000000000,  -1.600000000000,
  -0.560000000000,   0.200000000000,
};

const PylithScalar pylith::faults::CohesiveDynDataTri3d::_slipSlipE[] = {
  -0.002279036295,   0.000000000000,
  -0.000743287324,   0.000000000000,
   0.000567271195,   0.000000000000,
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
-10.8, 0.8, // 13
-10.0, 0.1, // 14
  1.2,-10.2, // 15
};

// Output
const PylithScalar pylith::faults::CohesiveDynDataTri3d::_fieldIncrOpenE[] = {
   1.100000000000,   2.100000000000,
   1.208740107962,   2.201462120349,
   1.304388525058,   2.303014075601,
   1.400000000000,   2.400000000000,
   1.502122067204,   2.503667076610,
   1.600000000000,   2.600000000000,
   1.191259892038,   2.198537879651,
   1.295611474942,   2.296985924399,
   1.497877932796,   2.496332923390,
   3.800000000000,  -4.800000000000,
  -3.000000000000,  -4.000000000000,
  -3.200000000000,  -4.200000000000,
};

const PylithScalar pylith::faults::CohesiveDynDataTri3d::_slipOpenE[] = {
  -0.010292628789,   0.014428129643,
   0.006028151203,   0.008777050116,
  -0.004244134409,   0.007334153220,
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

  constraintVertices = const_cast<int*>(_constraintVertices);
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

pylith::faults::CohesiveDynDataTri3d::~CohesiveDynDataTri3d(void)
{}


// End of file
