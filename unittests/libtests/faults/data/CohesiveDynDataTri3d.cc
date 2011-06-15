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
 * Cells are 0-3, 13-14, vertices are 4-12.
 *
 *         9
 *        / \
 *       /   \
 *      /     \
 *     /       \
 *    8---------  5
 * 15 |        13/|
 *   12--------10 |
 *     \       /| |\
 *      \     / | | \
 *       \   /  | |  \
 *        \ /   | |   \
 *         4    | |    7
 *          \   | |   /
 *           \  | |  /
 *            \ | | /
 *             \| |/
 *             11-6
 *               14
 */


#include "CohesiveDynDataTri3d.hh"

const char* pylith::faults::CohesiveDynDataTri3d::_meshFilename =
  "data/tri3d.mesh";

const int pylith::faults::CohesiveDynDataTri3d::_spaceDim = 2;

const int pylith::faults::CohesiveDynDataTri3d::_cellDim = 1;

const int pylith::faults::CohesiveDynDataTri3d::_numBasis = 2;

const int pylith::faults::CohesiveDynDataTri3d::_numQuadPts = 1;

const double pylith::faults::CohesiveDynDataTri3d::_quadPts[] = {
  0.0,
};

const double pylith::faults::CohesiveDynDataTri3d::_quadWts[] = {
  2.0,
};

const double pylith::faults::CohesiveDynDataTri3d::_basis[] = {
  0.5,
  0.5
};

const double pylith::faults::CohesiveDynDataTri3d::_basisDeriv[] = {
  -0.5,
   0.5
};

const double pylith::faults::CohesiveDynDataTri3d::_verticesRef[] = {
  -1.0, 1.0
};

const int pylith::faults::CohesiveDynDataTri3d::_id = 10;

const char* pylith::faults::CohesiveDynDataTri3d::_label = "fault";

const char* pylith::faults::CohesiveDynDataTri3d::_initialTractFilename = 
  "data/tri3d_initialtract.spatialdb";

const double pylith::faults::CohesiveDynDataTri3d::_fieldT[] = {
  6.1, 8.1,
  6.2, 8.2, // 5
  6.3, 8.3, // 6
  6.4, 8.4,
  6.5, 8.5, // 8
  6.6, 8.6,
  6.7, 8.7, // 10
  6.9, 8.9, // 11
  7.1, 9.1, // 12
  6.8, 8.8, // 13
  6.0, 8.0, // 14
  7.2, 9.2, // 15
};

// :TODO: Make sensible values for Jacobian for DOF on positive and
// negative sides of the fault. Add semi-random values for other DOF.
const double pylith::faults::CohesiveDynDataTri3d::_jacobian[] = {
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
  6.0, 7.3,
  6.9, 7.4,
  6.8, 7.5,
  6.7, 7.6,
  6.6, 7.7,
  6.5, 7.8,
  6.4, 7.9,
  6.3, 7.0,
 -0.70710678118654757, +0.70710678118654757, // 13
  6.2, 7.1,
  6.1, 7.2,
  3.0, 5.1, // 5y
  3.1, 5.2,
  3.2, 5.2,
  3.3, 5.1,
  3.4, 5.3,
  3.5, 5.0,
  3.6, 5.4,
  3.7, 5.9,
  3.8, 5.5,
 +0.70710678118654757, +0.70710678118654757, // 13
  4.2, 5.6,
  4.1, 5.8,
  3.0, 2.7, // 6x
  3.9, 2.7,
  3.8, 2.8,
  3.7, 2.6,
  3.6, 2.9,
  3.5, 2.5,
  3.4, 2.1,
  3.3, 2.4,
  3.2, 2.2,
  3.1, 2.3,
  0.0,+1.0, // 14
  3.0, 2.4,
  4.1, 6.1, // 6y
  4.1, 6.6,
  4.2, 6.2,
  4.2, 6.5,
  4.3, 6.3,
  4.3, 6.4,
  4.4, 6.4,
  4.4, 6.3,
  4.5, 6.5,
  4.5, 6.2,
 +1.0, 0.0, // 14
  4.6, 6.1,
  4.6, 6.6, // 7x
  5.7, 9.7,
  5.7, 9.9,
  5.8, 9.8,
  5.8, 9.8,
  5.9, 9.9,
  5.1, 9.7,
  5.2, 9.1,
  5.3, 9.6,
  5.1, 9.2,
  5.4, 9.5,
  5.2, 9.3,
  5.5, 9.4, // 7y
  6.3, 9.4,
  6.6, 9.3,
  6.4, 9.5,
  6.7, 9.2,
  6.5, 9.6,
  6.8, 9.1,
  6.6, 9.7,
  7.9, 8.9,
  7.7, 8.8,
  7.1, 8.8,
  7.8, 8.9,
  7.2, 8.7, // 8x
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
 -1.0, 0.0, // 15
  6.5, 3.3, // 8y
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
  0.0,+1.0, // 15
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
  1.0, 3.5, // 10x
  1.1, 4.4,
  1.2, 4.3,
  1.3, 4.2,
  1.4, 4.1,
  1.5, 4.0,
  1.6, 4.1,
  1.7, 4.2,
  1.8, 4.3,
 +0.70710678118654757, -0.70710678118654757, // 13
  1.9, 5.4,
  9.0, 5.5,
  8.0, 4.7, // 10y
  7.1, 4.5,
  6.2, 4.4,
  6.3, 4.6,
  6.4, 4.7,
  6.5, 4.4,
  4.6, 4.8,
  4.7, 4.4,
  4.8, 4.2,
 -0.70710678118654757, -0.70710678118654757, // 13
  6.1, 4.8,
  6.2, 4.7,
  6.3, 5.6, // 11x
  6.4, 5.5,
  6.5, 5.4,
  6.6, 5.3,
  6.7, 5.1,
  6.8, 5.2,
  6.9, 5.3,
  7.0, 5.9,
  7.7, 5.8,
  7.6, 5.7,
  0.0,-1.0, // 14
  7.5, 5.5,
  7.4, 6.4, // 11y
  7.3, 6.2,
  7.2, 6.0,
  7.2, 6.9,
  7.1, 6.8,
  7.0, 6.7,
  7.2, 6.6,
  7.3, 6.5,
  7.4, 6.3,
  7.5, 6.2,
 -1.0, 0.0, // 14
  8.6, 6.2,
  8.7, 6.3, // 12x
  8.8, 6.3,
  8.9, 7.4,
  8.0, 7.5,
  8.6, 7.6,
  8.5, 7.7,
  8.4, 7.8,
  8.0, 7.9,
  8.3, 7.1,
  8.0, 7.2,
  8.2, 7.3,
 +1.0, 0.0, // 15
  7.2, 8.6, // 12y
  6.6, 8.5,
  6.7, 8.4,
  6.9, 8.2,
  6.5, 8.3,
  6.4, 8.5,
  6.3, 8.6,
  6.5, 8.8,
  4.7, 8.7,
  4.9, 8.5,
  7.5, 8.3,
  6.0,-1.0, // 15
  3.2, 8.3, // 13x
 -0.70710678118654757, +0.70710678118654757, // 5
  5.4, 9.3,
  5.6, 9.7,
  3.7, 9.0,
  5.9, 9.9,
 +0.70710678118654757, -0.70710678118654757, // 10
  4.4, 4.8,
  4.6, 4.7,
  4.8, 4.6,
  4.9, 4.4,
  4.0, 4.2,
  4.2, 4.3, // 13y
 +0.70710678118654757, +0.70710678118654757, // 5
  7.5, 3.4,
  6.7, 3.5,
  6.4, 3.6,
  4.6, 3.9,
 -0.70710678118654757, -0.70710678118654757, // 10
  8.9, 2.8,
  7.6, 2.7,
  6.4, 2.6,
  5.3, 2.5,
  3.8, 2.3,
  4.5, 2.2, // 14x
  8.5, 2.4,
  0.0,+1.0, // 6
  7.4, 3.6,
  6.6, 3.5,
  4.7, 3.4,
  3.8, 3.5,
  0.0,-1.0, // 11
  5.9, 3.7,
  8.7, 4.6,
  7.6, 4.5,
  6.5, 4.4,
  5.5, 4.3, // 14y
  4.3, 4.8,
 +1.0, 0.0, // 6
  4.3, 4.7,
  6.5, 4.6,
  9.6, 4.5,
  8.7, 4.3,
 -1.0, 0.0, // 11
  7.9, 4.5,
  6.7, 5.3,
  5.6, 5.8,
  3.5, 5.7,
  2.4, 5.6, // 15x
  3.3, 5.5,
  4.2, 5.3,
  5.4, 5.6,
 -1.0, 6.0, // 8
  0.8, 6.6,
  9.8, 6.5,
  8.5, 6.5,
 +1.0, 0.0, // 12
  7.5, 7.3,
  6.4, 7.6,
  5.2, 7.8,
  4.5, 7.7, // 15y
  3.7, 8.6,
  3.8, 8.5,
  2.9, 8.3,
  0.0,+1.0, // 8
  2.0, 9.9,
  2.7, 9.8,
  1.6, 9.7,
  0.0,-1.0, // 12
  1.5, 5.5,
  1.2, 5.4,
  1.1, 5.3,
};


// ----------------------------------------------------------------------
// Computed values
// ----------------------------------------------------------------------

const double pylith::faults::CohesiveDynDataTri3d::_orientation[] = {
  +0.70710678118654757, -0.70710678118654757,  
  -0.70710678118654757, -0.70710678118654757,
  0.0, -1.0,  -1.0,  0.0,
 +1.0,  0.0,   0.0, -1.0
};

const double pylith::faults::CohesiveDynDataTri3d::_area[] = {
  2.0,
  1.0,
  1.0,
};

const int pylith::faults::CohesiveDynDataTri3d::_numConstraintVert = 3;
const int pylith::faults::CohesiveDynDataTri3d::_constraintVertices[] = {
  13, 14, 15
};

const double pylith::faults::CohesiveDynDataTri3d::_forcesInitial[] = {
  3.15*1.4142135623730951, 1.00*1.41421356237309,
  2.05, -1.05,
  1.10,  2.10,
};


// ----------------------------------------------------------------------
// Stick case
// ----------------------------------------------------------------------
// Input
const double pylith::faults::CohesiveDynDataTri3d::_fieldIncrStick[] = {
  1.1, 29.1,
  1.2, 29.2, // 5
  1.3, 29.3, // 6
  1.4, 29.4,
  1.5, 29.5, // 8
  1.6, 29.6,
  1.7, 29.7, // 10
  1.9, 29.9, // 11
  2.1, 29.1, // 12
  1.8, -29.8, // 13
  1.0, -29.0, // 14
  2.2, -29.2, // 15
};

// No change in fieldIncr
// Zero slip

// ----------------------------------------------------------------------
// Slip case
// ----------------------------------------------------------------------
// Input
const double pylith::faults::CohesiveDynDataTri3d::_fieldIncrSlip[] = {
  9.1, 10.1,
  9.2, 10.2, // 5
  9.3, 10.3, // 6
  9.4, 10.4,
  9.5, 10.5, // 8
  9.6, 10.6,
  9.7, 10.7, // 10
  9.9, 10.9, // 11
  9.1, 10.1, // 12
  9.8, -10.8, // 13
  9.0, -10.0, // 14
  9.2, -10.2, // 15
};

// Output
// TODO Update
const double pylith::faults::CohesiveDynDataTri3d::_fieldIncrSlipE[] = {
  9.100000000000,  10.100000000000,
  4.178047263424,  15.221952736576,
  9.300000000000,  10.050098043990,
  9.400000000000,  10.400000000000,
  9.655679817715,  10.500000000000,
  9.600000000000,  10.600000000000,
 14.721952736576,   5.678047263424,
  9.900000000000,  11.149901956010,
  8.944320182285,  10.100000000000,
 -5.600000000000, -10.800000000000,
 -4.800000000000, -10.000000000000,
 -6.600000000000, -10.200000000000,
};

const double pylith::faults::CohesiveDynDataTri3d::_slipSlipE[] = {
  14.204227339325,   0.0,
  -0.499803912020,   0.0,
  -0.311359635429,   0.0,
};

// ----------------------------------------------------------------------
// Open case
// ----------------------------------------------------------------------
// Input
const double pylith::faults::CohesiveDynDataTri3d::_fieldIncrOpen[] = {
  9.1, 10.1,
  9.2, 10.2, // 5
  9.3, 10.3, // 6
  9.4, 10.4,
  9.5, 10.5, // 8
  9.6, 10.6,
  9.7, 10.7, // 10
  9.9, 10.9, // 11
  9.1, 10.1, // 12
  9.8, 10.8, // 13
  9.0, 10.0, // 14
  9.2, 10.2, // 15
};

// Output
const double pylith::faults::CohesiveDynDataTri3d::_fieldIncrOpenE[] = {
   9.100000000000,  10.100000000000,
  17.946141606808,   1.453858393192,
   9.300000000000,  16.358051491522,
   9.400000000000,  10.400000000000,
   1.376240123890,  19.299429217321,
   9.600000000000,  10.600000000000,
   0.953858393192,  19.446141606808,
   9.900000000000,   4.841948508478,
  17.223759876110,   1.300570782679,
  -6.800000000000,  -8.800000000000,
  -6.000000000000,  -8.000000000000,
  -7.200000000000,  -9.200000000000,
};

const double pylith::faults::CohesiveDynDataTri3d::_slipOpenE[] = {
-24.737824157567,  0.0,
 12.116102983044,  0.0,
 16.247519752219,  17.598858434641,

};

// ----------------------------------------------------------------------
pylith::faults::CohesiveDynDataTri3d::CohesiveDynDataTri3d(void)
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
  forcesInitial = const_cast<double*>(_forcesInitial);

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

pylith::faults::CohesiveDynDataTri3d::~CohesiveDynDataTri3d(void)
{}


// End of file
