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

/* Original mesh using Sieve labels.
 *
 * Cells are 0-1, vertices are 2-5.
 *
 *              3
 *             /|\
 *            / | \
 *           /  |  \
 *          /   |   \
 *         2    |    5
 *          \   |   /
 *           \  |  /
 *            \ | /
 *             \|/
 *              4
 *
 *
 * After adding cohesive elements using Sieve labels.
 *
 * Cells are 0-1, 8, vertices are 2-7.
 *
 *              6 -8- 3
 *             /|     |\
 *            / |     | \
 *           /  |     |  \
 *          /   |     |   \
 *         2    |     |    5
 *          \   |     |   /
 *           \  |     |  /
 *            \ |     | /
 *             \|     |/
 *              7 -9- 4
 */

#include "CohesiveDynDataTri3.hh"

const char* pylith::faults::CohesiveDynDataTri3::_meshFilename =
  "data/tri3.mesh";

const int pylith::faults::CohesiveDynDataTri3::_spaceDim = 2;

const int pylith::faults::CohesiveDynDataTri3::_cellDim = 1;

const int pylith::faults::CohesiveDynDataTri3::_numBasis = 2;

const int pylith::faults::CohesiveDynDataTri3::_numQuadPts = 2;

const PylithScalar pylith::faults::CohesiveDynDataTri3::_quadPts[] = {
  -1.0, 1.0,
};

const PylithScalar pylith::faults::CohesiveDynDataTri3::_quadWts[] = {
  1.0, 1.0
};

const PylithScalar pylith::faults::CohesiveDynDataTri3::_basis[] = {
  1.0, 0.0,
  0.0, 1.0,
};

const PylithScalar pylith::faults::CohesiveDynDataTri3::_basisDeriv[] = {
  -0.5, 0.5,
  -0.5, 0.5,
};

const PylithScalar pylith::faults::CohesiveDynDataTri3::_verticesRef[] = {
  -1.0, 1.0
};

const int pylith::faults::CohesiveDynDataTri3::_id = 10;

const char* pylith::faults::CohesiveDynDataTri3::_label = "fault";

const char* pylith::faults::CohesiveDynDataTri3::_initialTractFilename = 
  "data/tri3_initialtract.spatialdb";

const PylithScalar pylith::faults::CohesiveDynDataTri3::_fieldT[] = {
  8.1, 9.1,
  8.2, 9.2, // 3
  8.3, 9.3, // 4
  8.4, 9.4,
  8.2, 9.2, // 6
  8.3, 9.3, // 7
  8.6, 9.6, // 8
  8.8, 9.8, // 9
};

const PylithScalar pylith::faults::CohesiveDynDataTri3::_jacobian[] = {
  1.0, 1.1, // 2x
  0.1, 1.2,
  0.2, 1.3,
  0.3, 1.4,
  0.4, 1.5,
  0.5, 1.6,
  0.6, 1.7,
  0.7, 1.8,
  2.1, 1.0, // 2y
  2.2, 3.1,
  2.3, 3.2,
  2.4, 3.3,
  2.5, 3.4,
  2.6, 3.5,
  2.7, 3.6,
  2.8, 3.7,
  4.1, 5.1, // 3x
 +4.0,-1.2, // 3
 -2.2,-2.3, // 4
  4.3, 5.4,
  4.4, 5.5,
  4.5, 5.6,
  4.6,+1.0, 
  4.7, 5.7,
  6.1, 7.1, // 3y
 -1.2,+5.0, // 3
 -1.3,-3.2, // 4
  6.4, 7.3,
  6.5, 7.4,
  6.6, 7.5,
  1.0, 7.6,
  6.7, 7.7,
  8.1, 9.1, // 4x
 -2.2,-1.3, // 3
 +4.1,-4.3, // 4
  8.3, 9.4,
  8.4, 9.5,
  8.5, 9.6,
  8.6, 9.7,
  8.7,+1.0,
  1.1, 1.1, // 4y
 -2.3,-3.2, // 3
 -4.3,+5.1, // 4
  1.4, 1.3,
  1.5, 1.4,
  1.6, 1.5,
  1.7, 1.6,
  1.0, 1.7,
  2.1, 3.1, // 5x
  2.2, 3.2,
  2.3, 3.3,
  1.0, 3.4,
  2.4, 3.5,
  2.5, 3.6,
  2.6, 3.7,
  2.7, 3.8,
  4.1, 5.1, // 5y
  4.2, 5.2,
  4.3, 5.3,
  4.4, 1.0,
  4.5, 5.4,
  4.6, 5.5,
  4.7, 5.6,
  4.8, 5.7,
  6.1, 7.1, // 6x
  6.2, 7.2,
  6.3, 7.3,
  6.4, 7.4,
 +5.0,-1.2, // 6
 -2.2,-2.3, // 7
  6.6, 1.0,
  6.7, 7.7,
  8.1, 9.1, // 6y
  8.2, 9.2,
  8.3, 9.3,
  8.4, 9.4,
 -1.2,+4.0, // 6
 -1.3,-3.2, // 7
  1.0, 9.6,
  8.7, 9.7,
  0.1, 1.1, // 7x
  0.2, 1.2,
  0.3, 1.3,
  0.4, 1.4,
 -2.2,-1.3, // 6
 +5.1,-4.3, // 7
  0.6, 1.7,
  0.7, 1.0,
  2.1, 3.1, // 7y
  2.2, 3.2,
  2.3, 3.3,
  2.4, 3.4,
 -2.3,-3.2, // 6
 -4.3,+4.1, // 7
  2.7, 3.6,
  1.0, 3.7, //  9

  24.1, 25.1, // 8x (rows associated with Lagrange multiplier vertex label 8)
  24.2,+1.0, //  3
  24.3, 25.2,
  24.4, 25.3,
  24.5,-1.0, //  6
  24.6, 25.4,
  24.7, 25.5,
  24.8, 25.6,
  26.1, 27.1, // 8y
 +1.0, 27.2, //  3
  26.2, 27.3,
  26.3, 27.4,
 -1.0, 27.5, //  6
  26.4, 27.6,
  26.5, 27.7,
  26.6, 27.8,

  29.1, 30.1, // 9x
  29.2, 30.2,
  29.3,+1.0, //  4
  29.4, 30.3,
  29.5, 30.4,
  29.6,-1.0, //  7
  29.7, 30.5,
  29.8, 30.6,
  31.1, 32.1, // 9y
  31.2, 32.2,
 +1.0, 32.3, //  4
  31.3, 32.4,
  31.4, 32.5,
 -1.0, 32.6, //  7
  31.5, 32.7,
  31.6, 32.8,
};

// ----------------------------------------------------------------------
// Computed values
// ----------------------------------------------------------------------

const PylithScalar pylith::faults::CohesiveDynDataTri3::_orientation[] = {
  0.0, -1.0,  -1.0, 0.0,
  0.0, -1.0,  -1.0, 0.0
};

const PylithScalar pylith::faults::CohesiveDynDataTri3::_area[] = {
  1.0,
  1.0,
};

const PylithScalar pylith::faults::CohesiveDynDataTri3::_initialTractions[] = {
  2.0, -1.0,
  2.1, -1.1,
};


const int pylith::faults::CohesiveDynDataTri3::_numConstraintVert = 2;
const int pylith::faults::CohesiveDynDataTri3::_constraintVertices[] = {
  8, 9
};

// ----------------------------------------------------------------------
// Stick case
// ----------------------------------------------------------------------
// Input
const PylithScalar pylith::faults::CohesiveDynDataTri3::_fieldIncrStick[] = {
  1.1, 2.1,
  1.2, 2.2, // 3
  1.3, 2.3, // 4
  1.4, 2.4,
  1.2, 2.2, // 6
  1.3, 2.3, // 7
  21.6, 2.6, // 8
  21.8, 2.8, // 9
};

// No slip
const PylithScalar pylith::faults::CohesiveDynDataTri3::_slipStickE[] = {
 0.0,  0.0,
 0.0,  0.0,
};

// ----------------------------------------------------------------------
// Slip case
// ----------------------------------------------------------------------
// Input
const PylithScalar pylith::faults::CohesiveDynDataTri3::_fieldIncrSlip[] = {
  9.1, 7.1,
  9.2, 7.2, // 3
  9.3, 7.3, // 4
  9.4, 7.4,
  9.2, 7.2, // 6
  9.3, 7.3, // 7
  1.6, 2.6, // 8
  1.8, 2.8, // 9
};

// Output
const PylithScalar pylith::faults::CohesiveDynDataTri3::_fieldIncrSlipE[] = {
   9.100000000000,   7.100000000000,
   9.200000000000,   7.390546440275,
   9.300000000000,   8.283993111958,
   9.400000000000,   7.400000000000,
   9.200000000000,   7.009453559725,
   9.300000000000,   6.316006888042,
   1.600000000000,  -3.480000000000,
   1.800000000000,  -3.440000000000,
};

const PylithScalar pylith::faults::CohesiveDynDataTri3::_slipSlipE[] = {
   0.381092880550,   0.000000000000,
   1.967986223916,   0.000000000000,
};

// ----------------------------------------------------------------------
// Open case
// ----------------------------------------------------------------------
// Input
const PylithScalar pylith::faults::CohesiveDynDataTri3::_fieldIncrOpen[] = {
  9.1, 7.1,
  9.2, 7.2, // 3
  9.3, 7.3, // 4
  9.4, 7.4,
  9.2, 7.2, // 6
  9.3, 7.3, // 7
  -10.6, 2.6, // 8
  -10.8, 2.8, // 9
};

// Output
const PylithScalar pylith::faults::CohesiveDynDataTri3::_fieldIncrOpenE[] = {
   9.100000000000,   7.100000000000,
  11.874337677762,   7.148344513504,
  12.357207850137,   8.757371967576,
   9.400000000000,   7.400000000000,
   6.525662322238,   7.251655486496,
   6.242792149863,   5.842628032424,
  -8.600000000000,  -9.600000000000,
  -8.800000000000,  -9.800000000000,
};

const PylithScalar pylith::faults::CohesiveDynDataTri3::_slipOpenE[] = {
  -0.103310972991,   5.348675355523,
   2.914743935152,   6.114415700274,
};

// ----------------------------------------------------------------------
pylith::faults::CohesiveDynDataTri3::CohesiveDynDataTri3(void)
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

pylith::faults::CohesiveDynDataTri3::~CohesiveDynDataTri3(void)
{}


// End of file
