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

const int pylith::faults::CohesiveDynDataTri3::_numQuadPts = 1;

const double pylith::faults::CohesiveDynDataTri3::_quadPts[] = {
  0.0,
};

const double pylith::faults::CohesiveDynDataTri3::_quadWts[] = {
  2.0,
};

const double pylith::faults::CohesiveDynDataTri3::_basis[] = {
  0.5,
  0.5
};

const double pylith::faults::CohesiveDynDataTri3::_basisDeriv[] = {
  -0.5,
   0.5
};

const double pylith::faults::CohesiveDynDataTri3::_verticesRef[] = {
  -1.0, 1.0
};

const int pylith::faults::CohesiveDynDataTri3::_id = 10;

const char* pylith::faults::CohesiveDynDataTri3::_label = "fault";

const char* pylith::faults::CohesiveDynDataTri3::_initialTractFilename = 
  "data/tri3_initialtract.spatialdb";

const double pylith::faults::CohesiveDynDataTri3::_fieldT[] = {
  8.1, 9.1,
  8.2, 9.2, // 3
  8.3, 9.3, // 4
  8.4, 9.4,
  8.5, 9.5, // 6
  8.7, 9.7, // 7
  8.6, 9.6, // 8
  8.8, 9.8, // 9
};

// :TODO: Make sensible values for Jacobian for DOF on positive and
// negative sides of the fault. Add semi-random values for other DOF.
const double pylith::faults::CohesiveDynDataTri3::_jacobian[] = {
  1.0, 1.1, // 2x (values for row associated with vertex with label 2, x DOF)
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
  4.1, 5.1, // 3x (values for row associated with vertex label 3, x DOF)
  1.0, 5.2,
  4.2, 5.3,
  4.3, 5.4,
  4.4, 5.5,
  4.5, 5.6,
  4.6,+1.0, // 8 (column for DOF associated with vertex label 8)
  4.7, 5.7,
  6.1, 7.1, // 3y
  6.2, 1.0,
  6.3, 7.2,
  6.4, 7.3,
  6.5, 7.4,
  6.6, 7.5,
 +1.0, 7.6, //  8
  6.7, 7.7,
  8.1, 9.1, // 4x
  8.2, 9.2,
  1.0, 9.3,
  8.3, 9.4,
  8.4, 9.5,
  8.5, 9.6,
  8.6, 9.7,
  8.7,+1.0, //  9
  10.1, 11.1, // 4y
  10.2, 11.2,
  10.3, 1.0,
  10.4, 11.3,
  10.5, 11.4,
  10.6, 11.5,
  10.7, 11.6,
 +1.0, 11.7, //  9
  12.1, 13.1, // 5x
  12.2, 13.2,
  12.3, 13.3,
  1.0, 13.4,
  12.4, 13.5,
  12.5, 13.6,
  12.6, 13.7,
  12.7, 13.8,
  14.1, 15.1, // 5y
  14.2, 15.2,
  14.3, 15.3,
  14.4, 1.0,
  14.5, 15.4,
  14.6, 15.5,
  14.7, 15.6,
  14.8, 15.7,
  16.1, 17.1, // 6x
  16.2, 17.2,
  16.3, 17.3,
  16.4, 17.4,
  1.0, 17.5,
  16.5, 17.6,
  16.6,-1.0, //  8
  16.7, 17.7,
  18.1, 19.1, // 6y
  18.2, 19.2,
  18.3, 19.3,
  18.4, 19.4,
  18.5, 1.0,
  18.6, 19.5,
 -1.0, 19.6, //  8
  18.7, 19.7,
  20.1, 21.1, // 7x
  20.2, 21.2,
  20.3, 21.3,
  20.4, 21.4,
  20.5, 21.5,
  1.0, 21.6,
  20.6, 21.7,
  20.7,-1.0, //  9
  22.1, 23.1, // 7y
  22.2, 23.2,
  22.3, 23.3,
  22.4, 23.4,
  22.5, 23.5,
  22.6, 1.0,
  22.7, 23.6,
 -1.0, 23.7, //  9

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

const double pylith::faults::CohesiveDynDataTri3::_orientation[] = {
  0.0, -1.0,  -1.0, 0.0,
  0.0, -1.0,  -1.0, 0.0
};

const double pylith::faults::CohesiveDynDataTri3::_area[] = {
  1.0,
  1.0,
};

const double pylith::faults::CohesiveDynDataTri3::_forcesInitial[] = {
  2.05, -1.05,
  2.05, -1.05,
};


const int pylith::faults::CohesiveDynDataTri3::_numConstraintVert = 2;
const int pylith::faults::CohesiveDynDataTri3::_constraintVertices[] = {
  8, 9
};

// ----------------------------------------------------------------------
// Stick case
// ----------------------------------------------------------------------
// Input
const double pylith::faults::CohesiveDynDataTri3::_fieldIncrStick[] = {
  1.1, 29.1,
  1.2, 29.2, // 3
  1.3, 29.3, // 4
  1.4, 29.4,
  1.5, 29.5, // 6
  1.7, 29.7, // 7
  1.6, -29.6, // 8
  1.8, -29.8, // 9
};

// No change in fieldIncr
// Zero slip

// ----------------------------------------------------------------------
// Slip case
// ----------------------------------------------------------------------
// Input
const double pylith::faults::CohesiveDynDataTri3::_fieldIncrSlip[] = {
  9.1, 10.1,
  9.2, 10.2, // 3
  9.3, 10.3, // 4
  9.4, 10.4,
  9.5, 10.5, // 6
  9.7, 10.7, // 7
  9.6, -10.6, // 8
  9.8, -10.8, // 9
};

// Output
const double pylith::faults::CohesiveDynDataTri3::_fieldIncrSlipE[] = {
  9.1,  10.100000000000,
  9.2,   9.304186565683,
  9.3,  10.066643380561,
  9.4,  10.400000000000,
  9.5,  11.395813434317,
  9.7,  10.933356619439,
 -8.0, -10.600000000000,
 -8.2, -10.800000000000,
};

const double pylith::faults::CohesiveDynDataTri3::_slipSlipE[] = {
  -1.791626868633122,                   0.0,
  -0.466713238878134,                   0.0,
};

// ----------------------------------------------------------------------
// Open case
// ----------------------------------------------------------------------
// Input
const double pylith::faults::CohesiveDynDataTri3::_fieldIncrOpen[] = {
  9.1, 10.1,
  9.2, 10.2, // 3
  9.3, 10.3, // 4
  9.4, 10.4,
  9.5, 10.5, // 6
  9.7, 10.7, // 7
  9.6, 10.6, // 8
  9.8, 10.8, // 9
};

// Output
const double pylith::faults::CohesiveDynDataTri3::_fieldIncrOpenE[] = {
 9.100000000000,  10.100000000000,
 9.200000000000,  10.851688943849,
10.191208845155,  11.487449638489,
 9.400000000000,  10.400000000000,
 9.500000000000,   9.848311056151,
 8.808791154845,   9.512550361511,
-8.600000000000,  -9.600000000000,
-8.800000000000,  -9.800000000000,
};

const double pylith::faults::CohesiveDynDataTri3::_slipOpenE[] = {
  1.303377887698464,  0.0,
  2.374899276978848,  1.782417690310881,
};

// ----------------------------------------------------------------------
pylith::faults::CohesiveDynDataTri3::CohesiveDynDataTri3(void)
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

pylith::faults::CohesiveDynDataTri3::~CohesiveDynDataTri3(void)
{}


// End of file
