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


#include "CohesiveKinDataTri3d.hh"

const char* pylith::faults::CohesiveKinDataTri3d::_meshFilename =
  "data/tri3d.mesh";

const int pylith::faults::CohesiveKinDataTri3d::_spaceDim = 2;

const int pylith::faults::CohesiveKinDataTri3d::_cellDim = 1;

const int pylith::faults::CohesiveKinDataTri3d::_numBasis = 2;

const int pylith::faults::CohesiveKinDataTri3d::_numQuadPts = 1;

const double pylith::faults::CohesiveKinDataTri3d::_quadPts[] = {
  0.0,
};

const double pylith::faults::CohesiveKinDataTri3d::_quadWts[] = {
  2.0,
};

const double pylith::faults::CohesiveKinDataTri3d::_basis[] = {
  0.5,
  0.5
};

const double pylith::faults::CohesiveKinDataTri3d::_basisDeriv[] = {
  -0.5,
   0.5
};

const double pylith::faults::CohesiveKinDataTri3d::_verticesRef[] = {
  -1.0, 1.0
};

const int pylith::faults::CohesiveKinDataTri3d::_id = 10;

const char* pylith::faults::CohesiveKinDataTri3d::_label = "fault";

const char* pylith::faults::CohesiveKinDataTri3d::_finalSlipFilename = 
  "data/tri3d_finalslip.spatialdb";

const char* pylith::faults::CohesiveKinDataTri3d::_slipTimeFilename = 
  "data/tri3d_sliptime.spatialdb";

const char* pylith::faults::CohesiveKinDataTri3d::_riseTimeFilename = 
  "data/tri3d_risetime.spatialdb";

const double pylith::faults::CohesiveKinDataTri3d::_fieldT[] = {
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


const double pylith::faults::CohesiveKinDataTri3d::_orientation[] = {
  +0.70710678118654757, -0.70710678118654757,  
  -0.70710678118654757, -0.70710678118654757,
  0.0, -1.0,  -1.0,  0.0,
 +1.0,  0.0,   0.0, -1.0
};

const double pylith::faults::CohesiveKinDataTri3d::_area[] = {
  2.0,
  1.0,
  1.0,
};

const int pylith::faults::CohesiveKinDataTri3d::_numFaultVertices = 3;
const int pylith::faults::CohesiveKinDataTri3d::_verticesFault[] = {
  4, 2, 3
};
const int pylith::faults::CohesiveKinDataTri3d::_verticesLagrange[] = {
  15, 13, 14
};
const int pylith::faults::CohesiveKinDataTri3d::_verticesNegative[] = {
  8, 5, 6
};
const int pylith::faults::CohesiveKinDataTri3d::_verticesPositive[] = {
  12, 10, 11
};

const int pylith::faults::CohesiveKinDataTri3d::_numCohesiveCells = 2;
const int pylith::faults::CohesiveKinDataTri3d::_cellMappingFault[] = {
  0, 1
};
const int pylith::faults::CohesiveKinDataTri3d::_cellMappingCohesive[] = {
  16, 17
};

const double pylith::faults::CohesiveKinDataTri3d::_residual[] = {
  0.0,  0.0,
 -1.4142135623730949, -11.030865786510143, // 5
 -8.0,  -6.0, // 6
  0.0,  0.0,
 +7.2,  -9.2, // 8
  0.0,  0.0,
 +1.4142135623730949, +11.030865786510143, // 10
 +8.0, +6.0, // 11
 -7.2, +9.2, // 12
  0.0+1.89546413727, +0.70710678118654757+0.08241148423, // 13
  0.6+1.77538035254, 0.6+0.14794836271, // 14
 -0.6+1.59887481971,  0.6+0.19186497837, // 15
};

const double pylith::faults::CohesiveKinDataTri3d::_residualIncr[] = {
  0.0,  0.0,
 -1.4142135623730949, -11.030865786510143, // 5
 -8.0,  -6.0, // 6
  0.0,  0.0,
 +7.2,  -9.2, // 8
  0.0,  0.0,
 +1.4142135623730949, +11.030865786510143, // 10
 +8.0, +6.0, // 11
 -7.2, +9.2, // 12
  0.0+1.89546413727, +0.70710678118654757+0.08241148423, // 13
  0.6+1.77538035254, 0.6+0.14794836271, // 14
 -0.6+1.59887481971,  0.6+0.19186497837, // 15
};

const double pylith::faults::CohesiveKinDataTri3d::_jacobian[] = {
  0.0, 0.0, // 4x
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0, // 4y
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0, // 5x
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
 -0.70710678118654757, +0.70710678118654757, // 13
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0, // 5y
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
 +0.70710678118654757, +0.70710678118654757, // 13
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0, // 6x
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0,+1.0, // 14
  0.0, 0.0,
  0.0, 0.0, // 6y
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
 +1.0, 0.0, // 14
  0.0, 0.0,
  0.0, 0.0, // 7x
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0, // 7y
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0, // 8x
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
 -1.0, 0.0, // 15
  0.0, 0.0, // 8y
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0,+1.0, // 15
  0.0, 0.0, // 9x
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0, // 9y
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0, // 10x
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
 +0.70710678118654757, -0.70710678118654757, // 13
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0, // 10y
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
 -0.70710678118654757, -0.70710678118654757, // 13
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0, // 11x
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0,-1.0, // 14
  0.0, 0.0,
  0.0, 0.0, // 11y
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
 -1.0, 0.0, // 14
  0.0, 0.0,
  0.0, 0.0, // 12x
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
 +1.0, 0.0, // 15
  0.0, 0.0, // 12y
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0,-1.0, // 15
  0.0, 0.0, // 13x
 -0.70710678118654757, +0.70710678118654757, // 5
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
 +0.70710678118654757, -0.70710678118654757, // 10
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0, // 13y
 +0.70710678118654757, +0.70710678118654757, // 5
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
 -0.70710678118654757, -0.70710678118654757, // 10
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0, // 14x
  0.0, 0.0,
  0.0,+1.0, // 6
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0,-1.0, // 11
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0, // 14y
  0.0, 0.0,
 +1.0, 0.0, // 6
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
 -1.0, 0.0, // 11
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0, // 15x
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
 -1.0, 0.0, // 8
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
 +1.0, 0.0, // 12
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0, // 15y
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0,+1.0, // 8
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0,-1.0, // 12
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
};

pylith::faults::CohesiveKinDataTri3d::CohesiveKinDataTri3d(void)
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
  finalSlipFilename = const_cast<char*>(_finalSlipFilename);
  slipTimeFilename = const_cast<char*>(_slipTimeFilename);
  riseTimeFilename = const_cast<char*>(_riseTimeFilename);
  fieldT = const_cast<double*>(_fieldT);
  orientation = const_cast<double*>(_orientation);
  area = const_cast<double*>(_area);
  residual = const_cast<double*>(_residual);
  residualIncr = const_cast<double*>(_residualIncr);
  jacobian = const_cast<double*>(_jacobian);
  verticesFault = const_cast<int*>(_verticesFault);
  verticesLagrange = const_cast<int*>(_verticesLagrange);
  verticesNegative = const_cast<int*>(_verticesNegative);
  verticesPositive = const_cast<int*>(_verticesPositive);
  numFaultVertices = _numFaultVertices;  
  cellMappingFault = const_cast<int*>(_cellMappingFault);
  cellMappingCohesive = const_cast<int*>(_cellMappingCohesive);
  numCohesiveCells = _numCohesiveCells;  
} // constructor

pylith::faults::CohesiveKinDataTri3d::~CohesiveKinDataTri3d(void)
{}


// End of file
