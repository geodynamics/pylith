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


#include "CohesiveKinDataTri3d.hh"

const char* pylith::faults::CohesiveKinDataTri3d::_meshFilename =
  "data/tri3d.mesh";

const int pylith::faults::CohesiveKinDataTri3d::_spaceDim = 2;

const int pylith::faults::CohesiveKinDataTri3d::_cellDim = 1;

const int pylith::faults::CohesiveKinDataTri3d::_numBasis = 2;

const int pylith::faults::CohesiveKinDataTri3d::_numQuadPts = 2;

const PylithScalar pylith::faults::CohesiveKinDataTri3d::_quadPts[] = {
  -1.0, 1.0,
};

const PylithScalar pylith::faults::CohesiveKinDataTri3d::_quadWts[] = {
  1.0, 1.0
};

const PylithScalar pylith::faults::CohesiveKinDataTri3d::_basis[] = {
  1.0, 0.0,
  0.0, 1.0,
};

const PylithScalar pylith::faults::CohesiveKinDataTri3d::_basisDeriv[] = {
  -0.5, 0.5,
  -0.5, 0.5,
};

const PylithScalar pylith::faults::CohesiveKinDataTri3d::_verticesRef[] = {
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

const PylithScalar pylith::faults::CohesiveKinDataTri3d::_fieldT[] = {
  6.1, 8.1,
  6.2, 8.2, // 5
  6.3, 8.3, // 6
  6.4, 8.4,
  6.5, 8.5, // 8
  6.6, 8.6,
  6.7, 8.7, // 10
  6.9, 8.9, // 11
  7.1, 9.1, // 12
  6.8, 8.8, // 26
  6.0, 8.0, // 27
  7.2, 9.2, // 28
};

const PylithScalar pylith::faults::CohesiveKinDataTri3d::_fieldIncr[] = {
  3.1, 7.1,
  3.2, 7.2, // 5
  3.3, 7.3, // 6
  3.4, 7.4,
  3.5, 7.5, // 8
  3.6, 7.6,
  3.7, 7.7, // 10
  3.9, 7.9, // 11
  3.1, 7.1, // 12
  3.8, 7.8, // 26
  3.0, 7.0, // 27
  2.2, 5.2, // 28
};

const PylithScalar pylith::faults::CohesiveKinDataTri3d::_jacobianLumped[] = {
  6.1, 8.1,
  6.2, 8.2, // 5
  6.3, 8.3, // 6
  6.4, 8.4,
  6.5, 8.5, // 8
  6.6, 8.6,
  6.7, 8.7, // 10
  6.9, 8.9, // 11
  7.1, 9.1, // 12
  1.0, 1.0, // 26
  1.0, 1.0, // 27
  1.0, 1.0, // 28
};


const PylithScalar pylith::faults::CohesiveKinDataTri3d::_orientation[] = {
  -0.70710678118654757, +0.70710678118654757,  
  +0.70710678118654757, +0.70710678118654757,
  0.0, +1.0,  +1.0,  0.0,
 -1.0,  0.0,   0.0, +1.0
};

const PylithScalar pylith::faults::CohesiveKinDataTri3d::_area[] = {
  2.0,
  1.0,
  1.0,
};

const int pylith::faults::CohesiveKinDataTri3d::_numFaultVertices = 3;
const int pylith::faults::CohesiveKinDataTri3d::_verticesFault[] = {
   4,  5,  6
};
const int pylith::faults::CohesiveKinDataTri3d::_edgesLagrange[] = {
  26, 27, 28
};
const int pylith::faults::CohesiveKinDataTri3d::_verticesNegative[] = {
   7,  8, 10
};
const int pylith::faults::CohesiveKinDataTri3d::_verticesPositive[] = {
  12, 13, 14
};

const int pylith::faults::CohesiveKinDataTri3d::_numCohesiveCells = 2;
const int pylith::faults::CohesiveKinDataTri3d::_cellMappingFault[] = {
  7, 8
};
const int pylith::faults::CohesiveKinDataTri3d::_cellMappingCohesive[] = {
  4, 5
};

const PylithScalar pylith::faults::CohesiveKinDataTri3d::_residual[] = {
  0.0,  0.0,
 +2.0*6.8, +2.0*8.8, // 5
 +1.0*6.0, +1.0*8.0, // 6
  0.0,  0.0,
 +1.0*7.2, +1.0*9.2, // 8
  0.0,  0.0,
 -2.0*6.8, -2.0*8.8, // 10
 -1.0*6.0, -1.0*8.0, // 6
 -1.0*7.2, -1.0*9.2, // 8
  -2.0*(6.7-6.2 +0.70710678118654757*(1.89546413727-0.08241148423)),
  -2.0*(8.7-8.2 +0.70710678118654757*(-1.89546413727-0.08241148423)), // 26
  -1.0*(6.9-6.3 -0.14794836271),
  -1.0*(8.9-8.3 -1.77538035254), // 27
  -1.0*(7.1-6.5 +1.59887481971),
  -1.0*(9.1-8.5 -0.19186497837), // 28
};

const PylithScalar pylith::faults::CohesiveKinDataTri3d::_jacobian[] = {
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
 -2.0, 0.0, // 13
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
  0.0,-2.0, // 13
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
 -1.0, 0.0, // 14
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
  0.0,-1.0, // 14
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
  0.0,-1.0, // 15
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
 +2.0, 0.0, // 13
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
  0.0,+2.0, // 13
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
 +1.0, 0.0, // 14
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
  0.0,+1.0, // 14
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
  0.0,+1.0, // 15
  0.0, 0.0, // 13x
 -2.0, 0.0, // 5
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
 +2.0, 0.0, // 10
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0, // 13y
  0.0,-2.0, // 5
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0,+2.0, // 10
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0, // 14x
  0.0, 0.0,
 -1.0, 0.0, // 6
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
 +1.0, 0.0, // 11
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0, // 14y
  0.0, 0.0,
  0.0,-1.0, // 6
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0,+1.0, // 11
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
  0.0,-1.0, // 8
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0,+1.0, // 12
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
  quadPts = const_cast<PylithScalar*>(_quadPts);
  quadWts = const_cast<PylithScalar*>(_quadWts);
  basis = const_cast<PylithScalar*>(_basis);
  basisDeriv = const_cast<PylithScalar*>(_basisDeriv);
  verticesRef = const_cast<PylithScalar*>(_verticesRef);
  id = _id;
  label = const_cast<char*>(_label);
  finalSlipFilename = const_cast<char*>(_finalSlipFilename);
  slipTimeFilename = const_cast<char*>(_slipTimeFilename);
  riseTimeFilename = const_cast<char*>(_riseTimeFilename);
  fieldT = const_cast<PylithScalar*>(_fieldT);
  fieldIncr = const_cast<PylithScalar*>(_fieldIncr);
  jacobianLumped = const_cast<PylithScalar*>(_jacobianLumped);
  orientation = const_cast<PylithScalar*>(_orientation);
  area = const_cast<PylithScalar*>(_area);
  residual = const_cast<PylithScalar*>(_residual);
  jacobian = const_cast<PylithScalar*>(_jacobian);
  verticesFault = const_cast<int*>(_verticesFault);
  edgesLagrange = const_cast<int*>(_edgesLagrange);
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
