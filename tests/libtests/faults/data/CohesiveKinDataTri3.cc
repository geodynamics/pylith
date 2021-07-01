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
 * After adding cohesive elements
 *
 * Cells are 0-1, 2, vertices are 3-10.
 *
 *              7 -9- 4
 *             /|     |\
 *            / |     | \
 *           /  |     |  \
 *          /   |     |   \
 *         3    |  2  |    6
 *          \   |     |   /
 *           \  |     |  /
 *            \ |     | /
 *             \|     |/
 *              8-10- 5
 */

#include "CohesiveKinDataTri3.hh"

const char* pylith::faults::CohesiveKinDataTri3::_meshFilename =
  "data/tri3.mesh";

const int pylith::faults::CohesiveKinDataTri3::_spaceDim = 2;

const int pylith::faults::CohesiveKinDataTri3::_cellDim = 1;

const int pylith::faults::CohesiveKinDataTri3::_numBasis = 2;

const int pylith::faults::CohesiveKinDataTri3::_numQuadPts = 2;

const PylithScalar pylith::faults::CohesiveKinDataTri3::_quadPts[] = {
  -1.0, 1.0,
};

const PylithScalar pylith::faults::CohesiveKinDataTri3::_quadWts[] = {
  1.0, 1.0
};

const PylithScalar pylith::faults::CohesiveKinDataTri3::_basis[] = {
  1.0, 0.0,
  0.0, 1.0,
};

const PylithScalar pylith::faults::CohesiveKinDataTri3::_basisDeriv[] = {
  -0.5, 0.5,
  -0.5, 0.5,
};

const PylithScalar pylith::faults::CohesiveKinDataTri3::_verticesRef[] = {
  -1.0, 1.0
};

const int pylith::faults::CohesiveKinDataTri3::_id = 10;

const char* pylith::faults::CohesiveKinDataTri3::_label = "fault";

const char* pylith::faults::CohesiveKinDataTri3::_finalSlipFilename = 
  "data/tri3_finalslip.spatialdb";

const char* pylith::faults::CohesiveKinDataTri3::_slipTimeFilename = 
  "data/tri3_sliptime.spatialdb";

const char* pylith::faults::CohesiveKinDataTri3::_riseTimeFilename = 
  "data/tri3_risetime.spatialdb";

const PylithScalar pylith::faults::CohesiveKinDataTri3::_fieldT[] = {
  8.1, 9.1,
  8.2, 9.2, // 4
  8.3, 9.3, // 5
  8.4, 9.4,
  8.5, 9.5, // 7
  8.7, 9.7, // 8
  8.6, 9.6, // 15
  8.8, 9.8, // 16
};

const PylithScalar pylith::faults::CohesiveKinDataTri3::_fieldIncr[] = {
  3.1, 4.1,
  3.2, 4.2, // 4
  3.3, 4.3, // 5
  3.4, 4.4,
  3.5, 4.5, // 7
  3.7, 4.7, // 8
  3.6, 4.6, // 15
  3.8, 4.8, // 16
};

const PylithScalar pylith::faults::CohesiveKinDataTri3::_jacobianLumped[] = {
  1.1, 1.1,
  1.2, 1.2, // 4
  1.3, 1.3, // 5
  1.4, 1.4,
  1.5, 1.5, // 7
  1.7, 1.7, // 8
  1.0, 1.0, // 15
  1.0, 1.0, // 16
};

const int pylith::faults::CohesiveKinDataTri3::_numFaultVertices = 2;
const int pylith::faults::CohesiveKinDataTri3::_verticesFault[] = {
  2, 3
};
const int pylith::faults::CohesiveKinDataTri3::_edgesLagrange[] = {
  15, 16
};
const int pylith::faults::CohesiveKinDataTri3::_verticesNegative[] = {
  4, 5
};
const int pylith::faults::CohesiveKinDataTri3::_verticesPositive[] = {
  7, 8
};

const int pylith::faults::CohesiveKinDataTri3::_numCohesiveCells = 1;
const int pylith::faults::CohesiveKinDataTri3::_cellMappingFault[] = {
  4
};
const int pylith::faults::CohesiveKinDataTri3::_cellMappingCohesive[] = {
  2
};


const PylithScalar pylith::faults::CohesiveKinDataTri3::_orientation[] = {
  0.0, +1.0,  +1.0, 0.0,
  0.0, +1.0,  +1.0, 0.0
};

const PylithScalar pylith::faults::CohesiveKinDataTri3::_area[] = {
  1.0,
  1.0,
};

const PylithScalar pylith::faults::CohesiveKinDataTri3::_residual[] = {
  0.0,  0.0,
 +8.6, +9.6, // 4
 +8.8, +9.8, // 5
  0.0,  0.0,
 -8.6, -9.6, // 7
 -8.8, -9.8, // 8
 -(8.5-8.2) + (0.08241148423),
 -(9.5-9.2) + (1.89546413727), // 15
 -(8.7-8.3) + (0.14794836271),
 -(9.7-9.3) + (1.77538035254), // 16
};

const PylithScalar pylith::faults::CohesiveKinDataTri3::_jacobian[] = {
  0.0, 0.0, // 3x
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0, // 3y
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0, // 4x
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
 -1.0, 0.0, // 15
  0.0, 0.0,
  0.0, 0.0, // 4y
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0,-1.0, // 15
  0.0, 0.0,
  0.0, 0.0, // 5x
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
 -1.0, 0.0, // 16
  0.0, 0.0, // 5y
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0,-1.0, // 16
  0.0, 0.0, // 6x
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0, // 6y
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0, // 7x
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
 +1.0, 0.0, // 15
  0.0, 0.0,
  0.0, 0.0, // 7y
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0,+1.0, // 15
  0.0, 0.0,
  0.0, 0.0, // 8x
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
 +1.0, 0.0, // 16
  0.0, 0.0, // 8y
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0,+1.0, // 16

  0.0, 0.0, // 15x
 -1.0, 0.0, // 4
  0.0, 0.0,
  0.0, 0.0,
 +1.0, 0.0, // 7
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0, // 15y
  0.0,-1.0, // 4
  0.0, 0.0,
  0.0, 0.0,
  0.0,+1.0, // 7
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,

  0.0, 0.0, // 16x
  0.0, 0.0,
 -1.0, 0.0, // 5
  0.0, 0.0,
  0.0, 0.0,
 +1.0, 0.0, // 8
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0, // 16y
  0.0, 0.0,
  0.0,-1.0, // 5
  0.0, 0.0,
  0.0, 0.0,
  0.0,+1.0, // 8
  0.0, 0.0,
  0.0, 0.0,
};

const PylithScalar pylith::faults::CohesiveKinDataTri3::_fieldIncrAdjusted[] = {
  3.1, 4.1,
  3.4875491754, 3.48029770152, // 4
  3.6694959278, 3.74728446689, // 5
  3.4, 4.4,
  3.26996065966, 5.07576183879, // 7
  3.41744429051, 5.12266481943, // 8
  0.345059010513, -0.86364275818, // 15
  0.480344706137, -0.718530193038, // 16
};

pylith::faults::CohesiveKinDataTri3::CohesiveKinDataTri3(void)
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
  fieldIncrAdjusted = const_cast<PylithScalar*>(_fieldIncrAdjusted);
  verticesFault = const_cast<int*>(_verticesFault);
  edgesLagrange = const_cast<int*>(_edgesLagrange);
  verticesNegative = const_cast<int*>(_verticesNegative);
  verticesPositive = const_cast<int*>(_verticesPositive);
  numFaultVertices = _numFaultVertices;  
  cellMappingFault = const_cast<int*>(_cellMappingFault);
  cellMappingCohesive = const_cast<int*>(_cellMappingCohesive);
  numCohesiveCells = _numCohesiveCells;  
} // constructor

pylith::faults::CohesiveKinDataTri3::~CohesiveKinDataTri3(void)
{}


// End of file
