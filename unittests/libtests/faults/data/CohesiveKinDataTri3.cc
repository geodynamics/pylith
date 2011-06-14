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
 * Cells are 0-1, 10, vertices are 2-7.
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

#include "CohesiveKinDataTri3.hh"

const char* pylith::faults::CohesiveKinDataTri3::_meshFilename =
  "data/tri3.mesh";

const int pylith::faults::CohesiveKinDataTri3::_spaceDim = 2;

const int pylith::faults::CohesiveKinDataTri3::_cellDim = 1;

const int pylith::faults::CohesiveKinDataTri3::_numBasis = 2;

const int pylith::faults::CohesiveKinDataTri3::_numQuadPts = 1;

const double pylith::faults::CohesiveKinDataTri3::_quadPts[] = {
  0.0,
};

const double pylith::faults::CohesiveKinDataTri3::_quadWts[] = {
  2.0,
};

const double pylith::faults::CohesiveKinDataTri3::_basis[] = {
  0.5,
  0.5
};

const double pylith::faults::CohesiveKinDataTri3::_basisDeriv[] = {
  -0.5,
   0.5
};

const double pylith::faults::CohesiveKinDataTri3::_verticesRef[] = {
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

const double pylith::faults::CohesiveKinDataTri3::_fieldT[] = {
  8.1, 9.1,
  8.2, 9.2, // 3
  8.3, 9.3, // 4
  8.4, 9.4,
  8.5, 9.5, // 6
  8.7, 9.7, // 7
  8.6, 9.6, // 8
  8.8, 9.8, // 9
};

const double pylith::faults::CohesiveKinDataTri3::_fieldIncr[] = {
  3.1, 4.1,
  3.2, 4.2, // 3
  3.3, 4.3, // 4
  3.4, 4.4,
  3.5, 4.5, // 6
  3.7, 4.7, // 7
  3.6, 4.6, // 8
  3.8, 4.8, // 9
};

const double pylith::faults::CohesiveKinDataTri3::_jacobianLumped[] = {
  1.1, 1.1,
  1.2, 1.2, // 3
  1.3, 1.3, // 4
  1.4, 1.4,
  1.5, 1.5, // 6
  1.7, 1.7, // 7
  1.6, 1.6, // 8
  1.8, 1.8, // 9
};

const int pylith::faults::CohesiveKinDataTri3::_numFaultVertices = 2;
const int pylith::faults::CohesiveKinDataTri3::_verticesFault[] = {
  1, 2
};
const int pylith::faults::CohesiveKinDataTri3::_verticesLagrange[] = {
  8, 9
};
const int pylith::faults::CohesiveKinDataTri3::_verticesNegative[] = {
  3, 4
};
const int pylith::faults::CohesiveKinDataTri3::_verticesPositive[] = {
  6, 7
};

const int pylith::faults::CohesiveKinDataTri3::_numCohesiveCells = 1;
const int pylith::faults::CohesiveKinDataTri3::_cellMappingFault[] = {
  0
};
const int pylith::faults::CohesiveKinDataTri3::_cellMappingCohesive[] = {
  10
};


const double pylith::faults::CohesiveKinDataTri3::_orientation[] = {
  0.0, -1.0,  -1.0, 0.0,
  0.0, -1.0,  -1.0, 0.0
};

const double pylith::faults::CohesiveKinDataTri3::_area[] = {
  1.0,
  1.0,
};

const double pylith::faults::CohesiveKinDataTri3::_residual[] = {
  0.0,  0.0,
 -9.6, -8.6, // 3
 -9.8, -8.8, // 4
  0.0,  0.0,
 +9.6, +8.6, // 6
 +9.8, +8.8, // 7
  0.3+1.89546413727,  0.3+0.08241148423, // 8
  0.4+1.77538035254,  0.4+0.14794836271, // 9
};

const double pylith::faults::CohesiveKinDataTri3::_residualIncr[] = {
  0.0,  0.0,
 -9.6, -8.6, // 3
 -9.8, -8.8, // 4
  0.0,  0.0,
 +9.6, +8.6, // 6
 +9.8, +8.8, // 7
  0.3+1.89546413727,  0.3+0.08241148423, // 8
  0.4+1.77538035254,  0.4+0.14794836271, // 9
};

const double pylith::faults::CohesiveKinDataTri3::_jacobian[] = {
  0.0, 0.0, // 2x
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0, // 2y
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0, // 3x
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0,+1.0, // 8
  0.0, 0.0,
  0.0, 0.0, // 3y
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
 +1.0, 0.0, //  8
  0.0, 0.0,
  0.0, 0.0, // 4x
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0,+1.0, //  9
  0.0, 0.0, // 4y
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
 +1.0, 0.0, //  9
  0.0, 0.0, // 5x
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
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
  0.0, 0.0, // 6x
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0,-1.0, //  8
  0.0, 0.0,
  0.0, 0.0, // 6y
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
 -1.0, 0.0, //  8
  0.0, 0.0,
  0.0, 0.0, // 7x
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0,-1.0, //  9
  0.0, 0.0, // 7y
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
 -1.0, 0.0, //  9

  0.0, 0.0, // 8x
  0.0,+1.0, //  3
  0.0, 0.0,
  0.0, 0.0,
  0.0,-1.0, //  6
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0, // 8y
 +1.0, 0.0, //  3
  0.0, 0.0,
  0.0, 0.0,
 -1.0, 0.0, //  6
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,

  0.0, 0.0, // 9x
  0.0, 0.0,
  0.0,+1.0, //  4
  0.0, 0.0,
  0.0, 0.0,
  0.0,-1.0, //  7
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0, // 9y
  0.0, 0.0,
 +1.0, 0.0, //  4
  0.0, 0.0,
  0.0, 0.0,
 -1.0, 0.0, //  7
  0.0, 0.0,
  0.0, 0.0,
};

const double pylith::faults::CohesiveKinDataTri3::_fieldIncrAdjusted[] = {
  3.1, 4.1,
  11.4124508246, 12.5863689652, // 3
  11.1489656107, 12.3019463023, // 4
  3.4, 4.4,
  -3.06996065966, -2.20909517212, // 6
  -2.30215017286, -1.41913540767, // 7
  -10.0636427582, -9.85494098949, // 8
  -10.402530193, -10.2036552939, // 9
};

pylith::faults::CohesiveKinDataTri3::CohesiveKinDataTri3(void)
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
  fieldIncr = const_cast<double*>(_fieldIncr);
  jacobianLumped = const_cast<double*>(_jacobianLumped);
  orientation = const_cast<double*>(_orientation);
  area = const_cast<double*>(_area);
  residualIncr = const_cast<double*>(_residualIncr);
  residual = const_cast<double*>(_residual);
  jacobian = const_cast<double*>(_jacobian);
  fieldIncrAdjusted = const_cast<double*>(_fieldIncrAdjusted);
  verticesFault = const_cast<int*>(_verticesFault);
  verticesLagrange = const_cast<int*>(_verticesLagrange);
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
