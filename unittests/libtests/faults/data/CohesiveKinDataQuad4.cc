// -*- C++ -*-
//
// ======================================================================
//
//                           Brad T. Aagaard
//                        U.S. Geological Survey
//
// {LicenseText}
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
 *       3 -------- 5 -11- 10 -------- 7
 *       |          |       |          |
 *       |          |       |          |
 *       |          |       |          |
 *       |          |       |          |
 *       |          |       |          |
 *       |          |       |          |
 *       |          |       |          |
 *       |          |       |          |
 *       2 -------- 4 --9-- 8 -------- 6
 */

#include "CohesiveKinDataQuad4.hh"

const char* pylith::faults::CohesiveKinDataQuad4::_meshFilename =
  "data/quad4.mesh";

const int pylith::faults::CohesiveKinDataQuad4::_spaceDim = 2;

const int pylith::faults::CohesiveKinDataQuad4::_cellDim = 1;

const int pylith::faults::CohesiveKinDataQuad4::_numBasis = 2;

const int pylith::faults::CohesiveKinDataQuad4::_numQuadPts = 1;

const double pylith::faults::CohesiveKinDataQuad4::_quadPts[] = {
  0.0,
};

const double pylith::faults::CohesiveKinDataQuad4::_quadWts[] = {
  2.0,
};

const double pylith::faults::CohesiveKinDataQuad4::_basis[] = {
  0.5,
  0.5
};

const double pylith::faults::CohesiveKinDataQuad4::_basisDeriv[] = {
  -0.5,
   0.5
};

const double pylith::faults::CohesiveKinDataQuad4::_verticesRef[] = {
  -1.0, 1.0
};

const int pylith::faults::CohesiveKinDataQuad4::_id = 10;

const char* pylith::faults::CohesiveKinDataQuad4::_label = "fault";

const char* pylith::faults::CohesiveKinDataQuad4::_finalSlipFilename = 
  "data/quad4_finalslip.spatialdb";

const char* pylith::faults::CohesiveKinDataQuad4::_slipTimeFilename = 
  "data/quad4_sliptime.spatialdb";

const char* pylith::faults::CohesiveKinDataQuad4::_peakRateFilename = 
  "data/quad4_peakrate.spatialdb";

const double pylith::faults::CohesiveKinDataQuad4::_fieldT[] = {
  8.1, 9.1,
  8.2, 9.2,
  8.3, 9.3,
  8.4, 9.4,
  8.5, 9.5,
  8.6, 9.6,
  8.7, 9.7,
  8.8, 9.8,
  8.9, 9.9,
  8.0, 9.0,
};


const int pylith::faults::CohesiveKinDataQuad4::_numConstraintVert = 2;

const double pylith::faults::CohesiveKinDataQuad4::_orientation[] = {
  0.0,  1.0,  +1.0, 0.0,
  0.0,  1.0,  +1.0, 0.0
};

const int pylith::faults::CohesiveKinDataQuad4::_constraintVertices[] = {
  9, 11
};

const int pylith::faults::CohesiveKinDataQuad4::_constraintCells[] = {
  12, 12
};

const double pylith::faults::CohesiveKinDataQuad4::_valsResidual[] = {
  0.0,  0.0,
  0.0,  0.0,
  0.0,  0.0,
  0.0,  0.0,
  0.0,  0.0,
  0.0,  0.0,
  0.0,  0.0,
  0.989535448086, 0.0824612873405, // 9
  0.0,  0.0,
  1.05057813143, 0.0456773100622, // 11
};

const double pylith::faults::CohesiveKinDataQuad4::_valsJacobian[] = {
  0.0, 0.0, // 2x
  0.0, 0.0,
  0.0, 0.0,
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
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0, // 3x
  0.0, 0.0,
  0.0, 0.0,
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
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0, // 4x
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0,-1.0, // 9
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0, // 4y
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
 -1.0, 0.0, //  9
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
  0.0,-1.0, //  11
  0.0, 0.0, // 5y
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
 -1.0, 0.0, //  11
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
  0.0, 0.0, // 8x
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0, 
  0.0, 0.0,
  0.0, 0.0, 
  0.0, 0.0,
  0.0,+1.0, //  9
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0, // 8y
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
 +1.0, 0.0, //  9
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0, // 9x
  0.0, 0.0,
  0.0,-1.0, //  4
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0,+1.0, //  8
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0, // 9y
  0.0, 0.0,
 -1.0, 0.0, //  4
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
 +1.0, 0.0, //  8
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
  0.0,+1.0, // 11
  0.0, 0.0, // 10y
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
 +1.0, 0.0, // 11
  0.0, 0.0, // 11x
  0.0, 0.0,
  0.0, 0.0,
  0.0,-1.0, //  5
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0,+1.0, //  10
  0.0, 0.0,
  0.0, 0.0, // 11y
  0.0, 0.0,
  0.0, 0.0,
 -1.0, 0.0, //  5
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
 +1.0, 0.0, //  10
  0.0, 0.0,
};

pylith::faults::CohesiveKinDataQuad4::CohesiveKinDataQuad4(void)
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
  peakRateFilename = const_cast<char*>(_peakRateFilename);
  fieldT = const_cast<double*>(_fieldT);
  orientation = const_cast<double*>(_orientation);
  constraintVertices = const_cast<int*>(_constraintVertices);
  constraintCells = const_cast<int*>(_constraintCells);
  valsResidual = const_cast<double*>(_valsResidual);
  valsJacobian = const_cast<double*>(_valsJacobian);
  numConstraintVert = _numConstraintVert;  
} // constructor

pylith::faults::CohesiveKinDataQuad4::~CohesiveKinDataQuad4(void)
{}


// End of file
