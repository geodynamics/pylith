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
 *   Cells are 0-1, vertices are 2-4.
 *
 *   2 -------- 3 -------- 4
 *
 * After adding cohesive elements
 *   2 -------- 3-6-5 -------- 4
 */

#include "CohesiveKinDataLine2.hh"

const char* pylith::faults::CohesiveKinDataLine2::_meshFilename =
  "data/line2.mesh";

const int pylith::faults::CohesiveKinDataLine2::_spaceDim = 1;

const int pylith::faults::CohesiveKinDataLine2::_cellDim = 0;

const int pylith::faults::CohesiveKinDataLine2::_numBasis = 1;

const int pylith::faults::CohesiveKinDataLine2::_numQuadPts = 1;

const double pylith::faults::CohesiveKinDataLine2::_quadPts[] = {
  0.0,
};

const double pylith::faults::CohesiveKinDataLine2::_quadWts[] = {
  1.0,
};

const double pylith::faults::CohesiveKinDataLine2::_basis[] = {
  1.0,
};

const double pylith::faults::CohesiveKinDataLine2::_basisDeriv[] = {
  1.0
};

const double pylith::faults::CohesiveKinDataLine2::_verticesRef[] = {
  0.0
};

const int pylith::faults::CohesiveKinDataLine2::_id = 10;

const char* pylith::faults::CohesiveKinDataLine2::_label = "fault";

const char* pylith::faults::CohesiveKinDataLine2::_finalSlipFilename = 
  "data/line2_finalslip.spatialdb";

const char* pylith::faults::CohesiveKinDataLine2::_slipTimeFilename = 
  "data/line2_sliptime.spatialdb";

const char* pylith::faults::CohesiveKinDataLine2::_peakRateFilename = 
  "data/line2_peakrate.spatialdb";

const double pylith::faults::CohesiveKinDataLine2::_fieldTpdt[] = {
  0.0,
  0.0,
  0.0,
  0.0,
  0.0
};

const double pylith::faults::CohesiveKinDataLine2::_fieldT[] = {
  0.0,
  0.0,
  0.0,
  0.0,
  0.0
};

const int pylith::faults::CohesiveKinDataLine2::_numConstraintVert = 1;

const double pylith::faults::CohesiveKinDataLine2::_orientation[] = {
  1.0
};

const int pylith::faults::CohesiveKinDataLine2::_constraintVertices[] = {
  6
};

const double pylith::faults::CohesiveKinDataLine2::_valsResidual[] = {
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0
};

const double pylith::faults::CohesiveKinDataLine2::_valsSlip[] = {
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  1.2
};

const double pylith::faults::CohesiveKinDataLine2::_valsJacobian[] = {
  0.0, 0.0, 0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, 0.0, 0.0,
};

pylith::faults::CohesiveKinDataLine2::CohesiveKinDataLine2(void)
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
  fieldTpdt = const_cast<double*>(_fieldTpdt);
  fieldT = const_cast<double*>(_fieldT);
  orientation = const_cast<double*>(_orientation);
  constraintVertices = const_cast<int*>(_constraintVertices);
  valsResidual = const_cast<double*>(_valsResidual);
  valsSlip = const_cast<double*>(_valsSlip);
  valsJacobian = const_cast<double*>(_valsJacobian);
  numConstraintVert = _numConstraintVert;  
} // constructor

pylith::faults::CohesiveKinDataLine2::~CohesiveKinDataLine2(void)
{}


// End of file
