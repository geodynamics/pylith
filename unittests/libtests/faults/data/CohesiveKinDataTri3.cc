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
 * Cells are 0-1, 8, vertices are 2-7.
 *
 *              3 -7- 6
 *             /|     |\
 *            / |     | \
 *           /  |     |  \
 *          /   |     |   \
 *         2    |     |    5
 *          \   |     |   /
 *           \  |     |  /
 *            \ |     | /
 *             \|     |/
 *              4 -9- 8
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

const char* pylith::faults::CohesiveKinDataTri3::_peakRateFilename = 
  "data/tri3_peakrate.spatialdb";

const double pylith::faults::CohesiveKinDataTri3::_fieldTpdt[] = {
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
};

const double pylith::faults::CohesiveKinDataTri3::_fieldT[] = {
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
};

const int pylith::faults::CohesiveKinDataTri3::_numConstraintVert = 2;

const double pylith::faults::CohesiveKinDataTri3::_orientation[] = {
  0.0, 1.0,
  0.0, 1.0
};

const int pylith::faults::CohesiveKinDataTri3::_constraintVertices[] = {
  7, 9
};

const double pylith::faults::CohesiveKinDataTri3::_valsResidual[] = {
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
};

const double pylith::faults::CohesiveKinDataTri3::_valsSlip[] = {
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
};

const double pylith::faults::CohesiveKinDataTri3::_valsJacobian[] = {
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
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

pylith::faults::CohesiveKinDataTri3::~CohesiveKinDataTri3(void)
{}


// End of file
