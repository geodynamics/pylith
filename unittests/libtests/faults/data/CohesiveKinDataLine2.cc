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

const char* pylith::faults::CohesiveKinDataLine2::_riseTimeFilename = 
  "data/line2_risetime.spatialdb";

// Don't expect these values to be used, so just use some values.
const double pylith::faults::CohesiveKinDataLine2::_fieldT[] = {
  7.1,
  7.2, // 3
  7.3,
  7.4, // 5
  7.5
};

const double pylith::faults::CohesiveKinDataLine2::_fieldIncr[] = {
  1.1,
  1.2, // 3
  1.3,
  1.4, // 5
  1.5
};

const double pylith::faults::CohesiveKinDataLine2::_jacobianLumped[] = {
  2.1,
  2.2, // 3
  2.3,
  2.4, // 5
  1.0
};

const int pylith::faults::CohesiveKinDataLine2::_numFaultVertices = 1;
const int pylith::faults::CohesiveKinDataLine2::_verticesFault[] = {
  1
};
const int pylith::faults::CohesiveKinDataLine2::_verticesLagrange[] = {
  6
};
const int pylith::faults::CohesiveKinDataLine2::_verticesPositive[] = {
  5
};
const int pylith::faults::CohesiveKinDataLine2::_verticesNegative[] = {
  3
};

const int pylith::faults::CohesiveKinDataLine2::_numCohesiveCells = 1;
const int pylith::faults::CohesiveKinDataLine2::_cellMappingFault[] = {
  0
};
const int pylith::faults::CohesiveKinDataLine2::_cellMappingCohesive[] = {
  7
};


const double pylith::faults::CohesiveKinDataLine2::_orientation[] = {
  1.0
};

const double pylith::faults::CohesiveKinDataLine2::_area[] = {
  1.0
};

const double pylith::faults::CohesiveKinDataLine2::_residualIncr[] = {
   0.0,
   7.5,
   0.0,
  -7.5,
  -0.2+1.89546413727,
};

const double pylith::faults::CohesiveKinDataLine2::_residual[] = {
   0.0,
   7.5, // 3
   0.0,
   -7.5, // 5
  -0.2+1.89546413727,
};

const double pylith::faults::CohesiveKinDataLine2::_jacobian[] = {
  0.0,  0.0,  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,  0.0, -1.0,
  0.0,  0.0,  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,  0.0, +1.0,
  0.0, -1.0,  0.0, +1.0,  0.0,
};

const double pylith::faults::CohesiveKinDataLine2::_fieldIncrAdjusted[] = {
  1.1,
  -3.09368089375, // 3
  1.3,
  5.33587415261, // 5
  -9.44609796626, // 6
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

pylith::faults::CohesiveKinDataLine2::~CohesiveKinDataLine2(void)
{}


// End of file
