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

#include "CohesiveKinSrcsDataLine2.hh"

const char* pylith::faults::CohesiveKinSrcsDataLine2::_meshFilename =
  "data/line2.mesh";

const int pylith::faults::CohesiveKinSrcsDataLine2::_spaceDim = 1;

const int pylith::faults::CohesiveKinSrcsDataLine2::_cellDim = 0;

const int pylith::faults::CohesiveKinSrcsDataLine2::_numBasis = 1;

const int pylith::faults::CohesiveKinSrcsDataLine2::_numQuadPts = 1;

const PylithScalar pylith::faults::CohesiveKinSrcsDataLine2::_quadPts[] = {
  0.0,
};

const PylithScalar pylith::faults::CohesiveKinSrcsDataLine2::_quadWts[] = {
  1.0,
};

const PylithScalar pylith::faults::CohesiveKinSrcsDataLine2::_basis[] = {
  1.0,
};

const PylithScalar pylith::faults::CohesiveKinSrcsDataLine2::_basisDeriv[] = {
  1.0
};

const PylithScalar pylith::faults::CohesiveKinSrcsDataLine2::_verticesRef[] = {
  0.0
};

const int pylith::faults::CohesiveKinSrcsDataLine2::_id = 10;

const char* pylith::faults::CohesiveKinSrcsDataLine2::_label = "fault";

const char* pylith::faults::CohesiveKinSrcsDataLine2::_finalSlipFilename = 
  "data/line2_finalslip.spatialdb";

const char* pylith::faults::CohesiveKinSrcsDataLine2::_slipTimeFilename = 
  "data/line2_sliptime.spatialdb";

const char* pylith::faults::CohesiveKinSrcsDataLine2::_riseTimeFilename = 
  "data/line2_risetime.spatialdb";

// Don't expect these values to be used, so just use some values.
const PylithScalar pylith::faults::CohesiveKinSrcsDataLine2::_fieldT[] = {
  7.1,
  7.2,
  7.3,
  7.4,
  7.5
};


const PylithScalar pylith::faults::CohesiveKinSrcsDataLine2::_orientation[] = {
  1.0
};

const PylithScalar pylith::faults::CohesiveKinSrcsDataLine2::_area[] = {
  1.0
};

const int pylith::faults::CohesiveKinSrcsDataLine2::_numFaultVertices = 1;
const int pylith::faults::CohesiveKinSrcsDataLine2::_verticesFault[] = {
  1
};
const int pylith::faults::CohesiveKinSrcsDataLine2::_verticesLagrange[] = {
  6
};
const int pylith::faults::CohesiveKinSrcsDataLine2::_verticesPositive[] = {
  5
};
const int pylith::faults::CohesiveKinSrcsDataLine2::_verticesNegative[] = {
  3
};

const int pylith::faults::CohesiveKinSrcsDataLine2::_numCohesiveCells = 1;
const int pylith::faults::CohesiveKinSrcsDataLine2::_cellMappingFault[] = {
  0
};
const int pylith::faults::CohesiveKinSrcsDataLine2::_cellMappingCohesive[] = {
  7
};


const PylithScalar pylith::faults::CohesiveKinSrcsDataLine2::_residualIncr[] = {
   0.0,
   7.5,
   0.0,
  -7.5,
  -0.2+1.89546413727+0.99414665414,
};

const PylithScalar pylith::faults::CohesiveKinSrcsDataLine2::_residual[] = {
   0.0,
   7.5,
   0.0,
  -7.5,
  -0.2+1.89546413727+0.99414665414,
};

const PylithScalar pylith::faults::CohesiveKinSrcsDataLine2::_jacobian[] = {
  0.0,  0.0,  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,  0.0, -1.0,
  0.0,  0.0,  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,  0.0, +1.0,
  0.0, -1.0,  0.0, +1.0,  0.0,
};

pylith::faults::CohesiveKinSrcsDataLine2::CohesiveKinSrcsDataLine2(void)
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
  orientation = const_cast<PylithScalar*>(_orientation);
  area = const_cast<PylithScalar*>(_area);
  residualIncr = const_cast<PylithScalar*>(_residualIncr);
  residual = const_cast<PylithScalar*>(_residual);
  jacobian = const_cast<PylithScalar*>(_jacobian);
  verticesFault = const_cast<int*>(_verticesFault);
  verticesLagrange = const_cast<int*>(_verticesLagrange);
  verticesNegative = const_cast<int*>(_verticesNegative);
  verticesPositive = const_cast<int*>(_verticesPositive);
  numFaultVertices = _numFaultVertices;  
  cellMappingFault = const_cast<int*>(_cellMappingFault);
  cellMappingCohesive = const_cast<int*>(_cellMappingCohesive);
  numCohesiveCells = _numCohesiveCells;  
} // constructor

pylith::faults::CohesiveKinSrcsDataLine2::~CohesiveKinSrcsDataLine2(void)
{}


// End of file
