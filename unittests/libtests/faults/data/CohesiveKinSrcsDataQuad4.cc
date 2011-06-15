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
 *       3 -------- 5 -11-- 9 -------- 7
 *       |          |       |          |
 *       |          |       |          |
 *       |          |       |          |
 *       |          |       |          |
 *       |          |       |          |
 *       |          |       |          |
 *       |          |       |          |
 *       |          |       |          |
 *       2 -------- 4 -10-- 8 -------- 6
 */

#include "CohesiveKinSrcsDataQuad4.hh"

const char* pylith::faults::CohesiveKinSrcsDataQuad4::_meshFilename =
  "data/quad4.mesh";

const int pylith::faults::CohesiveKinSrcsDataQuad4::_spaceDim = 2;

const int pylith::faults::CohesiveKinSrcsDataQuad4::_cellDim = 1;

const int pylith::faults::CohesiveKinSrcsDataQuad4::_numBasis = 2;

const int pylith::faults::CohesiveKinSrcsDataQuad4::_numQuadPts = 1;

const double pylith::faults::CohesiveKinSrcsDataQuad4::_quadPts[] = {
  0.0,
};

const double pylith::faults::CohesiveKinSrcsDataQuad4::_quadWts[] = {
  2.0,
};

const double pylith::faults::CohesiveKinSrcsDataQuad4::_basis[] = {
  0.5,
  0.5
};

const double pylith::faults::CohesiveKinSrcsDataQuad4::_basisDeriv[] = {
  -0.5,
   0.5
};

const double pylith::faults::CohesiveKinSrcsDataQuad4::_verticesRef[] = {
  -1.0, 1.0
};

const int pylith::faults::CohesiveKinSrcsDataQuad4::_id = 10;

const char* pylith::faults::CohesiveKinSrcsDataQuad4::_label = "fault";

const char* pylith::faults::CohesiveKinSrcsDataQuad4::_finalSlipFilename = 
  "data/quad4_finalslip.spatialdb";

const char* pylith::faults::CohesiveKinSrcsDataQuad4::_slipTimeFilename = 
  "data/quad4_sliptime.spatialdb";

const char* pylith::faults::CohesiveKinSrcsDataQuad4::_riseTimeFilename = 
  "data/quad4_risetime.spatialdb";

const double pylith::faults::CohesiveKinSrcsDataQuad4::_fieldT[] = {
  8.1, 9.1,
  8.2, 9.2,
  8.3, 9.3, // 4
  8.4, 9.4, // 5
  8.5, 9.5,
  8.6, 9.6,
  8.7, 9.7, // 8
  8.9, 9.9, // 9
  8.8, 9.8, // 10
  8.0, 9.0, // 11
};



const double pylith::faults::CohesiveKinSrcsDataQuad4::_orientation[] = {
  0.0,  1.0,  +1.0, 0.0,
  0.0,  1.0,  +1.0, 0.0
};

const double pylith::faults::CohesiveKinSrcsDataQuad4::_area[] = {
  1.0, 1.0,
};

const int pylith::faults::CohesiveKinSrcsDataQuad4::_numFaultVertices = 2;
const int pylith::faults::CohesiveKinSrcsDataQuad4::_verticesFault[] = {
  1, 2
};
const int pylith::faults::CohesiveKinSrcsDataQuad4::_verticesLagrange[] = {
  10, 11
};
const int pylith::faults::CohesiveKinSrcsDataQuad4::_verticesNegative[] = {
  4, 5
};
const int pylith::faults::CohesiveKinSrcsDataQuad4::_verticesPositive[] = {
  8, 9
};

const int pylith::faults::CohesiveKinSrcsDataQuad4::_numCohesiveCells = 1;
const int pylith::faults::CohesiveKinSrcsDataQuad4::_cellMappingFault[] = {
  0
};
const int pylith::faults::CohesiveKinSrcsDataQuad4::_cellMappingCohesive[] = {
  12
};


const double pylith::faults::CohesiveKinSrcsDataQuad4::_residual[] = {
  0.0,  0.0,
  0.0,  0.0,
  9.8,  8.8, // 4
  9.0,  8.0, // 5
  0.0,  0.0,
  0.0,  0.0,
 -9.8, -8.8, // 8
 -9.0, -8.0, // 9
 -0.4+1.77538035254+0.68377062865, -0.4+0.14794836271+0.05698088572, // 10
 -0.5+1.89546413727+0.99414665414, -0.5+0.08241148423+0.04322376757, // 11
};

const double pylith::faults::CohesiveKinSrcsDataQuad4::_residualIncr[] = {
  0.0,  0.0,
  0.0,  0.0,
  9.8,  8.8, // 4
  9.0,  8.0, // 5
  0.0,  0.0,
  0.0,  0.0,
 -9.8, -8.8, // 8
 -9.0, -8.0, // 9
 -0.4+1.77538035254+0.68377062865, -0.4+0.14794836271+0.05698088572, // 10
 -0.5+1.89546413727+0.99414665414, -0.5+0.08241148423+0.04322376757, // 11
};

const double pylith::faults::CohesiveKinSrcsDataQuad4::_jacobian[] = {
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
  0.0, 0.0,
  0.0,-1.0, // 10
  0.0, 0.0,
  0.0, 0.0, // 4y
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
 -1.0, 0.0, //  10
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
  0.0, 0.0,
  0.0,+1.0, //  10
  0.0, 0.0,
  0.0, 0.0, // 8y
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
 +1.0, 0.0, // 10
  0.0, 0.0,
  0.0, 0.0, // 9x
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0,+1.0, // 11
  0.0, 0.0, // 9y
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
 +1.0, 0.0, // 11
  0.0, 0.0, // 10x
  0.0, 0.0,
  0.0,-1.0, //  4
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0,+1.0, //  8
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0, // 10y
  0.0, 0.0,
 -1.0, 0.0, //  4
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
 +1.0, 0.0, //  8
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0, // 11x
  0.0, 0.0,
  0.0, 0.0,
  0.0,-1.0, //  5
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0,+1.0, //  9
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0, // 11y
  0.0, 0.0,
  0.0, 0.0,
 -1.0, 0.0, //  5
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
 +1.0, 0.0, // 9
  0.0, 0.0,
  0.0, 0.0,
};

pylith::faults::CohesiveKinSrcsDataQuad4::CohesiveKinSrcsDataQuad4(void)
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

pylith::faults::CohesiveKinSrcsDataQuad4::~CohesiveKinSrcsDataQuad4(void)
{}


// End of file
