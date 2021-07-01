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
 */


#include "CohesiveKinDataTri3g.hh"

const char* pylith::faults::CohesiveKinDataTri3g::_meshFilename = "data/tri3g.mesh";

const int pylith::faults::CohesiveKinDataTri3g::_spaceDim = 2;

const int pylith::faults::CohesiveKinDataTri3g::_cellDim = 1;

const int pylith::faults::CohesiveKinDataTri3g::_numBasis = 2;

const int pylith::faults::CohesiveKinDataTri3g::_numQuadPts = 2;

const PylithScalar pylith::faults::CohesiveKinDataTri3g::_quadPts[2] = {
  -1.0, 1.0,
};

const PylithScalar pylith::faults::CohesiveKinDataTri3g::_quadWts[2] = {
  1.0, 1.0
};

const PylithScalar pylith::faults::CohesiveKinDataTri3g::_basis[2*2] = {
  1.0, 0.0,
  0.0, 1.0,
};

const PylithScalar pylith::faults::CohesiveKinDataTri3g::_basisDeriv[2*2] = {
  -0.5, 0.5,
  -0.5, 0.5,
};

const PylithScalar pylith::faults::CohesiveKinDataTri3g::_verticesRef[2] = {
  -1.0, 1.0
};

const int pylith::faults::CohesiveKinDataTri3g::_id = 10;

const char* pylith::faults::CohesiveKinDataTri3g::_label = "fault";
const char* pylith::faults::CohesiveKinDataTri3g::_edge = "edge";

const char* pylith::faults::CohesiveKinDataTri3g::_finalSlipFilename = 
  "data/tri3g_finalslip.spatialdb";

const char* pylith::faults::CohesiveKinDataTri3g::_slipTimeFilename = 
  "data/tri3g_sliptime.spatialdb";

const char* pylith::faults::CohesiveKinDataTri3g::_riseTimeFilename = 
  "data/tri3g_risetime.spatialdb";

const PylithScalar pylith::faults::CohesiveKinDataTri3g::_fieldT[(8+1)*2] = {
  6.1, 8.1, // 7
  6.2, 8.2, // 8
  6.3, 8.3, // 9
  6.4, 8.4, // 10
  6.5, 8.5, // 11
  6.6, 8.6, // 12
  6.7, 8.7, // 13
  6.9, 8.9, // 14
  7.1, 9.1, // 28
};

const PylithScalar pylith::faults::CohesiveKinDataTri3g::_fieldIncr[(8+1)*2] = {
  3.1, 7.1, // 7
  3.2, 7.2, // 8
  3.3, 7.3, // 9
  3.4, 7.4, // 10
  3.5, 7.5, // 11
  3.6, 7.6, // 12
  3.7, 7.7, // 13
  3.9, 7.9, // 14
  3.1, 7.1, // 28
};

const PylithScalar pylith::faults::CohesiveKinDataTri3g::_jacobianLumped[(8+1)*2] = {
  6.1, 8.1, // 7
  6.2, 8.2, // 8
  6.3, 8.3, // 9
  6.4, 8.4, // 10
  6.5, 8.5, // 11
  6.6, 8.6, // 12
  6.7, 8.7, // 13
  6.9, 8.9, // 14
  7.1, 9.1, // 28
};


const PylithScalar pylith::faults::CohesiveKinDataTri3g::_orientation[2*2*2] = {
  0.0, +1.0,  +1.0,  0.0,
  0.0, +1.0,  +1.0,  0.0,
};

const PylithScalar pylith::faults::CohesiveKinDataTri3g::_area[2] = {
  1.0,
  1.0,
};

const int pylith::faults::CohesiveKinDataTri3g::_numFaultVertices = 2;
const int pylith::faults::CohesiveKinDataTri3g::_verticesFault[2] = {
  2,  3,
};
const int pylith::faults::CohesiveKinDataTri3g::_edgesLagrange[2] = {
  29, 28,
};
const int pylith::faults::CohesiveKinDataTri3g::_verticesNegative[2] = {
   9, 12,
};
const int pylith::faults::CohesiveKinDataTri3g::_verticesPositive[2] = {
  9, 14,
};

const int pylith::faults::CohesiveKinDataTri3g::_numCohesiveCells = 1;
const int pylith::faults::CohesiveKinDataTri3g::_cellMappingFault[1] = {
  4,
};
const int pylith::faults::CohesiveKinDataTri3g::_cellMappingCohesive[1] = {
  6,
};

const PylithScalar pylith::faults::CohesiveKinDataTri3g::_residual[(8+1)*2] = {
  0.0,  0.0, // 7
  0.0,  0.0, // 8
  0.0,  0.0, // 9
  0.0,  0.0, // 10
  0.0,  0.0, // 11
 +7.1, +9.1, // 12
  0.0,  0.0, // 13
 -7.1, -9.1, // 14
 -(6.9-6.6) + (0.14794836271),
 -(8.9-8.6) + (1.77538035254), // 28
};

const PylithScalar pylith::faults::CohesiveKinDataTri3g::_jacobian[(8+1)*2*(8+1)*2] = {
  0.0, 0.0, // 7x
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
  0.0, 0.0, // 8x
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0, // 8y
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
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
  0.0, 0.0, // 9y
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
  0.0, 0.0, // 10y
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
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
  0.0, 0.0, // 11y
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0, // 12x
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
 -1.0, 0.0, // 28
  0.0, 0.0, // 12y
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0,-1.0, // 28
  0.0, 0.0, // 13x
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0, // 13y
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0, // 14x
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
 +1.0, 0.0, // 28
  0.0, 0.0, // 14y
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0,+1.0, // 28
  0.0, 0.0, // 28x
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
 -1.0, 0.0, // 12
  0.0, 0.0,
 +1.0, 0.0, // 14
  0.0, 0.0,
  0.0, 0.0, // 28y
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0, 0.0,
  0.0,-1.0, // 12
  0.0, 0.0,
  0.0,+1.0, // 14
  0.0, 0.0,
};

pylith::faults::CohesiveKinDataTri3g::CohesiveKinDataTri3g(void)
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
  edge = const_cast<char*>(_edge);
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

pylith::faults::CohesiveKinDataTri3g::~CohesiveKinDataTri3g(void)
{}


// End of file
