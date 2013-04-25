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
// Copyright (c) 2010-2013 University of California, Davis
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
 * Cells are 0-1,2 vertices are 3-10,11-12.
 *
 *       4 -------- 6 -12--10 -------- 8
 *       |          |       |          |
 *       |          |       |          |
 *       |          |       |          |
 *       |          |       |          |
 *       |          |       |          |
 *       |          |       |          |
 *       |          |       |          |
 *       |          |       |          |
 *       3 -------- 5 -11-- 9 -------- 7
 */

#include "CohesiveImpulsesDataQuad4.hh"

const char* pylith::faults::CohesiveImpulsesDataQuad4::_meshFilename =
  "data/quad4.mesh";

const int pylith::faults::CohesiveImpulsesDataQuad4::_spaceDim = 2;

const int pylith::faults::CohesiveImpulsesDataQuad4::_cellDim = 1;

const int pylith::faults::CohesiveImpulsesDataQuad4::_numBasis = 2;

const int pylith::faults::CohesiveImpulsesDataQuad4::_numQuadPts = 2;

const PylithScalar pylith::faults::CohesiveImpulsesDataQuad4::_quadPts[] = {
  -1.0, 1.0,
};

const PylithScalar pylith::faults::CohesiveImpulsesDataQuad4::_quadWts[] = {
  1.0, 1.0
};

const PylithScalar pylith::faults::CohesiveImpulsesDataQuad4::_basis[] = {
  1.0, 0.0,
  0.0, 1.0,
};

const PylithScalar pylith::faults::CohesiveImpulsesDataQuad4::_basisDeriv[] = {
  -0.5, 0.5,
  -0.5, 0.5,
};

const PylithScalar pylith::faults::CohesiveImpulsesDataQuad4::_verticesRef[] = {
  -1.0, 1.0
};

const int pylith::faults::CohesiveImpulsesDataQuad4::_id = 10;

const char* pylith::faults::CohesiveImpulsesDataQuad4::_label = "fault";

const char* pylith::faults::CohesiveImpulsesDataQuad4::_impulseAmpFilename = 
  "data/quad4_impulses.spatialdb";

const int pylith::faults::CohesiveImpulsesDataQuad4::_impulseDOF[2] = {
  0,1,
};
const int pylith::faults::CohesiveImpulsesDataQuad4::_numComponents = 2;

const PylithScalar pylith::faults::CohesiveImpulsesDataQuad4::_fieldT[] = {
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

const PylithScalar pylith::faults::CohesiveImpulsesDataQuad4::_fieldIncr[] = {
  3.1, 4.1,
  3.2, 4.2,
  3.3, 4.3, // 4
  3.4, 4.4, // 5
  3.5, 4.5,
  3.6, 4.6,
  3.7, 4.7, // 8
  3.9, 4.9, // 9
  3.8, 4.8, // 10
  3.0, 4.0, // 11
};

// ----------------------------------------------------------------------
// Computed values
// ----------------------------------------------------------------------

const PylithScalar pylith::faults::CohesiveImpulsesDataQuad4::_orientation[] = {
  0.0,  1.0,  +1.0, 0.0,
  0.0,  1.0,  +1.0, 0.0
};

const PylithScalar pylith::faults::CohesiveImpulsesDataQuad4::_area[] = {
  1.0,
  1.0,
};

const PylithScalar pylith::faults::CohesiveImpulsesDataQuad4::_amplitude[] = {
  0.0,
  1.5,
};

const int pylith::faults::CohesiveImpulsesDataQuad4::_numImpulses = 2;

const int pylith::faults::CohesiveImpulsesDataQuad4::_numConstraintVert = 2;
const int pylith::faults::CohesiveImpulsesDataQuad4::_constraintVertices[] = {
  11, 12
};
const int pylith::faults::CohesiveImpulsesDataQuad4::_negativeVertices[] = {
   5,  6
};

const PylithScalar pylith::faults::CohesiveImpulsesDataQuad4::_residualIncr[] = {
  0.0,  0.0,
  0.0,  0.0,
 +8.8, +9.8, // 4
 +8.0, +9.0, // 5
  0.0,  0.0,
  0.0,  0.0,
 -8.8, -9.8, // 8
 -8.0, -9.0, // 9
  -(8.7-8.3) + 0.0,
  -(9.7-9.3) + 0.0, // 10
  -(8.9-8.4) + 1.5,
  -(9.9-9.4) + 0.0, // 11
};

pylith::faults::CohesiveImpulsesDataQuad4::CohesiveImpulsesDataQuad4(void)
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
  impulseAmpFilename = const_cast<char*>(_impulseAmpFilename);
  impulseDOF = const_cast<int*>(_impulseDOF);
  numComponents = _numComponents;
  fieldT = const_cast<PylithScalar*>(_fieldT);
  fieldIncr = const_cast<PylithScalar*>(_fieldIncr);
  orientation = const_cast<PylithScalar*>(_orientation);
  area = const_cast<PylithScalar*>(_area);
  amplitude = const_cast<PylithScalar*>(_amplitude);
  numImpulses = _numImpulses;
  residualIncr = const_cast<PylithScalar*>(_residualIncr);
  constraintVertices = const_cast<int*>(_constraintVertices);
  negativeVertices = const_cast<int*>(_negativeVertices);
  numConstraintVert = _numConstraintVert;  
} // constructor

pylith::faults::CohesiveImpulsesDataQuad4::~CohesiveImpulsesDataQuad4(void)
{}


// End of file
