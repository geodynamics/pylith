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
// Copyright (c) 2010-2012 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

/* Original mesh
 *
 * Cells are 0-1 and vertices are 2-13.
 *
 *       2,3,4,5 -------- 6,7,8,9 -------- 10,11,12,13
 *
 *                        ^^^^^^^ Vertices forming fault
 *
 * After adding cohesive elements
 *
 * Cells are 0-1,2 and vertices are 3-18,19-22.
 *
 *       3,4,5,6 -------- 7,8,9,10 -- 15,16,17,18 -------- 11,12,13,14
 *                                    19,20,21,22
 *                        ^^^^^^^^^^^^^^^^^^^^^^ Cohesive element
 *
 */

#include "CohesiveImpulsesDataHex8.hh"

const char* pylith::faults::CohesiveImpulsesDataHex8::_meshFilename =
  "data/hex8.mesh";

const int pylith::faults::CohesiveImpulsesDataHex8::_spaceDim = 3;

const int pylith::faults::CohesiveImpulsesDataHex8::_cellDim = 2;

const int pylith::faults::CohesiveImpulsesDataHex8::_numBasis = 4;

const int pylith::faults::CohesiveImpulsesDataHex8::_numQuadPts = 4;

const PylithScalar pylith::faults::CohesiveImpulsesDataHex8::_quadPts[] = {
  -1.0, -1.0,
  +1.0, -1.0,
  +1.0, +1.0,
  -1.0, +1.0
};

const PylithScalar pylith::faults::CohesiveImpulsesDataHex8::_quadWts[] = {
  1.0, 1.0, 1.0, 1.0
};

const PylithScalar pylith::faults::CohesiveImpulsesDataHex8::_basis[] = {
  1.0, 0.0, 0.0, 0.0,
  0.0, 1.0, 0.0, 0.0,
  0.0, 0.0, 1.0, 0.0,
  0.0, 0.0, 0.0, 1.0,
};

const PylithScalar pylith::faults::CohesiveImpulsesDataHex8::_basisDeriv[] = {
  -0.39433757, -0.39433757,
  +0.39433757, -0.10566243,
  +0.10566243, +0.10566243,
  -0.10566243, +0.39433757,

  -0.39433757, -0.10566243,
  +0.39433757, -0.39433757,
  +0.10566243, +0.39433757,
  -0.10566243, +0.10566243,

  -0.10566243, -0.10566243,
  +0.10566243, -0.39433757,
  +0.39433757, +0.39433757,
  -0.39433757, +0.10566243,

  -0.10566243, -0.39433757,
  +0.10566243, -0.10566243,
  +0.39433757, +0.10566243,
  -0.39433757, +0.39433757,
};

const PylithScalar pylith::faults::CohesiveImpulsesDataHex8::_verticesRef[] = {
  -1.0, -1.0,
  +1.0, -1.0,
  +1.0, +1.0,
  -1.0, +1.0
};

const int pylith::faults::CohesiveImpulsesDataHex8::_id = 10;

const char* pylith::faults::CohesiveImpulsesDataHex8::_label = "fault";

const char* pylith::faults::CohesiveImpulsesDataHex8::_impulseAmpFilename = 
  "data/hex8_impulses.spatialdb";

const int pylith::faults::CohesiveImpulsesDataHex8::_impulseDOF[1] = {
  1,
};
const int pylith::faults::CohesiveImpulsesDataHex8::_numComponents = 1;

const PylithScalar pylith::faults::CohesiveImpulsesDataHex8::_fieldT[] = {
  4.1, 6.1, 8.1,
  4.2, 6.2, 8.2,
  4.3, 6.3, 8.3,
  4.4, 6.4, 8.4,
  4.5, 6.5, 8.5, // 6
  4.6, 6.6, 8.6, // 7
  4.7, 6.7, 8.7, // 8
  4.8, 6.8, 8.8, // 9
  4.9, 6.9, 8.9,
  4.0, 6.0, 8.0,
  5.1, 7.1, 9.1,
  5.2, 7.2, 9.2,
  5.3, 7.3, 9.3, // 14
  5.5, 7.5, 9.5, // 15
  5.7, 7.7, 9.7, // 16
  5.9, 7.9, 9.9, // 17
  5.4, 7.4, 9.4, // 18
  5.6, 7.6, 9.6, // 19
  5.8, 7.8, 9.8, // 20
  5.0, 7.0, 9.0, // 21
};

const PylithScalar pylith::faults::CohesiveImpulsesDataHex8::_fieldIncr[] = {
  3.1, 4.1, 5.1,
  3.2, 4.2, 5.2,
  3.3, 4.3, 5.3,
  3.4, 4.4, 5.4,
  3.5, 4.5, 5.5, // 6
  3.6, 4.6, 5.6, // 7
  3.7, 4.7, 5.7, // 8
  3.8, 4.8, 5.8, // 9
  3.9, 4.9, 5.9,
  3.0, 4.0, 5.0,
  3.1, 4.1, 5.1,
  3.2, 4.2, 5.2,
  3.3, 4.3, 5.3, // 14
  3.5, 4.5, 5.5, // 15
  3.7, 4.7, 5.7, // 16
  3.9, 4.9, 5.9, // 17
  3.4, 4.4, 5.4, // 18
  3.6, 4.6, 5.6, // 19
  3.8, 4.8, 5.8, // 20
  3.0, 4.0, 5.0, // 21
};

// ----------------------------------------------------------------------
// Computed values
// ----------------------------------------------------------------------

const PylithScalar pylith::faults::CohesiveImpulsesDataHex8::_orientation[] = {
  0.0, -1.0, 0.0,    0.0, 0.0, +1.0,    -1.0, 0.0, 0.0,
  0.0, -1.0, 0.0,    0.0, 0.0, +1.0,    -1.0, 0.0, 0.0,
  0.0, -1.0, 0.0,    0.0, 0.0, +1.0,    -1.0, 0.0, 0.0,
  0.0, -1.0, 0.0,    0.0, 0.0, +1.0,    -1.0, 0.0, 0.0,
};

const PylithScalar pylith::faults::CohesiveImpulsesDataHex8::_area[] = {
  1.0, 1.0, 1.0, 1.0
};

const PylithScalar pylith::faults::CohesiveImpulsesDataHex8::_amplitude[4] = {
  0.8,
  1.6,
  0.0,
  1.2,
};

const int pylith::faults::CohesiveImpulsesDataHex8::_numImpulses = 3;

const int pylith::faults::CohesiveImpulsesDataHex8::_numConstraintVert = 4;
const int pylith::faults::CohesiveImpulsesDataHex8::_constraintVertices[4] = {
  19, 20, 21, 22
};

const PylithScalar pylith::faults::CohesiveImpulsesDataHex8::_residualIncr[] = {
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  +5.4, +7.4, +9.4, // 6
  +5.6, +7.6, +9.6, // 7
  +5.8, +7.8, +9.8, // 8
  +5.0, +7.0, +9.0, // 9
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  -5.4, -7.4, -9.4, // 14
  -5.6, -7.6, -9.6, // 15
  -5.8, -7.8, -9.8, // 16
  -5.0, -7.0, -9.0, // 17

  // 18 (constraint)
  -(5.3-4.5+0.0), -(7.3-6.5+0.0), -(9.3-8.5+0.0),

  // 19 (constraint)
  -(5.5-4.6+0.0), -(7.5-6.6+0.0), -(9.5-8.6-1.6),

  // 20 (constraint)
  -(5.7-4.7+0.0), -(7.7-6.7+0.0), -(9.7-8.7+0.0),

  // 21 (constraint)
  -(5.9-4.8+0.0), -(7.9-6.8+0.0), -(9.9-8.8+0.0),
};

pylith::faults::CohesiveImpulsesDataHex8::CohesiveImpulsesDataHex8(void)
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
  numConstraintVert = _numConstraintVert;  
} // constructor

pylith::faults::CohesiveImpulsesDataHex8::~CohesiveImpulsesDataHex8(void)
{}


// End of file
