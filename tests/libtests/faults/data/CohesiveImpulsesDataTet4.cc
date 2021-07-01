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
 * Cells are 0-1, vertices are 2-6.
 *
 * 2   3,4,5  6
 *
 *     ^^^^^ Face in x-y plane
 *
 * After adding cohesive elements
 *
 * Cells are 0-1,2, vertices are 3-10,11-13.
 *
 * 3   4,5,6  8,9,10   7
 *             11,12,13
 *     ^^^^^^^^^^^^ Cohesive element in x-y plane.
 */

#include "CohesiveImpulsesDataTet4.hh"

const char* pylith::faults::CohesiveImpulsesDataTet4::_meshFilename =
  "data/tet4.mesh";

const int pylith::faults::CohesiveImpulsesDataTet4::_spaceDim = 3;

const int pylith::faults::CohesiveImpulsesDataTet4::_cellDim = 2;

const int pylith::faults::CohesiveImpulsesDataTet4::_numBasis = 3;

const int pylith::faults::CohesiveImpulsesDataTet4::_numQuadPts = 3;

const PylithScalar pylith::faults::CohesiveImpulsesDataTet4::_quadPts[] = {
 -1.00000000e+00, -1.00000000e+00,
  1.00000000e+00, -1.00000000e+00,
 -1.00000000e+00,  1.00000000e+00,
};

const PylithScalar pylith::faults::CohesiveImpulsesDataTet4::_quadWts[] = {
  2.0/3.0, 2.0/3.0, 2.0/3.0,
};

const PylithScalar pylith::faults::CohesiveImpulsesDataTet4::_basis[] = {
  1.0, 0.0, 0.0,
  0.0, 1.0, 0.0,
  0.0, 0.0, 1.0,
};

const PylithScalar pylith::faults::CohesiveImpulsesDataTet4::_basisDeriv[] = {
 -0.50000000e+00, -0.50000000e+00,
  0.50000000e+00,  0.00000000e+00,
  0.00000000e+00,  0.50000000e+00,
 -0.50000000e+00, -0.50000000e+00,
  0.50000000e+00,  0.00000000e+00,
  0.00000000e+00,  0.50000000e+00,
 -0.50000000e+00, -0.50000000e+00,
  0.50000000e+00,  0.00000000e+00,
  0.00000000e+00,  0.50000000e+00,
};

const PylithScalar pylith::faults::CohesiveImpulsesDataTet4::_verticesRef[] = {
 -1.00000000e+00, -1.00000000e+00,
  1.00000000e+00, -1.00000000e+00,
 -1.00000000e+00,  1.00000000e+00,
};

const int pylith::faults::CohesiveImpulsesDataTet4::_id = 10;

const char* pylith::faults::CohesiveImpulsesDataTet4::_label = "fault";

const char* pylith::faults::CohesiveImpulsesDataTet4::_impulseAmpFilename = 
  "data/tet4_impulses.spatialdb";

const int pylith::faults::CohesiveImpulsesDataTet4::_impulseDOF[2] = {
  0, 2,
};
const int pylith::faults::CohesiveImpulsesDataTet4::_numComponents = 2;

const PylithScalar pylith::faults::CohesiveImpulsesDataTet4::_fieldT[] = {
  7.1, 8.1, 9.1,
  7.2, 8.2, 9.2, // 4
  7.3, 8.3, 9.3, // 5
  7.4, 8.4, 9.4, // 6
  7.5, 8.5, 9.5,
  7.6, 8.6, 9.6, // 8
  7.8, 8.8, 9.8, // 9
  7.0, 8.0, 9.0, // 10
  7.7, 8.7, 9.7, // 11
  7.9, 8.9, 9.9, // 12
  7.1, 8.1, 9.1, // 13
};

const PylithScalar pylith::faults::CohesiveImpulsesDataTet4::_fieldIncr[] = {
  3.1, 4.1, 5.1,
  3.2, 4.2, 5.2, // 4
  3.3, 4.3, 5.3, // 5
  3.4, 4.4, 5.4, // 6
  3.5, 4.5, 5.5,
  3.6, 4.6, 5.6, // 8
  3.8, 4.8, 5.8, // 9
  3.0, 4.0, 5.0, // 10
  3.7, 4.7, 5.7, // 11
  3.9, 4.9, 5.9, // 12
  3.1, 4.1, 5.1, // 13
};

// ----------------------------------------------------------------------
// Computed values
// ----------------------------------------------------------------------

const PylithScalar pylith::faults::CohesiveImpulsesDataTet4::_orientation[] = {
  0.0, +1.0, 0.0,    0.0, 0.0, +1.0,    +1.0, 0.0, 0.0,
  0.0, +1.0, 0.0,    0.0, 0.0, +1.0,    +1.0, 0.0, 0.0,
  0.0, +1.0, 0.0,    0.0, 0.0, +1.0,    +1.0, 0.0, 0.0,
};

const PylithScalar pylith::faults::CohesiveImpulsesDataTet4::_area[] = {
  1.0/3.0, 
  1.0/3.0, 
  1.0/3.0,
};

const PylithScalar pylith::faults::CohesiveImpulsesDataTet4::_amplitude[3] = {
  1.2,
  0.0,
  1.5,
};

const int pylith::faults::CohesiveImpulsesDataTet4::_numImpulses = 4;

const int pylith::faults::CohesiveImpulsesDataTet4::_numConstraintEdges = 3;
const int pylith::faults::CohesiveImpulsesDataTet4::_constraintEdges[] = {
  11, 12, 13,
};
const int pylith::faults::CohesiveImpulsesDataTet4::_negativeVertices[] = {
   4,  5,  6
};

const PylithScalar pylith::faults::CohesiveImpulsesDataTet4::_residual[] = {
  0.0,  0.0,  0.0,
   +7.7/3.0,  +8.7/3.0,  +9.7/3.0, // 4
  +7.9/3.0,  +8.9/3.0,  +9.9/3.0, // 5
  +7.1/3.0,  +8.1/3.0,  +9.1/3.0, // 6
  0.0,  0.0,  0.0,
  -7.7/3.0,  -8.7/3.0,  -9.7/3.0, // 8
  -7.9/3.0,  -8.9/3.0,  -9.9/3.0, // 9
  -7.1/3.0,  -8.1/3.0,  -9.1/3.0, // 10
  -1.0/3.0*(7.6-7.2 - 1.2), // 11
  -1.0/3.0*(8.6-8.2 + 0.0),
  -1.0/3.0*(9.6-9.2 + 0.0),
  -1.0/3.0*(7.8-7.3 + 0.0), // 12
  -1.0/3.0*(8.8-8.3 + 0.0),
  -1.0/3.0*(9.8-9.3 + 0.0),
  -1.0/3.0*(7.0-7.4 + 0.0), // 13
  -1.0/3.0*(8.0-8.4 + 0.0),
  -1.0/3.0*(9.0-9.4 + 0.0),
};

pylith::faults::CohesiveImpulsesDataTet4::CohesiveImpulsesDataTet4(void)
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
  residual = const_cast<PylithScalar*>(_residual);
  constraintEdges = const_cast<int*>(_constraintEdges);
  negativeVertices = const_cast<int*>(_negativeVertices);
  numConstraintEdges = _numConstraintEdges;  
} // constructor

pylith::faults::CohesiveImpulsesDataTet4::~CohesiveImpulsesDataTet4(void)
{}


// End of file
