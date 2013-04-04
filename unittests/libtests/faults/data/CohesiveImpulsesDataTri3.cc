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
 * Cells are 0-1, 2, vertices are 3-8,9-10.
 *
 *              7 -9- 4
 *             /|     |\
 *            / |     | \
 *           /  |     |  \
 *          /   |     |   \
 *         3    |     |    6
 *          \   |     |   /
 *           \  |     |  /
 *            \ |     | /
 *             \|     |/
 *              8-10- 5
 */

#include "CohesiveImpulsesDataTri3.hh"

const char* pylith::faults::CohesiveImpulsesDataTri3::_meshFilename =
  "data/tri3.mesh";

const int pylith::faults::CohesiveImpulsesDataTri3::_spaceDim = 2;

const int pylith::faults::CohesiveImpulsesDataTri3::_cellDim = 1;

const int pylith::faults::CohesiveImpulsesDataTri3::_numBasis = 2;

const int pylith::faults::CohesiveImpulsesDataTri3::_numQuadPts = 2;

const PylithScalar pylith::faults::CohesiveImpulsesDataTri3::_quadPts[] = {
  -1.0, 1.0,
};

const PylithScalar pylith::faults::CohesiveImpulsesDataTri3::_quadWts[] = {
  1.0, 1.0
};

const PylithScalar pylith::faults::CohesiveImpulsesDataTri3::_basis[] = {
  1.0, 0.0,
  0.0, 1.0,
};

const PylithScalar pylith::faults::CohesiveImpulsesDataTri3::_basisDeriv[] = {
  -0.5, 0.5,
  -0.5, 0.5,
};

const PylithScalar pylith::faults::CohesiveImpulsesDataTri3::_verticesRef[] = {
  -1.0, 1.0
};

const int pylith::faults::CohesiveImpulsesDataTri3::_id = 10;

const char* pylith::faults::CohesiveImpulsesDataTri3::_label = "fault";

const char* pylith::faults::CohesiveImpulsesDataTri3::_impulseAmpFilename = 
  "data/tri3_impulses.spatialdb";

const int pylith::faults::CohesiveImpulsesDataTri3::_impulseDOF[1] = {
  1,
};
const int pylith::faults::CohesiveImpulsesDataTri3::_numComponents = 1;

const PylithScalar pylith::faults::CohesiveImpulsesDataTri3::_fieldT[] = {
  8.1, 9.1,
  8.2, 9.2, // 3
  8.3, 9.3, // 4
  8.4, 9.4,
  8.5, 9.5, // 6
  8.7, 9.7, // 7
  8.6, 9.6, // 8
  8.8, 9.8, // 9
};

const PylithScalar pylith::faults::CohesiveImpulsesDataTri3::_fieldIncr[] = {
  3.1, 4.1,
  3.2, 4.2, // 3
  3.3, 4.3, // 4
  3.4, 4.4,
  3.5, 4.5, // 6
  3.7, 4.7, // 7
  3.6, 4.6, // 8
  3.8, 4.8, // 9
};

// ----------------------------------------------------------------------
// Computed values
// ----------------------------------------------------------------------

const PylithScalar pylith::faults::CohesiveImpulsesDataTri3::_orientation[] = {
  0.0, -1.0,  -1.0, 0.0,
  0.0, -1.0,  -1.0, 0.0
};

const PylithScalar pylith::faults::CohesiveImpulsesDataTri3::_area[] = {
  1.0,
  1.0,
};

const PylithScalar pylith::faults::CohesiveImpulsesDataTri3::_amplitude[] = {
  2.0,
  2.1,
};

const int pylith::faults::CohesiveImpulsesDataTri3::_numImpulses = 2;

const int pylith::faults::CohesiveImpulsesDataTri3::_numConstraintVert = 2;
const int pylith::faults::CohesiveImpulsesDataTri3::_constraintVertices[] = {
  9, 10
};
const int pylith::faults::CohesiveImpulsesDataTri3::_negativeVertices[] = {
  4, 5
};

const PylithScalar pylith::faults::CohesiveImpulsesDataTri3::_residualIncr[] = {
  0.0,  0.0,
 +8.6, +9.6, // 3
 +8.8, +9.8, // 4
  0.0,  0.0,
 -8.6, -9.6, // 6
 -8.8, -9.8, // 7
 -(8.5-8.2) - (0),
 -(9.5-9.2) - (0), // 8
 -(8.7-8.3) - (2.1),
 -(9.7-9.3) - (0), // 9
};

pylith::faults::CohesiveImpulsesDataTri3::CohesiveImpulsesDataTri3(void)
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

pylith::faults::CohesiveImpulsesDataTri3::~CohesiveImpulsesDataTri3(void)
{}


// End of file
