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

/* Mesh: meshTet4.txt
 *
 * DirichletPoints BC A at vertices 1 and 3.
 *
 * Fixed DOF: { 0 }
 *
 * Initial values
 *   0: 0.1
 *   2: 0.2
 *   3: 0.3
 *   4: 0.4
 * tRef = 0.0
 * Rate of change of values
 *   0: 1.0
 *   2: 2.0
 *   3: 3.0
 *   4: 4.0
 *
 * DirichletPoints BC B at vertices 2 and 3.
 *
 * Fixed DOF: { 2 }
 *
 * Initial values
 *   0: 0.01
 *   2: 0.02
 *   3: 0.03
 * tRef = 0.0
 * Rate of change of values
 *   0: -1.0
 *   2: -2.0
 *   3: -3.0
 *
 * DirichletPoints BC C at vertices 2 and 3.
 *
 * Fixed DOF: { 1 }
 *
 * Initial values
 *   1: 11.0
 *   2: 22.0
 *   3: 33.0
 * tRef = 0.0
 * Rate of change of values
 *   1: 10.0
 *   2: 20.0
 *   3: 30.0
 */

#include "DirichletDataMultiTet4.hh"

const int pylith::bc::DirichletDataMultiTet4::_numDOF = 3;

const int pylith::bc::DirichletDataMultiTet4::_idA = 0;
const char* pylith::bc::DirichletDataMultiTet4::_labelA = "bc4";
const int pylith::bc::DirichletDataMultiTet4::_numFixedDOFA = 1;
const int pylith::bc::DirichletDataMultiTet4::_fixedDOFA[] = { 0 };
const int pylith::bc::DirichletDataMultiTet4::_numConstrainedPtsA = 4;
const int pylith::bc::DirichletDataMultiTet4::_constrainedPointsA[] = { 0, 2, 3, 4 };

const char* pylith::bc::DirichletDataMultiTet4::_dbFilenameA =
  "data/tet4_disp2.spatialdb";
const char* pylith::bc::DirichletDataMultiTet4::_dbFilenameARate =
  "data/tet4_vel2.spatialdb";
const PylithScalar pylith::bc::DirichletDataMultiTet4::_tRefA = 0.0;

const int pylith::bc::DirichletDataMultiTet4::_idB = 1;
const char* pylith::bc::DirichletDataMultiTet4::_labelB = "bc2";
const int pylith::bc::DirichletDataMultiTet4::_numFixedDOFB = 1;
const int pylith::bc::DirichletDataMultiTet4::_fixedDOFB[] = { 2 };
const int pylith::bc::DirichletDataMultiTet4::_numConstrainedPtsB = 3;
const int pylith::bc::DirichletDataMultiTet4::_constrainedPointsB[] = { 0, 2, 3 };

const char* pylith::bc::DirichletDataMultiTet4::_dbFilenameB =
  "data/tet4_disp2.spatialdb";
const char* pylith::bc::DirichletDataMultiTet4::_dbFilenameBRate =
  "data/tet4_vel2.spatialdb";
const PylithScalar pylith::bc::DirichletDataMultiTet4::_tRefB = 0.0;

const int pylith::bc::DirichletDataMultiTet4::_idC = 1;
const char* pylith::bc::DirichletDataMultiTet4::_labelC = "bc3";
const int pylith::bc::DirichletDataMultiTet4::_numFixedDOFC = 1;
const int pylith::bc::DirichletDataMultiTet4::_fixedDOFC[] = { 1 };
const int pylith::bc::DirichletDataMultiTet4::_numConstrainedPtsC = 3;
const int pylith::bc::DirichletDataMultiTet4::_constrainedPointsC[] = { 1, 2, 3 };

const char* pylith::bc::DirichletDataMultiTet4::_dbFilenameC =
  "data/tet4_disp2.spatialdb";
const char* pylith::bc::DirichletDataMultiTet4::_dbFilenameCRate =
  "data/tet4_vel2.spatialdb";
const PylithScalar pylith::bc::DirichletDataMultiTet4::_tRefC = 0.0;

const int pylith::bc::DirichletDataMultiTet4::_constraintSizes[] = {
  2,
  1,
  3,
  3,
  1
};

const int pylith::bc::DirichletDataMultiTet4::_constrainedDOF[] = {
  0, 2, 
  1,
  0, 1, 2,
  0, 1, 2,
  0
};

// Values at t=10.0
const PylithScalar pylith::bc::DirichletDataMultiTet4::_field[] = {
  10.1, 0.0, -9.99,
   0.0, 111.0, 0.0,
  20.2, 222.0, -19.98,
  30.3, 333.0, -29.97,
  40.4, 0.0, 0.0
};

// Increment values from t=10.0 to t=14.0
const PylithScalar pylith::bc::DirichletDataMultiTet4::_fieldIncr[] = {
   4.0,   0.0,  -4.0,
   0.0,  40.0,   0.0,
   8.0,  80.0,  -8.0,
  12.0, 120.0, -12.0,
  16.0,   0.0,   0.0
};

const char* pylith::bc::DirichletDataMultiTet4::_meshFilename = 
  "data/tet4.mesh";

pylith::bc::DirichletDataMultiTet4::DirichletDataMultiTet4(void)
{ // constructor
  numDOF = _numDOF;

  idA = _idA;
  labelA = const_cast<char*>(_labelA);
  numFixedDOFA = _numFixedDOFA;
  fixedDOFA = const_cast<int*>(_fixedDOFA);
  numConstrainedPtsA = _numConstrainedPtsA;
  constrainedPointsA = const_cast<int*>(_constrainedPointsA);

  dbFilenameA = const_cast<char*>(_dbFilenameA);
  dbFilenameARate = const_cast<char*>(_dbFilenameARate);
  tRefA = _tRefA;

  idB = _idB;
  labelB = const_cast<char*>(_labelB);
  numFixedDOFB = _numFixedDOFB;
  fixedDOFB = const_cast<int*>(_fixedDOFB);
  numConstrainedPtsB = _numConstrainedPtsB;
  constrainedPointsB = const_cast<int*>(_constrainedPointsB);

  dbFilenameB = const_cast<char*>(_dbFilenameB);
  dbFilenameBRate = const_cast<char*>(_dbFilenameBRate);
  tRefB = _tRefB;

  idC = _idC;
  labelC = const_cast<char*>(_labelC);
  numFixedDOFC = _numFixedDOFC;
  fixedDOFC = const_cast<int*>(_fixedDOFC);
  numConstrainedPtsC = _numConstrainedPtsC;
  constrainedPointsC = const_cast<int*>(_constrainedPointsC);

  dbFilenameC = const_cast<char*>(_dbFilenameC);
  dbFilenameCRate = const_cast<char*>(_dbFilenameCRate);
  tRefC = _tRefC;

  field = const_cast<PylithScalar*>(_field);
  fieldIncr = const_cast<PylithScalar*>(_fieldIncr);
  constraintSizes = const_cast<int*>(_constraintSizes);
  constrainedDOF = const_cast<int*>(_constrainedDOF);

  meshFilename = const_cast<char*>(_meshFilename);
} // constructor

pylith::bc::DirichletDataMultiTet4::~DirichletDataMultiTet4(void)
{}


// End of file
