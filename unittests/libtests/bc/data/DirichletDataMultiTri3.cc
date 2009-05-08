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

/* Mesh: meshTri3.txt
 *
 * DirichletPoints BC A at vertices 1 and 3.
 *
 * Fixed DOF: { 1 }
 *
 * Initial values
 *   1: 0.3
 *   3: 0.7
 * tRef = 3.2
 * Rate of change of values
 *   1: 0.2
 *   3: 0.8
 *
 * DirichletPoints BC B at vertices 2 and 3.
 *
 * Fixed DOF: { 0 }
 *
 * Initial values
 *   2: 0.9
 *   3: 0.5
 * tRef = 0.4
 * Rate of change of values
 *   2: -0.4
 *   3: 0.2
 */

#include "DirichletDataMultiTri3.hh"

const int pylith::bc::DirichletDataMultiTri3::_numDOF = 2;

const int pylith::bc::DirichletDataMultiTri3::_idA = 0;
const char* pylith::bc::DirichletDataMultiTri3::_labelA = "bc";
const int pylith::bc::DirichletDataMultiTri3::_numFixedDOFA = 1;
const int pylith::bc::DirichletDataMultiTri3::_fixedDOFA[] = { 1 };
const int pylith::bc::DirichletDataMultiTri3::_numConstrainedPtsA = 2;
const int pylith::bc::DirichletDataMultiTri3::_constrainedPointsA[] = { 1, 3 };

const char* pylith::bc::DirichletDataMultiTri3::_dbFilenameA =
  "data/tri3.spatialdb";
const char* pylith::bc::DirichletDataMultiTri3::_dbFilenameARate =
  "data/tri3_rate.spatialdb";
const double pylith::bc::DirichletDataMultiTri3::_tRefA = 3.2;

const int pylith::bc::DirichletDataMultiTri3::_idB = 1;
const char* pylith::bc::DirichletDataMultiTri3::_labelB = "bc2";
const int pylith::bc::DirichletDataMultiTri3::_numFixedDOFB = 1;
const int pylith::bc::DirichletDataMultiTri3::_fixedDOFB[] = { 0 };
const int pylith::bc::DirichletDataMultiTri3::_numConstrainedPtsB = 2;
const int pylith::bc::DirichletDataMultiTri3::_constrainedPointsB[] = { 2, 3 };

const char* pylith::bc::DirichletDataMultiTri3::_dbFilenameB =
  "data/tri3_b.spatialdb";
const char* pylith::bc::DirichletDataMultiTri3::_dbFilenameBRate =
  "data/tri3_b_rate.spatialdb";
const double pylith::bc::DirichletDataMultiTri3::_tRefB = 0.4;

const int pylith::bc::DirichletDataMultiTri3::_constraintSizes[] = {
  0,
  1,
  1,
  2
};

const int pylith::bc::DirichletDataMultiTri3::_constrainedDOF[] = {
  1, 
  0,
  0, 1
};

// Values at t=10.0
const double pylith::bc::DirichletDataMultiTri3::_field[] = {
  0.0, 0.0,
  0.0, 1.66,
  -2.94, 0.0,
  2.42, 6.14
};

// Values from t=10.0 to t=14.0.
const double pylith::bc::DirichletDataMultiTri3::_fieldIncr[] = {
  0.0, 0.0,
  0.0, 0.8,
 -1.6, 0.0,
  0.8, 3.2
};

const char* pylith::bc::DirichletDataMultiTri3::_meshFilename = 
  "data/tri3.mesh";

pylith::bc::DirichletDataMultiTri3::DirichletDataMultiTri3(void)
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

  idB = _idA;
  labelB = const_cast<char*>(_labelB);
  numFixedDOFB = _numFixedDOFB;
  fixedDOFB = const_cast<int*>(_fixedDOFB);
  numConstrainedPtsB = _numConstrainedPtsB;
  constrainedPointsB = const_cast<int*>(_constrainedPointsB);

  dbFilenameB = const_cast<char*>(_dbFilenameB);
  dbFilenameBRate = const_cast<char*>(_dbFilenameBRate);
  tRefB = _tRefB;

  field = const_cast<double*>(_field);
  fieldIncr = const_cast<double*>(_fieldIncr);
  constraintSizes = const_cast<int*>(_constraintSizes);
  constrainedDOF = const_cast<int*>(_constrainedDOF);

  meshFilename = const_cast<char*>(_meshFilename);
} // constructor

pylith::bc::DirichletDataMultiTri3::~DirichletDataMultiTri3(void)
{}


// End of file
