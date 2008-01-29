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
 * DirichletPoints BC B at vertex 0 and 3.
 *
 * Fixed DOF: { 0 }
 *
 * Initial values
 *   0: 0.9
 *   3: 0.5
 * tRef = 0.4
 * Rate of change of values
 *   1: -0.4
 *   3: 0.2
 */

#include "DirichletPointsDataMultiTri3.hh"

const int pylith::bc::DirichletPointsDataMultiTri3::_numDOF = 2;

const int pylith::bc::DirichletPointsDataMultiTri3::_idA = 0;
const char* pylith::bc::DirichletPointsDataMultiTri3::_labelA = "bc";
const int pylith::bc::DirichletPointsDataMultiTri3::_numFixedDOFA = 1;
const int pylith::bc::DirichletPointsDataMultiTri3::_fixedDOFA[] = { 1 };
const int pylith::bc::DirichletPointsDataMultiTri3::_numConstrainedPtsA = 2;
const int pylith::bc::DirichletPointsDataMultiTri3::_constrainedPointsA[] = { 1, 3 };

const char* pylith::bc::DirichletPointsDataMultiTri3::_dbFilenameA =
  "data/tri3.spatialdb";
const char* pylith::bc::DirichletPointsDataMultiTri3::_dbFilenameARate =
  "data/tri3_rate.spatialdb";
const double pylith::bc::DirichletPointsDataMultiTri3::_tRefA = 3.2;

const int pylith::bc::DirichletPointsDataMultiTri3::_idB = 1;
const char* pylith::bc::DirichletPointsDataMultiTri3::_labelB = "bc2";
const int pylith::bc::DirichletPointsDataMultiTri3::_numFixedDOFB = 1;
const int pylith::bc::DirichletPointsDataMultiTri3::_fixedDOFB[] = { 0 };
const int pylith::bc::DirichletPointsDataMultiTri3::_numConstrainedPtsB = 1;
const int pylith::bc::DirichletPointsDataMultiTri3::_constrainedPointsB[] = { 0, 3 };

const char* pylith::bc::DirichletPointsDataMultiTri3::_dbFilenameB =
  "data/tri3_b.spatialdb";
const char* pylith::bc::DirichletPointsDataMultiTri3::_dbFilenameBRate =
  "data/tri3_b_rate.spatialdb";
const double pylith::bc::DirichletPointsDataMultiTri3::_tRefB = 0.4;

const int pylith::bc::DirichletPointsDataMultiTri3::_constraintSizes[] = {
  1,
  1,
  0,
  2
};

const int pylith::bc::DirichletPointsDataMultiTri3::_constrainedDOF[] = {
  0, 
  1,
  0, 1
};

// Values at t=10.0
const double pylith::bc::DirichletPointsDataMultiTri3::_field[] = {
  -2.94, 0.0,
  0.0, 1.66,
  0.0, 0.0,
  2.42, 6.14
};

const char* pylith::bc::DirichletPointsDataMultiTri3::_meshFilename = 
  "data/tri3.mesh";

pylith::bc::DirichletPointsDataMultiTri3::DirichletPointsDataMultiTri3(void)
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
  constraintSizes = const_cast<int*>(_constraintSizes);
  constrainedDOF = const_cast<int*>(_constrainedDOF);

  meshFilename = const_cast<char*>(_meshFilename);
} // constructor

pylith::bc::DirichletPointsDataMultiTri3::~DirichletPointsDataMultiTri3(void)
{}


// End of file
