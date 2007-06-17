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
 * Dirichlet BC A at vertices 1 and 3.
 *
 * Fixed DOF: { 1 }
 *
 * Values
 *   1: 0.3
 *   3: 0.7
 *
 * Dirichlet BC B at vertex 0 and 3.
 *
 * Fixed DOF: { 0 }
 *
 * Values
 *   0: 0.9
 *   3: 0.5
 */

#include "DirichletDataMultiTri3.hh"

const int pylith::bc::DirichletDataMultiTri3::_numDOF = 2;

const int pylith::bc::DirichletDataMultiTri3::_idA = 0;
const char* pylith::bc::DirichletDataMultiTri3::_labelA = "bc";
const int pylith::bc::DirichletDataMultiTri3::_numFixedDOFA = 1;
const int pylith::bc::DirichletDataMultiTri3::_fixedDOFA[] = { 1 };
const int pylith::bc::DirichletDataMultiTri3::_numConstrainedPtsA = 2;
const int pylith::bc::DirichletDataMultiTri3::_constrainedPointsA[] = { 1, 3 };
const double pylith::bc::DirichletDataMultiTri3::_valuesA[] = { 0.3, 0.7 };
const char* pylith::bc::DirichletDataMultiTri3::_dbFilenameA =
  "data/tri3.spatialdb";

const int pylith::bc::DirichletDataMultiTri3::_idB = 1;
const char* pylith::bc::DirichletDataMultiTri3::_labelB = "bc2";
const int pylith::bc::DirichletDataMultiTri3::_numFixedDOFB = 1;
const int pylith::bc::DirichletDataMultiTri3::_fixedDOFB[] = { 0 };
const int pylith::bc::DirichletDataMultiTri3::_numConstrainedPtsB = 1;
const int pylith::bc::DirichletDataMultiTri3::_constrainedPointsB[] = { 0, 3 };
const double pylith::bc::DirichletDataMultiTri3::_valuesB[] = { 0.9, 0.5 };
const char* pylith::bc::DirichletDataMultiTri3::_dbFilenameB =
  "data/tri3_b.spatialdb";

const int pylith::bc::DirichletDataMultiTri3::_constraintSizes[] = {
  1,
  1,
  0,
  2
};

const int pylith::bc::DirichletDataMultiTri3::_constrainedDOF[] = {
  0, 
  1,
  0, 1
};

const double pylith::bc::DirichletDataMultiTri3::_field[] = {
  0.9, 0.0,
  0.0, 0.3,
  0.0, 0.0,
  0.5, 0.7
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
  valuesA = const_cast<double*>(_valuesA);
  dbFilenameA = const_cast<char*>(_dbFilenameA);

  idB = _idA;
  labelB = const_cast<char*>(_labelB);
  numFixedDOFB = _numFixedDOFB;
  fixedDOFB = const_cast<int*>(_fixedDOFB);
  numConstrainedPtsB = _numConstrainedPtsB;
  constrainedPointsB = const_cast<int*>(_constrainedPointsB);
  valuesB = const_cast<double*>(_valuesB);
  dbFilenameB = const_cast<char*>(_dbFilenameB);

  field = const_cast<double*>(_field);
  constraintSizes = const_cast<int*>(_constraintSizes);
  constrainedDOF = const_cast<int*>(_constrainedDOF);

  meshFilename = const_cast<char*>(_meshFilename);
} // constructor

pylith::bc::DirichletDataMultiTri3::~DirichletDataMultiTri3(void)
{}


// End of file
