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

/* Mesh: meshHex8.txt
 *
 * Dirichlet BC at vertices 0, 1, 6, 7.
 *
 * Fixed DOF: { 0, 2 }
 *
 * Values
 *   0: -0.2, 0.3
 *   1:  0.1, 0.7
 *   6:  0.5, 0.4
 *   7:  3.2, 6.1
 */

#include "DirichletDataHex8.hh"

const int pylith::bc::DirichletDataHex8::_id = 0;

const char* pylith::bc::DirichletDataHex8::_label = "bc";

const int pylith::bc::DirichletDataHex8::_numDOF = 3;
const int pylith::bc::DirichletDataHex8::_numFixedDOF = 2;
const int pylith::bc::DirichletDataHex8::_fixedDOF[] = { 0, 2 };

const int pylith::bc::DirichletDataHex8::_numConstrainedPts = 4;
const int pylith::bc::DirichletDataHex8::_constrainedPoints[] = { 0, 1, 6, 7 };
const double pylith::bc::DirichletDataHex8::_values[] = {
  -0.2, 0.3,
   0.1, 0.7,
   0.5, 0.4,
   3.2, 6.1,
};

const char* pylith::bc::DirichletDataHex8::_meshFilename = 
  "data/hex8.mesh";
const char* pylith::bc::DirichletDataHex8::_dbFilename =
  "data/hex8.spatialdb";

pylith::bc::DirichletDataHex8::DirichletDataHex8(void)
{ // constructor
  id = _id;
  label = const_cast<char*>(_label);

  numDOF = _numDOF;
  numFixedDOF = _numFixedDOF;
  fixedDOF = const_cast<int*>(_fixedDOF);

  numConstrainedPts = _numConstrainedPts;
  constrainedPoints = const_cast<int*>(_constrainedPoints);
  values = const_cast<double*>(_values);

  meshFilename = const_cast<char*>(_meshFilename);
  dbFilename = const_cast<char*>(_dbFilename);
} // constructor

pylith::bc::DirichletDataHex8::~DirichletDataHex8(void)
{}


// End of file
