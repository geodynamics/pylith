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
 * DirichletPoints BC at vertices 0, 1, 6, 7.
 *
 * Fixed DOF: { 0, 2 }
 *
 * Initial values
 *   0: -0.2, 0.3
 *   1:  0.1, 0.7
 *   6:  0.5, 0.4
 *   7:  3.2, 6.1
 * tref = 0.2
 * Rate of change
 *   +0.4
 */

#include "DirichletPointsDataHex8.hh"

const int pylith::bc::DirichletPointsDataHex8::_id = 0;

const char* pylith::bc::DirichletPointsDataHex8::_label = "bc";

const int pylith::bc::DirichletPointsDataHex8::_numDOF = 3;
const int pylith::bc::DirichletPointsDataHex8::_numFixedDOF = 2;
const int pylith::bc::DirichletPointsDataHex8::_fixedDOF[] = { 0, 2 };

const int pylith::bc::DirichletPointsDataHex8::_numConstrainedPts = 4;
const int pylith::bc::DirichletPointsDataHex8::_constrainedPoints[] = { 0, 1, 6, 7 };

const double pylith::bc::DirichletPointsDataHex8::_tRef = 0.2;
const double pylith::bc::DirichletPointsDataHex8::_valueRate = 0.4;
const double pylith::bc::DirichletPointsDataHex8::_valuesInitial[] = {
  -0.2, 0.3,
   0.1, 0.7,
   0.5, 0.4,
   3.2, 6.1,
};

const char* pylith::bc::DirichletPointsDataHex8::_meshFilename = 
  "data/hex8.mesh";
const char* pylith::bc::DirichletPointsDataHex8::_dbFilename =
  "data/hex8.spatialdb";

pylith::bc::DirichletPointsDataHex8::DirichletPointsDataHex8(void)
{ // constructor
  id = _id;
  label = const_cast<char*>(_label);

  tRef = _tRef;
  valueRate = _valueRate;

  numDOF = _numDOF;
  numFixedDOF = _numFixedDOF;
  fixedDOF = const_cast<int*>(_fixedDOF);

  numConstrainedPts = _numConstrainedPts;
  constrainedPoints = const_cast<int*>(_constrainedPoints);
  valuesInitial = const_cast<double*>(_valuesInitial);

  meshFilename = const_cast<char*>(_meshFilename);
  dbFilename = const_cast<char*>(_dbFilename);
} // constructor

pylith::bc::DirichletPointsDataHex8::~DirichletPointsDataHex8(void)
{}


// End of file
