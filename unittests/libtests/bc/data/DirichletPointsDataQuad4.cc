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

/* Mesh: meshQuad4.txt
 *
 * DirichletPoints BC at vertices 0, 1, 4.
 *
 * Fixed DOF: { 0, 1 }
 *
 * Values
 *   0: 0.1, 0.6
 *   1: 0.5, 0.3
 *   4: 0.4, 0.2
 * tRef = 3.0
 * Rate of change = -0.5
 */

#include "DirichletPointsDataQuad4.hh"

const int pylith::bc::DirichletPointsDataQuad4::_id = 0;

const char* pylith::bc::DirichletPointsDataQuad4::_label = "bc";

const int pylith::bc::DirichletPointsDataQuad4::_numDOF = 2;
const int pylith::bc::DirichletPointsDataQuad4::_numFixedDOF = 2;
const int pylith::bc::DirichletPointsDataQuad4::_fixedDOF[] = { 0, 1 };

const int pylith::bc::DirichletPointsDataQuad4::_numConstrainedPts = 3;
const int pylith::bc::DirichletPointsDataQuad4::_constrainedPoints[] = { 0, 1, 4 };

const double pylith::bc::DirichletPointsDataQuad4::_tRef = 3.0;
const double pylith::bc::DirichletPointsDataQuad4::_valueRate = -0.5;
const double pylith::bc::DirichletPointsDataQuad4::_valuesInitial[] =
  { 0.1, 0.6, 0.5, 0.3, 0.4, 0.2 };

const char* pylith::bc::DirichletPointsDataQuad4::_meshFilename = 
  "data/quad4.mesh";
const char* pylith::bc::DirichletPointsDataQuad4::_dbFilename =
  "data/quad4.spatialdb";

pylith::bc::DirichletPointsDataQuad4::DirichletPointsDataQuad4(void)
{ // constructor
  id = _id;
  label = const_cast<char*>(_label);

  numDOF = _numDOF;
  numFixedDOF = _numFixedDOF;
  fixedDOF = const_cast<int*>(_fixedDOF);

  numConstrainedPts = _numConstrainedPts;
  constrainedPoints = const_cast<int*>(_constrainedPoints);

  tRef = _tRef;
  valueRate = _valueRate;
  valuesInitial = const_cast<double*>(_valuesInitial);

  meshFilename = const_cast<char*>(_meshFilename);
  dbFilename = const_cast<char*>(_dbFilename);
} // constructor

pylith::bc::DirichletPointsDataQuad4::~DirichletPointsDataQuad4(void)
{}


// End of file
