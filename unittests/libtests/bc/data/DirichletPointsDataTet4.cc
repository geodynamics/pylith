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

/* Mesh: meshTet4.txt
 *
 * DirichletPoints BC at vertices 2.
 *
 * Fixed DOF: { 1, 2 }
 *
 * Values
 *   2: 0.7, 0.2
 * tRef = 1.2
 * Rate of change = 4.0
 */

#include "DirichletPointsDataTet4.hh"

const int pylith::bc::DirichletPointsDataTet4::_id = 0;

const char* pylith::bc::DirichletPointsDataTet4::_label = "bc";

const int pylith::bc::DirichletPointsDataTet4::_numDOF = 3;
const int pylith::bc::DirichletPointsDataTet4::_numFixedDOF = 2;
const int pylith::bc::DirichletPointsDataTet4::_fixedDOF[] = { 1, 2 };

const int pylith::bc::DirichletPointsDataTet4::_numConstrainedPts = 1;
const int pylith::bc::DirichletPointsDataTet4::_constrainedPoints[] = { 2 };


const double pylith::bc::DirichletPointsDataTet4::_tRef = 1.2;
const double pylith::bc::DirichletPointsDataTet4::_valueRate = 4.0;
const double pylith::bc::DirichletPointsDataTet4::_valuesInitial[] =
  { 0.7, 0.2 };

const char* pylith::bc::DirichletPointsDataTet4::_meshFilename = 
  "data/tet4.mesh";
const char* pylith::bc::DirichletPointsDataTet4::_dbFilename =
  "data/tet4.spatialdb";

pylith::bc::DirichletPointsDataTet4::DirichletPointsDataTet4(void)
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

pylith::bc::DirichletPointsDataTet4::~DirichletPointsDataTet4(void)
{}


// End of file
