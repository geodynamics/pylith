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

/* Mesh: meshLine2.txt
 *
 * DirichletPoints BC at vertices 0 and 2.
 *
 * Fixed DOF: { 0 }
 *
 * Values
 *   0: 1.1 [constrained]
 *   1: 0.8 [solution]
 *   2: 2.2 [constrained]
 */

#include "DirichletPointsDataLine2.hh"

const int pylith::bc::DirichletPointsDataLine2::_id = 0;

const char* pylith::bc::DirichletPointsDataLine2::_label = "bc0";

const int pylith::bc::DirichletPointsDataLine2::_numDOF = 1;
const int pylith::bc::DirichletPointsDataLine2::_numFixedDOF = 1;
const int pylith::bc::DirichletPointsDataLine2::_fixedDOF[] = { 0 };

const int pylith::bc::DirichletPointsDataLine2::_numConstrainedPts = 2;
const int pylith::bc::DirichletPointsDataLine2::_constrainedPoints[] = { 0, 2 };
const double pylith::bc::DirichletPointsDataLine2::_values[] = { 1.1, 2.2 };

const char* pylith::bc::DirichletPointsDataLine2::_meshFilename = 
  "data/line2.mesh";
const char* pylith::bc::DirichletPointsDataLine2::_dbFilename =
  "data/line2.spatialdb";

pylith::bc::DirichletPointsDataLine2::DirichletPointsDataLine2(void)
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

pylith::bc::DirichletPointsDataLine2::~DirichletPointsDataLine2(void)
{}


// End of file
