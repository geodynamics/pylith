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
 * Fixed DOF: None
 */

#include "DirichletPointsDataLine2b.hh"

const int pylith::bc::DirichletPointsDataLine2b::_id = 0;

const char* pylith::bc::DirichletPointsDataLine2b::_label = "bc0";

const int pylith::bc::DirichletPointsDataLine2b::_numDOF = 1;
const int pylith::bc::DirichletPointsDataLine2b::_numFixedDOF = 0;
const int pylith::bc::DirichletPointsDataLine2b::_fixedDOF[] = {};

const int pylith::bc::DirichletPointsDataLine2b::_numConstrainedPts = 2;
const int pylith::bc::DirichletPointsDataLine2b::_constrainedPoints[] = { 0, 2 };
const double pylith::bc::DirichletPointsDataLine2b::_values[] = {0};

const char* pylith::bc::DirichletPointsDataLine2b::_meshFilename = 
  "data/line2.mesh";
const char* pylith::bc::DirichletPointsDataLine2b::_dbFilename =
  "data/line2.spatialdb";

pylith::bc::DirichletPointsDataLine2b::DirichletPointsDataLine2b(void)
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

pylith::bc::DirichletPointsDataLine2b::~DirichletPointsDataLine2b(void)
{}


// End of file
