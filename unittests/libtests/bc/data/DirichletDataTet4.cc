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
 * Dirichlet BC at vertices 2.
 *
 * Fixed DOF: { 1, 2 }
 *
 * Values
 *   2: 0.7, 0.2
 */

#include "DirichletDataTet4.hh"

const int pylith::bc::DirichletDataTet4::_id = 0;

const char* pylith::bc::DirichletDataTet4::_label = "bc";

const int pylith::bc::DirichletDataTet4::_numDOF = 3;
const int pylith::bc::DirichletDataTet4::_numFixedDOF = 2;
const int pylith::bc::DirichletDataTet4::_fixedDOF[] = { 1, 2 };

const int pylith::bc::DirichletDataTet4::_numConstrainedPts = 1;
const int pylith::bc::DirichletDataTet4::_constrainedPoints[] = { 2 };
const double pylith::bc::DirichletDataTet4::_values[] = { 0.7, 0.2 };

const char* pylith::bc::DirichletDataTet4::_meshFilename = 
  "data/tet4.mesh";
const char* pylith::bc::DirichletDataTet4::_dbFilename =
  "data/tet4.spatialdb";

pylith::bc::DirichletDataTet4::DirichletDataTet4(void)
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

pylith::bc::DirichletDataTet4::~DirichletDataTet4(void)
{}


// End of file
