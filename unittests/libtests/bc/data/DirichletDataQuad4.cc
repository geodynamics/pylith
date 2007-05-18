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
 * Dirichlet BC at vertices 0, 1, 4.
 *
 * Fixed DOF: { 0, 1 }
 *
 * Values
 *   0: 0.1, 0.6
 *   1: 0.5, 0.3
 *   4: 0.4, 0.2
 */

#include "DirichletDataQuad4.hh"

const int pylith::bc::DirichletDataQuad4::_id = 0;

const char* pylith::bc::DirichletDataQuad4::_label = "bc";

const int pylith::bc::DirichletDataQuad4::_numDOF = 2;
const int pylith::bc::DirichletDataQuad4::_numFixedDOF = 2;
const int pylith::bc::DirichletDataQuad4::_fixedDOF[] = { 0, 1 };

const int pylith::bc::DirichletDataQuad4::_numConstrainedPts = 3;
const int pylith::bc::DirichletDataQuad4::_constrainedPoints[] = { 0, 1, 4 };
const double pylith::bc::DirichletDataQuad4::_values[] =
  { 0.1, 0.6, 0.5, 0.3, 0.4, 0.2 };

const char* pylith::bc::DirichletDataQuad4::_meshFilename = 
  "data/quad4.mesh";
const char* pylith::bc::DirichletDataQuad4::_dbFilename =
  "data/quad4.spatialdb";

pylith::bc::DirichletDataQuad4::DirichletDataQuad4(void)
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

pylith::bc::DirichletDataQuad4::~DirichletDataQuad4(void)
{}


// End of file
