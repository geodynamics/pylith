// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University at Buffalo
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2021 University of California, Davis
//
// See LICENSE.md for license information.
//
// ======================================================================
//

/* Mesh: meshTri3.txt
 *
 * Dirichlet BC at vertices 1 and 3.
 *
 * Fixed DOF: { 1 }
 *
 * Values
 *   1: 0.3
 *   3: 0.7
 * tRef = 0.7
 * Rate of change = -0.2
 */

#include "DirichletDataTri3.hh"

const int pylith::bc::DirichletDataTri3::_id = 0;

const char* pylith::bc::DirichletDataTri3::_label = "bc";

const int pylith::bc::DirichletDataTri3::_numDOF = 2;
const int pylith::bc::DirichletDataTri3::_numFixedDOF = 1;
const int pylith::bc::DirichletDataTri3::_fixedDOF[] = { 1 };

const int pylith::bc::DirichletDataTri3::_numConstrainedPts = 2;
const int pylith::bc::DirichletDataTri3::_constrainedPoints[] = { 1, 3 };

const PylithScalar pylith::bc::DirichletDataTri3::_tRef = 0.7;
const PylithScalar pylith::bc::DirichletDataTri3::_valueRate = -0.2;
const PylithScalar pylith::bc::DirichletDataTri3::_valuesInitial[] =
  { 0.3, 0.7 };

const char* pylith::bc::DirichletDataTri3::_meshFilename = 
  "data/tri3.mesh";
const char* pylith::bc::DirichletDataTri3::_dbFilename =
  "data/tri3_disp.spatialdb";

pylith::bc::DirichletDataTri3::DirichletDataTri3(void)
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
  valuesInitial = const_cast<PylithScalar*>(_valuesInitial);

  meshFilename = const_cast<char*>(_meshFilename);
  dbFilename = const_cast<char*>(_dbFilename);
} // constructor

pylith::bc::DirichletDataTri3::~DirichletDataTri3(void)
{}


// End of file
