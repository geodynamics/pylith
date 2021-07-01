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

/* Mesh: meshQuad4.txt
 *
 * Dirichlet BC at vertices 0, 2, 4.
 *
 * Fixed DOF: { 0, 1 }
 *
 * Values
 *   0: 0.1, 0.6
 *   2: 0.5, 0.3
 *   4: 0.4, 0.2
 * tRef = 3.0
 * Rate of change = -0.5
 */

#include "DirichletDataQuad4.hh"

const int pylith::bc::DirichletDataQuad4::_id = 0;

const char* pylith::bc::DirichletDataQuad4::_label = "bc3";

const int pylith::bc::DirichletDataQuad4::_numDOF = 2;
const int pylith::bc::DirichletDataQuad4::_numFixedDOF = 2;
const int pylith::bc::DirichletDataQuad4::_fixedDOF[] = { 0, 1 };

const int pylith::bc::DirichletDataQuad4::_numConstrainedPts = 3;
const int pylith::bc::DirichletDataQuad4::_constrainedPoints[] = { 0, 2, 4 };

const PylithScalar pylith::bc::DirichletDataQuad4::_tRef = 3.0;
const PylithScalar pylith::bc::DirichletDataQuad4::_valueRate = -0.5;
const PylithScalar pylith::bc::DirichletDataQuad4::_valuesInitial[] =
  { 0.1, 0.6, 0.5, 0.3, 0.4, 0.2 };

const char* pylith::bc::DirichletDataQuad4::_meshFilename = 
  "data/quad4.mesh";
const char* pylith::bc::DirichletDataQuad4::_dbFilename =
  "data/quad4_disp.spatialdb";

pylith::bc::DirichletDataQuad4::DirichletDataQuad4(void)
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

pylith::bc::DirichletDataQuad4::~DirichletDataQuad4(void)
{}


// End of file
