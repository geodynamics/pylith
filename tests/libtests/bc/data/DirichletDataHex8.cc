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

/* Mesh: meshHex8.txt
 *
 * Dirichlet BC at vertices 0, 1, 6, 7.
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

#include "DirichletDataHex8.hh"

const int pylith::bc::DirichletDataHex8::_id = 0;

const char* pylith::bc::DirichletDataHex8::_label = "bc";

const int pylith::bc::DirichletDataHex8::_numDOF = 3;
const int pylith::bc::DirichletDataHex8::_numFixedDOF = 2;
const int pylith::bc::DirichletDataHex8::_fixedDOF[] = { 0, 2 };

const int pylith::bc::DirichletDataHex8::_numConstrainedPts = 4;
const int pylith::bc::DirichletDataHex8::_constrainedPoints[] = { 0, 1, 6, 7 };

const PylithScalar pylith::bc::DirichletDataHex8::_tRef = 0.2;
const PylithScalar pylith::bc::DirichletDataHex8::_valueRate = 0.4;
const PylithScalar pylith::bc::DirichletDataHex8::_valuesInitial[] = {
  -0.2, 0.3,
   0.1, 0.7,
   0.5, 0.4,
   3.2, 6.1,
};

const char* pylith::bc::DirichletDataHex8::_meshFilename = 
  "data/hex8.mesh";
const char* pylith::bc::DirichletDataHex8::_dbFilename =
  "data/hex8_disp.spatialdb";

pylith::bc::DirichletDataHex8::DirichletDataHex8(void)
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
  valuesInitial = const_cast<PylithScalar*>(_valuesInitial);

  meshFilename = const_cast<char*>(_meshFilename);
  dbFilename = const_cast<char*>(_dbFilename);
} // constructor

pylith::bc::DirichletDataHex8::~DirichletDataHex8(void)
{}


// End of file
