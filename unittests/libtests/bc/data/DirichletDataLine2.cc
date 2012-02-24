// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2012 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

/* Mesh: meshLine2.txt
 *
 * Dirichlet BC at vertices 0 and 2.
 *
 * Fixed DOF: { 0 }
 *
 * Values
 *   0: 1.1 [constrained]
 *   1: 0.8 [solution]
 *   2: 2.2 [constrained]
 * tref = 0.6
 * Rate of change
 *   +0.3
 */

#include "DirichletDataLine2.hh"

const int pylith::bc::DirichletDataLine2::_id = 0;

const char* pylith::bc::DirichletDataLine2::_label = "bc0";

const int pylith::bc::DirichletDataLine2::_numDOF = 1;
const int pylith::bc::DirichletDataLine2::_numFixedDOF = 1;
const int pylith::bc::DirichletDataLine2::_fixedDOF[] = { 0 };

const int pylith::bc::DirichletDataLine2::_numConstrainedPts = 2;
const int pylith::bc::DirichletDataLine2::_constrainedPoints[] = { 0, 2 };

const PylithScalar pylith::bc::DirichletDataLine2::_tRef = 0.6;
const PylithScalar pylith::bc::DirichletDataLine2::_valueRate = 0.3;
const PylithScalar pylith::bc::DirichletDataLine2::_valuesInitial[] =
  { 1.1, 2.2 };

const char* pylith::bc::DirichletDataLine2::_meshFilename = 
  "data/line2.mesh";
const char* pylith::bc::DirichletDataLine2::_dbFilename =
  "data/line2_disp.spatialdb";

pylith::bc::DirichletDataLine2::DirichletDataLine2(void)
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

pylith::bc::DirichletDataLine2::~DirichletDataLine2(void)
{}


// End of file
