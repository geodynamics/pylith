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
// Copyright (c) 2010-2013 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

/* Mesh: meshLine2.txt
 *
 * DirichletPoints BC at vertices 0 and 2.
 *
 * Fixed DOF: None
 */

#include "DirichletDataLine2b.hh"

const int pylith::bc::DirichletDataLine2b::_id = 0;

const char* pylith::bc::DirichletDataLine2b::_label = "bc0";

const int pylith::bc::DirichletDataLine2b::_numDOF = 1;
const int pylith::bc::DirichletDataLine2b::_numFixedDOF = 0;
const int pylith::bc::DirichletDataLine2b::_fixedDOF[] = {};

const int pylith::bc::DirichletDataLine2b::_numConstrainedPts = 2;
const int pylith::bc::DirichletDataLine2b::_constrainedPoints[] = { 0, 2 };

const PylithScalar pylith::bc::DirichletDataLine2b::_tRef = -0.2;
const PylithScalar pylith::bc::DirichletDataLine2b::_valueRate = 1.0;
const PylithScalar pylith::bc::DirichletDataLine2b::_valuesInitial[] = {0};

const char* pylith::bc::DirichletDataLine2b::_meshFilename = 
  "data/line2.mesh";
const char* pylith::bc::DirichletDataLine2b::_dbFilename =
  "data/line2_disp2.spatialdb";

pylith::bc::DirichletDataLine2b::DirichletDataLine2b(void)
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

pylith::bc::DirichletDataLine2b::~DirichletDataLine2b(void)
{}


// End of file
