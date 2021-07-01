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

/* Mesh: meshTet4.txt
 *
 * Dirichlet BC at vertices 2.
 *
 * Fixed DOF: { 1, 2 }
 *
 * Values
 *   1: 0.7, 0.2
 *   2: 0.7, 0.2
 *   3: 0.7, 0.2
 * tRef = 1.2
 * Rate of change = 4.0
 */

#include "DirichletDataTet4.hh"

const int pylith::bc::DirichletDataTet4::_id = 0;

const char* pylith::bc::DirichletDataTet4::_label = "bc3";

const int pylith::bc::DirichletDataTet4::_numDOF = 3;
const int pylith::bc::DirichletDataTet4::_numFixedDOF = 2;
const int pylith::bc::DirichletDataTet4::_fixedDOF[] = { 1, 2 };

const int pylith::bc::DirichletDataTet4::_numConstrainedPts = 3;
const int pylith::bc::DirichletDataTet4::_constrainedPoints[] = { 1, 2, 3 };


const PylithScalar pylith::bc::DirichletDataTet4::_tRef = 1.2;
const PylithScalar pylith::bc::DirichletDataTet4::_valueRate = 4.0;
const PylithScalar pylith::bc::DirichletDataTet4::_valuesInitial[] = {
  0.7, 0.2,
  0.7, 0.2,
  0.7, 0.2,
};

const char* pylith::bc::DirichletDataTet4::_meshFilename = 
  "data/tet4.mesh";
const char* pylith::bc::DirichletDataTet4::_dbFilename =
  "data/tet4_disp.spatialdb";

pylith::bc::DirichletDataTet4::DirichletDataTet4(void)
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

pylith::bc::DirichletDataTet4::~DirichletDataTet4(void)
{}


// End of file
