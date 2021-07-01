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
 * Point force BC at vertices 0, 2, 4.
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

#include "PointForceDataQuad4.hh"

const int pylith::bc::PointForceDataQuad4::_id = 0;

const char* pylith::bc::PointForceDataQuad4::_label = "bc3";

const int pylith::bc::PointForceDataQuad4::_numDOF = 2;
const int pylith::bc::PointForceDataQuad4::_numForceDOF = 2;
const int pylith::bc::PointForceDataQuad4::_forceDOF[] = { 0, 1 };

const int pylith::bc::PointForceDataQuad4::_numForcePts = 3;
const int pylith::bc::PointForceDataQuad4::_forcePoints[] = { 0, 2, 4 };

const PylithScalar pylith::bc::PointForceDataQuad4::_tRef = 3.0;
const PylithScalar pylith::bc::PointForceDataQuad4::_forceRate = -0.5;
const PylithScalar pylith::bc::PointForceDataQuad4::_forceInitial[] =
  { 0.1, 0.6,
    0.5, 0.3,
    0.4, 0.2,
  };

const PylithScalar pylith::bc::PointForceDataQuad4::_tResidual = 4.5;
const PylithScalar pylith::bc::PointForceDataQuad4::_residual[] =
  { -0.65, -0.15,
     0.0, 0.0,
    -0.25, -0.45,
     0.0, 0.0,
    -0.35, -0.55,
     0.0, 0.0,
  };

const char* pylith::bc::PointForceDataQuad4::_meshFilename = 
  "data/quad4.mesh";
const char* pylith::bc::PointForceDataQuad4::_dbFilename =
  "data/quad4_force.spatialdb";

pylith::bc::PointForceDataQuad4::PointForceDataQuad4(void)
{ // constructor
  id = _id;
  label = const_cast<char*>(_label);

  numDOF = _numDOF;
  numForceDOF = _numForceDOF;
  forceDOF = const_cast<int*>(_forceDOF);

  numForcePts = _numForcePts;
  forcePoints = const_cast<int*>(_forcePoints);

  tRef = _tRef;
  forceRate = _forceRate;
  forceInitial = const_cast<PylithScalar*>(_forceInitial);

  tResidual = _tResidual;
  residual = const_cast<PylithScalar*>(_residual);

  meshFilename = const_cast<char*>(_meshFilename);
  dbFilename = const_cast<char*>(_dbFilename);
} // constructor

pylith::bc::PointForceDataQuad4::~PointForceDataQuad4(void)
{}


// End of file
