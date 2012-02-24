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
 * Point force BC at vertices 0 and 2.
 *
 * Fixed DOF: { 0 }
 *
 * Values
 *   0: 1.1
 *   2: 2.2
 * tref = 0.6
 * Rate of change
 *   +0.3
 */

#include "PointForceDataLine2.hh"

const int pylith::bc::PointForceDataLine2::_id = 0;

const char* pylith::bc::PointForceDataLine2::_label = "bc0";

const int pylith::bc::PointForceDataLine2::_numDOF = 1;
const int pylith::bc::PointForceDataLine2::_numForceDOF = 1;
const int pylith::bc::PointForceDataLine2::_forceDOF[] = { 0 };

const int pylith::bc::PointForceDataLine2::_numForcePts = 2;
const int pylith::bc::PointForceDataLine2::_forcePoints[] = { 0, 2 };

const PylithScalar pylith::bc::PointForceDataLine2::_tRef = 0.6;
const PylithScalar pylith::bc::PointForceDataLine2::_forceRate = 0.3;
const PylithScalar pylith::bc::PointForceDataLine2::_forceInitial[] =
  { 1.1, 2.2 };

const PylithScalar pylith::bc::PointForceDataLine2::_tResidual = 1.5;
const PylithScalar pylith::bc::PointForceDataLine2::_residual[] =
  { 1.37,
    0.0,
    2.47,
  };

const char* pylith::bc::PointForceDataLine2::_meshFilename = 
  "data/line2.mesh";
const char* pylith::bc::PointForceDataLine2::_dbFilename =
  "data/line2_force.spatialdb";

pylith::bc::PointForceDataLine2::PointForceDataLine2(void)
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

pylith::bc::PointForceDataLine2::~PointForceDataLine2(void)
{}


// End of file
