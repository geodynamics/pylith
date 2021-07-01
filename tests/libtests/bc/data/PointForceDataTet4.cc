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
 * Point force BC at vertices 2.
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

#include "PointForceDataTet4.hh"

const int pylith::bc::PointForceDataTet4::_id = 0;

const char* pylith::bc::PointForceDataTet4::_label = "bc3";

const int pylith::bc::PointForceDataTet4::_numDOF = 3;
const int pylith::bc::PointForceDataTet4::_numForceDOF = 2;
const int pylith::bc::PointForceDataTet4::_forceDOF[] = { 1, 2 };

const int pylith::bc::PointForceDataTet4::_numForcePts = 3;
const int pylith::bc::PointForceDataTet4::_forcePoints[] = { 1, 2, 3 };

const PylithScalar pylith::bc::PointForceDataTet4::_tRef = 1.2;
const PylithScalar pylith::bc::PointForceDataTet4::_forceRate = 4.0;
const PylithScalar pylith::bc::PointForceDataTet4::_forceInitial[] = {
  0.7, 0.2,
  0.7, 0.2,
  0.7, 0.2,
};

const PylithScalar pylith::bc::PointForceDataTet4::_tResidual = 1.3;
const PylithScalar pylith::bc::PointForceDataTet4::_residual[] =
  { 0.0, 0.0, 0.0,
    0.0, 1.1, 0.6,
    0.0, 1.1, 0.6,
    0.0, 1.1, 0.6,
    0.0, 0.0, 0.0,
  };

const char* pylith::bc::PointForceDataTet4::_meshFilename = 
  "data/tet4.mesh";
const char* pylith::bc::PointForceDataTet4::_dbFilename =
  "data/tet4_force.spatialdb";

pylith::bc::PointForceDataTet4::PointForceDataTet4(void)
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

pylith::bc::PointForceDataTet4::~PointForceDataTet4(void)
{}


// End of file
