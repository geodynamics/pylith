// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information. 
// =================================================================================================
 Values
 *   1: 0.3
 *   3: 0.7
 * tRef = 0.7
 * Rate of change = -0.2
 */

#include "PointForceDataTri3.hh"

const int pylith::bc::PointForceDataTri3::_id = 0;

const char* pylith::bc::PointForceDataTri3::_label = "bc";

const int pylith::bc::PointForceDataTri3::_numDOF = 2;
const int pylith::bc::PointForceDataTri3::_numForceDOF = 1;
const int pylith::bc::PointForceDataTri3::_forceDOF[] = { 1 };

const int pylith::bc::PointForceDataTri3::_numForcePts = 2;
const int pylith::bc::PointForceDataTri3::_forcePoints[] = { 1, 3 };

const PylithScalar pylith::bc::PointForceDataTri3::_tRef = 0.7;
const PylithScalar pylith::bc::PointForceDataTri3::_forceRate = -0.2;
const PylithScalar pylith::bc::PointForceDataTri3::_forceInitial[] =
{ 0.3, 0.7 };

const PylithScalar pylith::bc::PointForceDataTri3::_tResidual = 1.5;
const PylithScalar pylith::bc::PointForceDataTri3::_residual[] =
{ 0.0, 0.0,
  0.0, 0.14,
  0.0, 0.0,
  0.0, 0.54,
};

const char* pylith::bc::PointForceDataTri3::_meshFilename =
    "data/tri3.mesh";
const char* pylith::bc::PointForceDataTri3::_dbFilename =
    "data/tri3_force.spatialdb";

pylith::bc::PointForceDataTri3::PointForceDataTri3(void) { // constructor
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


pylith::bc::PointForceDataTri3::~PointForceDataTri3(void) {}


// End of file
