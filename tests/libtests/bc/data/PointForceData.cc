// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

#include "PointForceData.hh"

// ----------------------------------------------------------------------
// Constructor
pylith::bc::PointForceData::PointForceData(void) :
    tRef(0),
    forceRate(0),
    tResidual(0),
    numDOF(0),
    numForceDOF(0),
    numForcePts(0),
    id(0),
    label(0),
    forceDOF(0),
    forcePoints(0),
    forceInitial(0),
    residual(0),
    meshFilename(0),
    dbFilename(0),
    setLengthScale(1.0),
    setRigidityScale(2.0e+6),
    setTimeScale(2.0) { // constructor
    const PylithScalar velScale = lengthScale / timeScale;
    densityScale = rigidityScale / (velScale*velScale);
} // constructor


// ----------------------------------------------------------------------
// Destructor
pylith::bc::PointForceData::~PointForceData(void) { // destructor
} // destructor


// End of file
