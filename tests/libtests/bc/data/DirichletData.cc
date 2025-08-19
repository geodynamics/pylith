// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

#include "DirichletData.hh"

// ----------------------------------------------------------------------
// Constructor
pylith::bc::DirichletData::DirichletData(void) :
    tRef(0),
    valueRate(0),
    numDOF(0),
    numFixedDOF(0),
    numConstrainedPts(0),
    id(0),
    label(0),
    fixedDOF(0),
    constrainedPoints(0),
    valuesInitial(0),
    meshFilename(0),
    dbFilename(0),
    setLengthScale(1.0),
    setPressureScale(2.0e+6),
    setTimeScale(2.0) {
    const PylithScalar velScale = lengthScale / timeScale;
    densityScale = pressureScale / (velScale*velScale);
} // constructor


// ----------------------------------------------------------------------
// Destructor
pylith::bc::DirichletData::~DirichletData(void) { // destructor
} // destructor


// End of file
