// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

#include "CohesiveImpulsesData.hh"

// ----------------------------------------------------------------------
// Constructor
pylith::faults::CohesiveImpulsesData::CohesiveImpulsesData(void) :
    meshFilename(0),
    setLengthScale(1.0),
    setRigidityScale(2.5e+6),
    setTimeScale(2.0),
    spaceDim(0),
    cellDim(0),
    numBasis(0),
    numQuadPts(0),
    quadPts(0),
    quadWts(0),
    basis(0),
    basisDeriv(0),
    verticesRef(0),
    id(0),
    label(0),
    impulseAmpFilename(0),
    impulseDOF(0),
    numComponents(0),
    fieldT(0),
    fieldIncr(0),
    orientation(0),
    area(0),
    amplitude(0),
    numImpulses(0),
    residual(0),
    constraintEdges(0),
    numConstraintEdges(0) { // constructor
    const PylithScalar velScale = lengthScale / timeScale;
    densityScale = rigidityScale / (velScale*velScale);
} // constructor


// ----------------------------------------------------------------------
// Destructor
pylith::faults::CohesiveImpulsesData::~CohesiveImpulsesData(void) { // destructor
} // destructor


// End of file
