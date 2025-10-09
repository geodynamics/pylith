// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

#include "CohesiveKinData.hh"

// ----------------------------------------------------------------------
// Constructor
pylith::faults::CohesiveKinData::CohesiveKinData(void) :
    meshFilename(0),
    setLengthScale(1.0),
    setPressureScale(2.0e+6),
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
    edge(0),
    finalSlipFilename(0),
    slipTimeFilename(0),
    riseTimeFilename(0),
    fieldT(0),
    fieldIncr(0),
    jacobianLumped(0),
    orientation(0),
    area(0),
    residual(0),
    jacobian(0),
    fieldIncrAdjusted(0),
    verticesFault(0),
    edgesLagrange(0),
    verticesPositive(0),
    verticesNegative(0),
    numFaultVertices(0),
    numCohesiveCells(0),
    cellMappingFault(0),
    cellMappingCohesive(0) { // constructor
    const PylithScalar velScale = lengthScale / timeScale;
    densityScale = pressureScale / (velScale*velScale);
} // constructor


// ----------------------------------------------------------------------
// Destructor
pylith::faults::CohesiveKinData::~CohesiveKinData(void) { // destructor
} // destructor


// End of file
