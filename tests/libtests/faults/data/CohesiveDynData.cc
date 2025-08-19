// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

#include "CohesiveDynData.hh"

// ----------------------------------------------------------------------
// Constructor
pylith::faults::CohesiveDynData::CohesiveDynData(void) :
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
    initialTractFilename(0),
    fieldT(0),
    fieldIncrStick(0),
    fieldIncrSlip(0),
    fieldIncrOpen(0),
    jacobian(0),
    orientation(0),
    initialTractions(0),
    area(0),
    slipStickE(0),
    fieldIncrSlipE(0),
    slipSlipE(0),
    fieldIncrOpenE(0),
    slipOpenE(0),
    constraintEdges(0),
    numConstraintEdges(0) { // constructor
    const PylithScalar velScale = lengthScale / timeScale;
    densityScale = pressureScale / (velScale*velScale);
} // constructor


// ----------------------------------------------------------------------
// Destructor
pylith::faults::CohesiveDynData::~CohesiveDynData(void) { // destructor
} // destructor


// End of file
