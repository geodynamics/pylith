// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

#include "AbsorbingDampersData.hh"

// ----------------------------------------------------------------------
// Constructor
pylith::bc::AbsorbingDampersData::AbsorbingDampersData(void) :
    meshFilename(0),
    setLengthScale(1.0),
    setPressureScale(2.0e+6),
    setTimeScale(2.0),
    numBasis(0),
    numQuadPts(0),
    quadPts(0),
    quadWts(0),
    basis(0),
    basisDerivRef(0),
    spatialDBFilename(0),
    id(0),
    label(""),
    dt(0),
    fieldTIncr(0),
    fieldT(0),
    fieldTmdt(0),
    spaceDim(0),
    cellDim(0),
    numVertices(0),
    numCells(0),
    numCorners(0),
    cells(0),
    dampingConsts(0),
    valsResidual(0),
    valsJacobian(0) { // constructor
    const PylithScalar velScale = lengthScale / timeScale;
    densityScale = pressureScale / (velScale*velScale);
} // constructor


// ----------------------------------------------------------------------
// Destructor
pylith::bc::AbsorbingDampersData::~AbsorbingDampersData(void) { // destructor
} // destructor


// End of file
