// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2023, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

#include "NeumannData.hh"

// ----------------------------------------------------------------------
// Constructor
pylith::bc::NeumannData::NeumannData(void) :
    meshFilename(0),
    setLengthScale(1.0e+3),
    setPressureScale(2.25e+10),
    setDensityScale(1.0),
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
    spaceDim(0),
    cellDim(0),
    numVertices(0),
    numCells(0),
    numCorners(0),
    cells(0),
    tractionsCell(0),
    valsResidual(0) { // constructor
    const PylithScalar velScale = lengthScale / timeScale;
    densityScale = pressureScale / (velScale*velScale);
} // constructor


// ----------------------------------------------------------------------
// Destructor
pylith::bc::NeumannData::~NeumannData(void) { // destructor
} // destructor


// End of file
