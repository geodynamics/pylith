// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

#include "OutputSolnPointsData.hh"

// ----------------------------------------------------------------------
// Constructor
pylith::meshio::OutputSolnPointsData::OutputSolnPointsData(void) :
    meshFilename(0),
    spaceDim(0),
    numPoints(0),
    points(0),
    names(0),
    fiberDim(0),
    coefs(0) { // constructor
} // constructor


// ----------------------------------------------------------------------
// Destructor
pylith::meshio::OutputSolnPointsData::~OutputSolnPointsData(void) { // destructor
} // destructor


// End of file
