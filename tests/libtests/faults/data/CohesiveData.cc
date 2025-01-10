// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

#include "CohesiveData.hh"

// ----------------------------------------------------------------------
// Constructor
pylith::faults::CohesiveData::CohesiveData(void) :
    numVertices(0),
    spaceDim(0),
    numCells(0),
    cellDim(0),
    numCorners(0),
    materialIds(0),
    groupSizes(0),
    groupNames(0),
    groupTypes(0),
    numGroups(0),
    filename(0),
    fault("fault"),
    edge(0) { // constructor
} // constructor


// ----------------------------------------------------------------------
// Destructor
pylith::faults::CohesiveData::~CohesiveData(void) { // destructor
} // destructor


// End of file
