// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2023, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

#include "MeshData.hh"

// ----------------------------------------------------------------------
// Constructor
pylith::meshio::MeshData::MeshData(void) :
    numVertices(0),
    spaceDim(0),
    numCells(0),
    cellDim(0),
    numCorners(0),
    vertices(0),
    cells(0),
    materialIds(0),
    groups(0),
    groupSizes(0),
    groupNames(0),
    groupTypes(0),
    numGroups(0),
    useIndexZero(true) { // constructor
} // constructor


// ----------------------------------------------------------------------
// Destructor
pylith::meshio::MeshData::~MeshData(void) { // destructor
} // destructor


// End of file
