// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2023, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================
#pragma once

#include "pylith/utils/types.hh" // USES PylithScalar

namespace pylith {
    namespace meshio {
        class MeshData;
    } // pylith
} // meshio

class pylith::meshio::MeshData {
    // PUBLIC METHODS ///////////////////////////////////////////////////////
public:

    /// Constructor
    MeshData(void);

    /// Destructor
    ~MeshData(void);

    // PUBLIC MEMBERS ///////////////////////////////////////////////////////
public:

    int numVertices; ///< Number of vertices
    int spaceDim; ///< Number of dimensions in vertex coordinates
    int numCells; ///< Number of cells
    int cellDim; ///< Number of dimensions associated with cell
    int numCorners; ///< Number of vertices in cell

    PylithScalar* vertices; ///< Pointer to coordinates of vertices
    int* cells; ///< Pointer to indices of vertices in cells
    int* materialIds; ///< Pointer to cell material identifiers

    int* groups; ///< Array of pointers to indices of points in groups
    int* groupSizes; ///< Array of sizes of each group
    char** groupNames; ///< Array of group names
    char** groupTypes; ///< Array of group types
    int numGroups; ///< Number of groups

    bool useIndexZero; ///< Indices start with 0 if true, 1 if false

};

// End of file
