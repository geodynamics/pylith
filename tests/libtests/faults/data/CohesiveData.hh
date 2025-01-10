// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================
#pragma once

#include "pylith/utils/types.hh" // HASA PylithScalar

namespace pylith {
    namespace faults {
        class CohesiveData;
    } // pylith
} // faults

class pylith::faults::CohesiveData {
    // PUBLIC METHODS ///////////////////////////////////////////////////////
public:

    /// Constructor
    CohesiveData(void);

    /// Destructor
    ~CohesiveData(void);

    // PUBLIC MEMBERS ///////////////////////////////////////////////////////
public:

    int numVertices; ///< Number of vertices
    int spaceDim; ///< Number of dimensions in vertex coordinates
    int numCells; ///< Number of cells
    int cellDim; ///< Number of dimensions associated with cell

    int* numCorners; ///< Number of vertices in cell
    int* materialIds; ///< Pointer to cell material identifiers

    int* groupSizes; ///< Array of sizes of each group
    char** groupNames; ///< Array of group names
    char** groupTypes; ///< Array of group types
    int numGroups; ///< Number of groups

    char* filename; ///< Filename for input mesh
    char* fault; ///< Label for fault.
    char* edge; ///< Label for fault edge.
};

// End of file
