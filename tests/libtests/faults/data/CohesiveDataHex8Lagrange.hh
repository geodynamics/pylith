// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2024, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================
#pragma once

#include "CohesiveData.hh"

namespace pylith {
    namespace faults {
        class CohesiveDataHex8Lagrange;
    } // pylith
} // faults

class pylith::faults::CohesiveDataHex8Lagrange : public CohesiveData {
    // PUBLIC METHODS ///////////////////////////////////////////////////////
public:

    /// Constructor
    CohesiveDataHex8Lagrange(void);

    /// Destructor
    ~CohesiveDataHex8Lagrange(void);

    // PRIVATE MEMBERS //////////////////////////////////////////////////////
private:

    static const int _numVertices; ///< Number of vertices
    static const int _spaceDim; ///< Number of dimensions in vertex coordinates
    static const int _numCells; ///< Number of cells
    static const int _cellDim; ///< Number of dimensions associated with cell

    static const int _numCorners[]; ///< Number of vertices in cell
    static const int _materialIds[]; ///< Pointer to cell material identifiers

    static const int _groupSizes[]; ///< Sizes of groups
    static const char* _groupNames[]; ///< Array of group names
    static const char* _groupTypes[]; ///< Array of group types
    static const int _numGroups; ///< Number of groups

    static const char* _filename; ///< Filename of input mesh
};

// End of file
