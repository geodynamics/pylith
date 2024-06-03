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

#include "pylith/utils/types.hh" // HASA PylithScalar

namespace pylith {
    namespace bc {
        class NeumannData;
    } // pylith
} // bc

class pylith::bc::NeumannData {
    // PUBLIC METHODS ///////////////////////////////////////////////////////
public:

    /// Constructor
    NeumannData(void);

    /// Destructor
    ~NeumannData(void);

    // PUBLIC MEMBERS ///////////////////////////////////////////////////////
public:

    const char* meshFilename; ///< Filename for input mesh.

    /// @name Scales information for nondimensionalization.
    //@{
    PylithScalar lengthScale; ///< Length scale.
    PylithScalar pressureScale; ///< Pressure scale.
    PylithScalar timeScale; ///< Time scale.
    PylithScalar densityScale; ///< Density scale.
    //@}

    /// @name Quadrature information
    //@{
    int numBasis; ///< Number of basis functions for surface cells.
    int numQuadPts; ///< Number of quadrature points per boundary cell.
    PylithScalar* quadPts; ///< Coordinates of quadrature points in ref cell.
    PylithScalar* quadWts; ///< Weights of quadrature points.
    PylithScalar* basis; ///< Cell basis functions at quad points.
    PylithScalar* basisDerivRef; ///< Derivatives of basis functions at quad points.
    //@}

    /// @name Parameter information
    //@{
    const char* spatialDBFilename; ///< Filename of simple spatial database.
    int id; ///< Boundary condition identifier
    const char* label; ///< Label for boundary condition group
    //@}

    /// @name Boundary mesh information
    //@{
    int spaceDim; ///< Number of dimensions in vertex coordinates
    int cellDim; ///< Number of dimensions associated with cell
    int numVertices; ///< Number of vertices
    int numCells; ///< Number of cells
    int numCorners; ///< Number of vertices in cell
    int* cells; ///< Indices of vertices in cells
    //@}

    /// @name Calculated values.
    //@{
    PylithScalar* tractionsCell; ///< Expected traction values at quadrature points.
    PylithScalar* valsResidual; ///< Expected residual at each vertex.
    //@}

};

// End of file
