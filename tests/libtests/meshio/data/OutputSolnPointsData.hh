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

#include "pylith/utils/types.hh" // HASA PylithScalar

namespace pylith {
    namespace meshio {
        class OutputSolnPointsData;
    } // pylith
} // meshio

class pylith::meshio::OutputSolnPointsData {
    // PUBLIC METHODS ///////////////////////////////////////////////////////
public:

    /// Constructor
    OutputSolnPointsData(void);

    /// Destructor
    ~OutputSolnPointsData(void);

    // PUBLIC MEMBERS ///////////////////////////////////////////////////////
public:

    char* meshFilename;   ///< Filename for input mesh

    /// @name Point information.
    //@{
    int spaceDim;   ///< Number of dimensions in vertex coordinates
    int numPoints;   ///< Number of points.
    PylithScalar* points;   ///< Coordinates of points [numPoints*spaceDim].
    const char** names;   ///< Names of points (e.g., station names).
    //@}

    /// @name Field data.
    //@{
    int fiberDim;
    PylithScalar* coefs;   ///< Polynomial coefficients [fiberDim*spaceDim] for computing field.
    //@}

};

// End of file
