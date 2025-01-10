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

#include "OutputSolnPointsData.hh"

namespace pylith {
    namespace meshio {
        class OutputSolnPointsDataTet4;
    } // pylith
} // meshio

class pylith::meshio::OutputSolnPointsDataTet4 : public OutputSolnPointsData {
    // PUBLIC METHODS ///////////////////////////////////////////////////////
public:

    /// Constructor
    OutputSolnPointsDataTet4(void);

    /// Destructor
    ~OutputSolnPointsDataTet4(void);

    // PRIVATE MEMBERS //////////////////////////////////////////////////////
private:

    static const char* _meshFilename; ///< Filename of input mesh

    static const int _spaceDim; ///< Number of dimensions in vertex coordinates
    static const int _numPoints; ///< Number of points.
    static const PylithScalar _points[]; ///< Coordinates of points.
    static const char* _names[]; /// Names of stations.

    static const int _fiberDim; ///< Fiber dimension of field.
    static const PylithScalar _coefs[]; ///< Polynomial coefficients for field.

};

// End of file
