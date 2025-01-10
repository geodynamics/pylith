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

#include "NeumannData.hh"

namespace pylith {
    namespace bc {
        class NeumannDataHex8;
    } // pylith
} // bc

class pylith::bc::NeumannDataHex8 : public NeumannData {
    // PUBLIC METHODS ///////////////////////////////////////////////////////
public:

    /// Constructor
    NeumannDataHex8(void);

    /// Destructor
    ~NeumannDataHex8(void);

    // PRIVATE MEMBERS //////////////////////////////////////////////////////
private:

    static const char* _meshFilename;

    // Quadrature information
    static const int _numBasis;
    static const int _numQuadPts;
    static const PylithScalar _quadPts[];
    static const PylithScalar _quadWts[];
    static const PylithScalar _basis[];
    static const PylithScalar _basisDerivRef[];

    // BC information
    static const char* _spatialDBFilename;
    static const int _id;
    static const char* _label;

    // Mesh information
    static const int _spaceDim;
    static const int _cellDim;
    static const int _numVertices;
    static const int _numCells;
    static const int _numCorners;
    static const int _cells[];

    // Calculated values.
    static const PylithScalar _tractionsCell[];
    static const PylithScalar _valsResidual[];

};

// End of file
