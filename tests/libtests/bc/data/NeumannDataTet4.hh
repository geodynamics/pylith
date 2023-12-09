// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2023, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

#if !defined(pylith_bc_neumanndatatet4_hh)
#define pylith_bc_neumanndatatet4_hh

#include "NeumannData.hh"

namespace pylith {
    namespace bc {
        class NeumannDataTet4;
    } // pylith
} // bc

class pylith::bc::NeumannDataTet4 : public NeumannData {
    // PUBLIC METHODS ///////////////////////////////////////////////////////
public:

    /// Constructor
    NeumannDataTet4(void);

    /// Destructor
    ~NeumannDataTet4(void);

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

#endif // pylith_bc_neumanndatatet4_hh

// End of file
