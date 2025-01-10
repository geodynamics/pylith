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

#include "CohesiveImpulsesData.hh"

namespace pylith {
    namespace faults {
        class CohesiveImpulsesDataTri3;
    } // pylith
} // faults

class pylith::faults::CohesiveImpulsesDataTri3 : public CohesiveImpulsesData {
    // PUBLIC METHODS ///////////////////////////////////////////////////////
public:

    /// Constructor
    CohesiveImpulsesDataTri3(void);

    /// Destructor
    ~CohesiveImpulsesDataTri3(void);

    // PRIVATE MEMBERS //////////////////////////////////////////////////////
private:

    static const char* _meshFilename; ///< Filename of input mesh

    static const int _spaceDim; ///< Number of dimensions in vertex coordinates
    static const int _cellDim; ///< Number of dimensions associated with cell

    static const int _numBasis; ///< Number of vertices in cell
    static const int _numQuadPts; ///< Number of quadrature points
    static const PylithScalar _quadPts[]; ///< Coordinates of quad pts in ref cell
    static const PylithScalar _quadWts[]; ///< Weights of quadrature points
    static const PylithScalar _basis[]; ///< Basis fns at quadrature points
    static const PylithScalar _basisDeriv[]; ///< Derivatives of basis fns at quad pts
    static const PylithScalar _verticesRef[]; ///< Coordinates of vertices in ref cell (dual basis)

    static const int _id; ///< Fault material identifier
    static const char* _label; ///< Label for fault
    static const char* _impulseAmpFilename; ///< Name of db for impulse amplitude
    static const int _impulseDOF[]; ///< Fault components associated with impulses.
    static const int _numComponents; ///< Number of components;
    //@}

    static const PylithScalar _fieldT[]; ///< Field over domain at time t.
    static const PylithScalar _fieldIncr[]; ///< Solution increment field over domain at time t.

    static const PylithScalar _orientation[]; ///< Expected values for fault orientation.
    static const PylithScalar _area[]; ///< Expected values for fault area.
    static const PylithScalar _amplitude[]; ///< Expected values for impulse amplitude.
    static const int _numImpulses; ///< Number of impulses.
    static const PylithScalar _residual[]; ///< Expected values from residual calculation.

    static const int _constraintEdges[]; ///< Expected points for constraint edges
    static const int _negativeVertices[]; ///< Expected points for negative-side fault vertices
    static const int _numConstraintEdges; ///< Number of constraint edges

};

// End of file
