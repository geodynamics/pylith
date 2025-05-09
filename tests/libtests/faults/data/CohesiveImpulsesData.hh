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
        class CohesiveImpulsesData;
    } // pylith
} // faults

class pylith::faults::CohesiveImpulsesData {
    // PUBLIC METHODS ///////////////////////////////////////////////////////
public:

    /// Constructor
    CohesiveImpulsesData(void);

    /// Destructor
    ~CohesiveImpulsesData(void);

    // PUBLIC MEMBERS ///////////////////////////////////////////////////////
public:

    char* meshFilename; ///< Filename for input mesh

    /// @name Scales information for nondimensionalization.
    //@{
    PylithScalar lengthScale; ///< Length scale.
    PylithScalar pressureScale; ///< Pressure scale.
    PylithScalar timeScale; ///< Time scale.
    PylithScalar densityScale; ///< Density scale.
    //@}

    /// @name Quadrature information
    //@{
    int spaceDim; ///< Number of dimensions in vertex coordinates
    int cellDim; ///< Number of dimensions associated with cell
    int numBasis; ///< Number of vertices in cell
    int numQuadPts; ///< Number of quadrature points
    PylithScalar* quadPts; ///< Coordinates of quad pts in ref cell
    PylithScalar* quadWts; ///< Weights of quadrature points
    PylithScalar* basis; ///< Basis fns at quadrature points
    PylithScalar* basisDeriv; ///< Derivatives of basis fns at quad pts
    PylithScalar* verticesRef; ///< Coordinates of vertices in ref cell (dual basis)
    //@}

    /// @name Fault information
    //@{
    int id; ///< Fault material identifier
    char* label; ///< Label for fault
    char* impulseAmpFilename; ///< Name of db for final slip
    int* impulseDOF; ///< Fault DOF associated with impulses.
    int numComponents; ///< Number of components in impulses.
    //@}

    /// @name Input fields
    //@{
    PylithScalar* fieldT; ///< Field over domain at time t.
    PylithScalar* fieldIncr; ///< Solution increment field over domain at time t.
    //@}

    /// @name Calculated values.
    //@{
    PylithScalar* orientation; ///< Expected values for fault orientation.
    PylithScalar* area; ///< Expected values for fault area.
    PylithScalar* amplitude; ///< Expected values for impulse amplitude.
    int numImpulses; ///< Number of impulses.

    /// Expected values from residual calculation.
    PylithScalar* residual;

    int* constraintEdges; ///< Expected points for constraint edges
    int* negativeVertices; ///< Expected points for negative side fault vertices
    int numConstraintEdges; ///< Number of constraint edges
    //@}

};

// End of file
