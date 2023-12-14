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
    namespace bc {
        class PointForceData;
    } // pylith
} // bc

class pylith::bc::PointForceData {
    // PUBLIC METHODS ///////////////////////////////////////////////////////
public:

    /// Constructor
    PointForceData(void);

    /// Destructor
    ~PointForceData(void);

    // PUBLIC MEMBERS ///////////////////////////////////////////////////////
public:

    PylithScalar tRef; ///< Reference time for rate of change of forces.
    PylithScalar forceRate; ///< Rate of change of force.
    PylithScalar tResidual; ///< Time for computing residual.

    int numDOF; ///< Number of degrees of freedom at each point.
    int numForceDOF; ///< Number of forces at points.
    int numForcePts; ///< Number of points with forces.

    int id; ///< Boundary condition identifier
    char* label; ///< Label for boundary condition group

    int* forceDOF; ///< Degrees of freedom that are constrained at each point
    int* forcePoints; ///< Array of indices of points with forces.
    PylithScalar* forceInitial; ///< Forces at points.
    PylithScalar* residual; ///< Residual field.

    char* meshFilename; ///< Filename for input mesh.
    char* dbFilename; ///< Filename of simple spatial database.

    /// @name Scales information for nondimensionalization.
    //@{
    PylithScalar lengthScale; ///< Length scale.
    PylithScalar pressureScale; ///< Pressure scale.
    PylithScalar timeScale; ///< Time scale.
    PylithScalar densityScale; ///< Density scale.
    //@}

};

// End of file
