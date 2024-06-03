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
        class DirichletData;
    } // pylith
} // bc

class pylith::bc::DirichletData {
    // PUBLIC METHODS ///////////////////////////////////////////////////////
public:

    /// Constructor
    DirichletData(void);

    /// Destructor
    ~DirichletData(void);

    // PUBLIC MEMBERS ///////////////////////////////////////////////////////
public:

    PylithScalar tRef; ///< Reference time for rate of change of values
    PylithScalar valueRate; ///< Rate of change of value at constrained points.

    int numDOF; ///< Number of degrees of freedom at each point.
    int numFixedDOF; ///< Number of fixedDOF at constrained points.
    int numConstrainedPts; ///< Number of points constrained.

    int id; ///< Boundary condition identifier
    char* label; ///< Label for boundary condition group

    int* fixedDOF; ///< Degrees of freedom that are constrained at each point
    int* constrainedPoints; ///< Array of indices of constrained points.
    PylithScalar* valuesInitial; ///< Values at constrained points.

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
