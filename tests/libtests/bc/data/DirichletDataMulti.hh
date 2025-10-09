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
    namespace bc {
        class DirichletDataMulti;
    } // pylith
} // bc

class pylith::bc::DirichletDataMulti {
    // PUBLIC METHODS ///////////////////////////////////////////////////////
public:

    /// Constructor
    DirichletDataMulti(void);

    /// Destructor
    ~DirichletDataMulti(void);

    // PUBLIC MEMBERS ///////////////////////////////////////////////////////
public:

    int numDOF; ///< Number of degrees of freedom at each point.

    //@{ Boundary condition A
    int numFixedDOFA; ///< Number of fixedDOF at constrained points.
    int numConstrainedPtsA; ///< Number of points constrained.
    int idA; ///< Boundary condition identifier
    char* labelA; ///< Label for boundary condition group
    int* fixedDOFA; ///< Degrees of freedom that are constrained at each point
    int* constrainedPointsA; ///< Array of indices of constrained points.
    char* dbFilenameA; ///< Filename for db with initial values.
    char* dbFilenameARate; ///< Filename for db with rate of change of values.
    PylithScalar tRefA; ///< Reference time for rate of change of values.
    //@}

    //@{ Boundary condition B
    int numFixedDOFB; ///< Number of fixedDOF at constrained points.
    int numConstrainedPtsB; ///< Number of points constrained.
    int idB; ///< Boundary condition identifier
    char* labelB; ///< Label for boundary condition group
    int* fixedDOFB; ///< Degrees of freedom that are constrained at each point
    int* constrainedPointsB; ///< Array of indices of constrained points.
    char* dbFilenameB; ///< Filename for db with initial values.
    char* dbFilenameBRate; ///< Filename for db with rate of change of values.
    PylithScalar tRefB; ///< Reference time for rate of change of values.
    //@}

    //@{ Boundary condition C
    int numFixedDOFC; ///< Number of fixedDOF at constrained points.
    int numConstrainedPtsC; ///< Number of points constrained.
    int idC; ///< Boundary condition identifier
    char* labelC; ///< Label for boundary condition group
    int* fixedDOFC; ///< Degrees of freedom that are constrained at each point
    int* constrainedPointsC; ///< Array of indices of constrained points.
    char* dbFilenameC; ///< Filename for db with initial values.
    char* dbFilenameCRate; ///< Filename for db with rate of change of values.
    PylithScalar tRefC; ///< Reference time for rate of change of values.
    //@}

    PylithScalar* field; ///< Values in field
    PylithScalar* fieldIncr; ///< Increment values in field
    int* constraintSizes; ///< Number of constrained DOF at each vertex
    int* constrainedDOF; ///< Indices of constrained DOF at each constrained vertex

    char* meshFilename; ///< Filename for input mesh.

    /// @name Scales information for nondimensionalization.
    //@{
    PylithScalar lengthScale; ///< Length scale.
    PylithScalar rigidityScale; ///< Pressure scale.
    PylithScalar timeScale; ///< Time scale.
    PylithScalar densityScale; ///< Density scale.
    //@}
};

// End of file
