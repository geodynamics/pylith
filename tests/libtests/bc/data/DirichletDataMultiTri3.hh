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

#include "DirichletDataMulti.hh"

namespace pylith {
    namespace bc {
        class DirichletDataMultiTri3;
    } // pylith
} // bc

class pylith::bc::DirichletDataMultiTri3 : public DirichletDataMulti {
    // PUBLIC METHODS ///////////////////////////////////////////////////////
public:

    /// Constructor
    DirichletDataMultiTri3(void);

    /// Destructor
    ~DirichletDataMultiTri3(void);

    // PRIVATE MEMBERS //////////////////////////////////////////////////////
private:

    static const int _numDOF; ///< Number of degrees of freedom at each point.

    static const int _numFixedDOFA; ///< Number of fixedDOF at constrained points.
    static const int _numConstrainedPtsA; ///< Number of points constrained.
    static const int _idA; ///< Boundary condition identifier
    static const char* _labelA; /// Label for boundary condition group
    static const int _fixedDOFA[]; ///< Degrees of freedom constrained at points
    static const int _constrainedPointsA[]; ///< Array of indices of constrained pts.
    static const char* _dbFilenameA; ///< Filename for db of initial values.
    static const char* _dbFilenameARate; ///< Filename for db of rate of change.
    static const PylithScalar _tRefA; ///< Reference time for rate of change.

    static const int _numFixedDOFB; ///< Number of fixedDOF at constrained points.
    static const int _numConstrainedPtsB; ///< Number of points constrained.
    static const int _idB; ///< Boundary condition identifier
    static const char* _labelB; /// Label for boundary condition group
    static const int _fixedDOFB[]; ///< Degrees of freedom constrained at points
    static const int _constrainedPointsB[]; ///< Array of indices of constrained pts.
    static const char* _dbFilenameB; ///< Filename for db of initial values.
    static const char* _dbFilenameBRate; ///< Filename for db of rate of change.
    static const PylithScalar _tRefB; ///< Reference time for rate of change.

    static const PylithScalar _field[]; ///< Values in field
    static const PylithScalar _fieldIncr[]; ///< Increment values in field
    static const int _constraintSizes[]; ///< Number of constrained DOF at each vertex
    static const int _constrainedDOF[]; ///< Indices of constrained DOF at each constrained vertex

    static const char* _meshFilename; ///< Filename of input mesh.
};

// End of file
