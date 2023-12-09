// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2023, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

#if !defined(pylith_bc_dirichletdataquad4_hh)
#define pylith_bc_dirichletdataquad4_hh

#include "DirichletData.hh"

namespace pylith {
    namespace bc {
        class DirichletDataQuad4;
    } // pylith
} // bc

class pylith::bc::DirichletDataQuad4 : public DirichletData {
    // PUBLIC METHODS ///////////////////////////////////////////////////////
public:

    /// Constructor
    DirichletDataQuad4(void);

    /// Destructor
    ~DirichletDataQuad4(void);

    // PRIVATE MEMBERS //////////////////////////////////////////////////////
private:

    static const int _numDOF; ///< Number of degrees of freedom at each point.

    static const int _numFixedDOF; ///< Number of fixedDOF at constrained points.
    static const int _numConstrainedPts; ///< Number of points constrained.

    static const int _id; ///< Boundary condition identifier
    static const char* _label; /// Label for boundary condition group

    static const int _fixedDOF[]; ///< Degrees of freedom constrained at points

    static const int _constrainedPoints[]; ///< Array of indices of constrained pts.
    static const PylithScalar _tRef; ///< Reference time for rate of change of value
    static const PylithScalar _valueRate; ///< Rate of change of values.
    static const PylithScalar _valuesInitial[]; ///< Initial values.

    static const char* _meshFilename; ///< Filename of input mesh.
    static const char* _dbFilename; ///< Filename of simple spatial database.
};

#endif // pylith_bc_dirichletdataquad4_hh

// End of file
