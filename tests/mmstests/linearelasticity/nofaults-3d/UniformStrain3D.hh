// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2024, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

#include <portinfo>

#include "TestLinearElasticity.hh" // USES TestLinearElasticity_Data

namespace pylith {
    class UniformStrain3D;
}

class pylith::UniformStrain3D {
public:

    // Data factory methods
    static TestLinearElasticity_Data* TetP1(void);

    static TestLinearElasticity_Data* TetP2(void);

    static TestLinearElasticity_Data* TetP3(void);

    static TestLinearElasticity_Data* TetP4(void);

    static TestLinearElasticity_Data* HexQ1(void);

    static TestLinearElasticity_Data* HexQ2(void);

    static TestLinearElasticity_Data* HexQ3(void);

    static TestLinearElasticity_Data* HexQ4(void);

private:

    UniformStrain3D(void); ///< Not implemented
}; // UniformStrain3D

// End of file
