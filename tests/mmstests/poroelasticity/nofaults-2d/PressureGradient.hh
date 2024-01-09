// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2023, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

#include <portinfo>

#include "TestLinearPoroelasticity.hh" // USES TestLinearPoroelasticity_Data

namespace pylith {
    class PressureGradient;
}

class pylith::PressureGradient {
public:

    // Data factory methods

    static TestLinearPoroelasticity_Data* TriP2P1P1(void);

    static TestLinearPoroelasticity_Data* TriP3P2P2(void);

    static TestLinearPoroelasticity_Data* QuadQ2Q1Q1(void);

    static TestLinearPoroelasticity_Data* QuadQ3Q2Q2(void);

    static TestLinearPoroelasticity_Data* TriP2P1P1_StateVars(void);

    static TestLinearPoroelasticity_Data* TriP3P2P2_StateVars(void);

    static TestLinearPoroelasticity_Data* QuadQ2Q1Q1_StateVars(void);

    static TestLinearPoroelasticity_Data* QuadQ3Q2Q2_StateVars(void);

private:

    PressureGradient(void); ///< Not implemented
}; // PressureGradient

// End of file
