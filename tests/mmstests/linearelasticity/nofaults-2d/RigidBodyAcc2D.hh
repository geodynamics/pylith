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
    class RigidBodyAcc2D;
}

class pylith::RigidBodyAcc2D {
public:

    // Data factory methods
    static TestLinearElasticity_Data* TriP1(void);

    static TestLinearElasticity_Data* TriP2(void);

    static TestLinearElasticity_Data* TriP3(void);

    static TestLinearElasticity_Data* TriP4(void);

    static TestLinearElasticity_Data* QuadQ1(void);

    static TestLinearElasticity_Data* QuadQ1Distorted(void);

    static TestLinearElasticity_Data* QuadQ2(void);

    static TestLinearElasticity_Data* QuadQ3(void);

    static TestLinearElasticity_Data* QuadQ4(void);

private:

    RigidBodyAcc2D(void); ///< Not implemented
}; // RigidBodyAcc2D

// End of file
