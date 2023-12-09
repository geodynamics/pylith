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

#include "TestLinearElasticity.hh" // USES TestLinearElasticity_Data

namespace pylith {
    class Gravity2D;
}

class pylith::Gravity2D {
public:

    // Data factory methods

    static TestLinearElasticity_Data* TriP2(void);

    static TestLinearElasticity_Data* TriP3(void);

    static TestLinearElasticity_Data* QuadQ2(void);

    static TestLinearElasticity_Data* QuadQ3(void);

private:

    Gravity2D(void); ///< Not implemented
}; // Gravity2D

// End of file
