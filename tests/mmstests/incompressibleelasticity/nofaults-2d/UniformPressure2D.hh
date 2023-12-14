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

#include "TestIncompressibleElasticity.hh" // USES TestIncompressibleElasticity_Data

namespace pylith {
    class UniformPressure2D;
}

class pylith::UniformPressure2D {
public:

    // Data factory methods
    static TestIncompressibleElasticity_Data* TriP1(void);

    static TestIncompressibleElasticity_Data* TriP2(void);

    static TestIncompressibleElasticity_Data* TriP3(void);

    static TestIncompressibleElasticity_Data* TriP4(void);

    static TestIncompressibleElasticity_Data* QuadQ1(void);

    static TestIncompressibleElasticity_Data* QuadQ1Distorted(void);

    static TestIncompressibleElasticity_Data* QuadQ2(void);

    static TestIncompressibleElasticity_Data* QuadQ3(void);

    static TestIncompressibleElasticity_Data* QuadQ4(void);

private:

    UniformPressure2D(void); ///< Not implemented
}; // UniformPressure2D

// End of file
