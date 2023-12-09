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

#include "TestFault.hh" // USES TestFault_Data

namespace pylith {
    class PlanePWave;
}

class pylith::PlanePWave {
public:

    // Data factory methods
    static TestFault_Data* TriP1(void);

    static TestFault_Data* TriP2(void);

    static TestFault_Data* TriP3(void);

    static TestFault_Data* TriP4(void);

    static TestFault_Data* QuadQ1(void);

    static TestFault_Data* QuadQ2(void);

    static TestFault_Data* QuadQ3(void);

    static TestFault_Data* QuadQ4(void);

private:

    PlanePWave(void); ///< Not implemented
}; // PlanePWave

// End of file
