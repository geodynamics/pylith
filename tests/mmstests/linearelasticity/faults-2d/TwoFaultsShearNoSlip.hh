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

#include "TestFaultKin.hh" // USES TestFaultKin_Data

namespace pylith {
    class TwoFaultsShearNoSlip;
}

class pylith::TwoFaultsShearNoSlip {
public:

    // Data factory methods
    static TestFaultKin_Data* TriP1(void);

    static TestFaultKin_Data* TriP2(void);

    static TestFaultKin_Data* TriP3(void);

    static TestFaultKin_Data* TriP4(void);

    static TestFaultKin_Data* QuadQ1(void);

    static TestFaultKin_Data* QuadQ2(void);

    static TestFaultKin_Data* QuadQ3(void);

    static TestFaultKin_Data* QuadQ4(void);

private:

    TwoFaultsShearNoSlip(void); ///< Not implemented
}; // TwoFaultsShearNoSlip

// End of file
