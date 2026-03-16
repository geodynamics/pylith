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

#include "TestMomentTensorSource.hh"

namespace pylith {
    class MomentTensorSource2D;
}

class pylith::MomentTensorSource2D {
public:

    // Data factory methods

    static TestMomentTensorSource_Data* TriP1(void);

    static TestMomentTensorSource_Data* TriP2(void);

    static TestMomentTensorSource_Data* QuadQ1(void);

    static TestMomentTensorSource_Data* QuadQ2(void);

}; // class MomentTensorSource2D

// End of file
