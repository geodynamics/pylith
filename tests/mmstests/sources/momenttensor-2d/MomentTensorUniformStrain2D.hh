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
    class MomentTensorUniformStrain2D;
}

/** MMS test with nonzero uniform strain solution for moment tensor source.
 *
 * The manufactured solution uses a linear displacement field (uniform strain)
 * with the moment tensor source contribution.
 */
class pylith::MomentTensorUniformStrain2D {
public:

    // Data factory methods

    static TestMomentTensorSource_Data* TriP1(void);

    static TestMomentTensorSource_Data* TriP2(void);

    static TestMomentTensorSource_Data* QuadQ1(void);

    static TestMomentTensorSource_Data* QuadQ2(void);

}; // class MomentTensorUniformStrain2D

// End of file
