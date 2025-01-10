// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

#include <portinfo>

#include "pylith/utils/GenericComponent.hh" // ISA GenericComponent

#include "catch2/catch_test_macros.hpp"

// ------------------------------------------------------------------------------------------------
namespace pylith {
    namespace bc {
        class TestPlaceholder;
    }
}

class pylith::bc::TestPlaceholder : public pylith::utils::GenericComponent {
public:

    /// Test methodCount().
    static
    void emptyTest(void);

}; // class TestPlaceholder

// ------------------------------------------------------------------------------------------------
TEST_CASE("TestPlaceholder::emptyTest", "[TestPlaceholder]") {
    pylith::bc::TestPlaceholder::emptyTest();
}

// ------------------------------------------------------------------------------------------------
void
pylith::bc::TestPlaceholder::emptyTest(void) {
}


// End of file
