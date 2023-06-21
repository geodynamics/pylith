// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University at Buffalo
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2022 University of California, Davis
//
// See LICENSE.md for license information.
//
// ----------------------------------------------------------------------
//

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
pylith::bc::TestPlaceholder::emptyTest(void) {}


// End of file
