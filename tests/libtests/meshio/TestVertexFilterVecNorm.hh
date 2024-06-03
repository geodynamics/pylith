// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2024, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================
#pragma once

#include <cppunit/extensions/HelperMacros.h>

/// Namespace for pylith package
namespace pylith {
    namespace meshio {
        class TestVertexFilterVecNorm;
    } // meshio
} // pylith

/// C++ unit testing for VertexFilterVecNorm
class pylith::meshio::TestVertexFilterVecNorm : public CppUnit::TestFixture { // class TestVertexFilterVecNorm
                                                                              // CPPUNIT TEST SUITE
                                                                              // /////////////////////////////////////////////////
    CPPUNIT_TEST_SUITE(TestVertexFilterVecNorm);

    CPPUNIT_TEST(testConstructor);
    CPPUNIT_TEST(testFilter);

    CPPUNIT_TEST_SUITE_END();

    // PUBLIC METHODS /////////////////////////////////////////////////////
public:

    /// Test constructor
    void testConstructor(void);

    /// Test filter()
    void testFilter(void);

}; // class TestVertexFilterVecNorm

// End of file
