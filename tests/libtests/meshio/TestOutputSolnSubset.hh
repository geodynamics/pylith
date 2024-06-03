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
        class TestOutputSolnSubset;
    } // meshio
} // pylith

/// C++ unit testing for OutputSolnSubset
class pylith::meshio::TestOutputSolnSubset : public CppUnit::TestFixture { // class TestOutputSolnSubset
                                                                           // CPPUNIT TEST SUITE
                                                                           // /////////////////////////////////////////////////
    CPPUNIT_TEST_SUITE(TestOutputSolnSubset);

    CPPUNIT_TEST(testConstructor);
    CPPUNIT_TEST(testLabel);
    CPPUNIT_TEST(testSubdomainMesh);

    CPPUNIT_TEST_SUITE_END();

    // PUBLIC METHODS /////////////////////////////////////////////////////
public:

    /// Test constructor
    void testConstructor(void);

    /// Test getLabel()
    void testLabel(void);

    /// Test subdomainMesh()
    void testSubdomainMesh(void);

}; // class TestOutputSolnSubset

// End of file
