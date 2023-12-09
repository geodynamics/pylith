// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2023, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

/**
 * @file tests/libtests/meshio/TestOutputSolnSubset.hh
 *
 * @brief C++ TestOutputSolnSubset object
 *
 * C++ unit testing for OutputSolnSubset.
 */

#if !defined(pylith_meshio_testoutputsolnsubset_hh)
#define pylith_meshio_testoutputsolnsubset_hh

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

#endif // pylith_meshio_testoutputsolnsubset_hh

// End of file
