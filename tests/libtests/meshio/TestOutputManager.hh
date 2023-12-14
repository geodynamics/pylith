// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2023, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================
#pragma once

#include <cppunit/extensions/HelperMacros.h>

/// Namespace for pylith package
namespace pylith {
    namespace meshio {
        class TestOutputManager;
    } // meshio
} // pylith

/// C++ unit testing for OutputManager
class pylith::meshio::TestOutputManager : public CppUnit::TestFixture { // class TestOutputManager
                                                                        // CPPUNIT TEST SUITE
                                                                        // /////////////////////////////////////////////////
    CPPUNIT_TEST_SUITE(TestOutputManager);

    CPPUNIT_TEST(testConstructor);
    CPPUNIT_TEST(testCoordsys);
    CPPUNIT_TEST(testWriter);
    CPPUNIT_TEST(testVertexFilter);
    CPPUNIT_TEST(testCellFilter);
    CPPUNIT_TEST(testOpenClose);
    CPPUNIT_TEST(testOpenCloseTimeStep);
    CPPUNIT_TEST(testAppendVertexField);
    CPPUNIT_TEST(testAppendCellField);

    CPPUNIT_TEST_SUITE_END();

    // PUBLIC METHODS /////////////////////////////////////////////////////
public:

    /// Test constructor
    void testConstructor(void);

    /// Test coordsys()
    void testCoordsys(void);

    /// Test writer()
    void testWriter(void);

    /// Test vertexFilter()
    void testVertexFilter(void);

    /// Test cellFilter().
    void testCellFilter(void);

    /// Test open() and close().
    void testOpenClose(void);

    /// Test openTimeStep() and closeTimeStep().
    void testOpenCloseTimeStep(void);

    /// Test appendVertexField().
    void testAppendVertexField(void);

    /// Test appendCellField().
    void testAppendCellField(void);

}; // class TestOutputManager

// End of file
