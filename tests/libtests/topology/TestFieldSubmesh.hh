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

// Include directives ---------------------------------------------------
#include <cppunit/extensions/HelperMacros.h>

#include "pylith/topology/topologyfwd.hh" // forward declarations

// Forward declarations -------------------------------------------------
/// Namespace for pylith package
namespace pylith {
    namespace topology {
        class TestFieldSubmesh;
    } // topology
} // pylith

// TestFieldSubmesh -----------------------------------------------------
/// C++ unit testing for Field.
class pylith::topology::TestFieldSubmesh : public CppUnit::TestFixture { // class TestFieldSubmesh
                                                                         // CPPUNIT TEST SUITE
                                                                         // /////////////////////////////////////////////////
    CPPUNIT_TEST_SUITE(TestFieldSubmesh);

    CPPUNIT_TEST(testConstructor);
    CPPUNIT_TEST(testSection);
    CPPUNIT_TEST(testMesh);
    CPPUNIT_TEST(testSpaceDim);
    CPPUNIT_TEST(testNewSectionPoints);
    CPPUNIT_TEST(testNewSectionDomain);
    CPPUNIT_TEST(testNewSectionField);
    CPPUNIT_TEST(testCloneSection);
    CPPUNIT_TEST(testClear);
    CPPUNIT_TEST(testAllocate);
    CPPUNIT_TEST(testZero);
    CPPUNIT_TEST(testComplete);
    CPPUNIT_TEST(testCopy);
    CPPUNIT_TEST(testOperatorAdd);
    CPPUNIT_TEST(testDimensionalize);
    CPPUNIT_TEST(testView);
    CPPUNIT_TEST(testCreateScatter);
    CPPUNIT_TEST(testCreateScatterWithBC);
    CPPUNIT_TEST(testVector);
    CPPUNIT_TEST(testScatterLocalToGlobal);
    CPPUNIT_TEST(testScatterGlobalToLocal);

    CPPUNIT_TEST_SUITE_END();

    // PUBLIC METHODS /////////////////////////////////////////////////////
public:

    /// Test constructor.
    void testConstructor(void);

    /// Test section().
    void testSection(void);

    /// Test mesh().
    void testMesh(void);

    /// Test getSpaceDim().
    void testSpaceDim(void);

    /// Test newSection(points).
    void testNewSectionPoints(void);

    /// Test newSection(domain).
    void testNewSectionDomain(void);

    /// Test newSection(field).
    void testNewSectionField(void);

    /// Test cloneSection().
    void testCloneSection(void);

    /// Test clear().
    void testClear(void);

    /// Test allocate().
    void testAllocate(void);

    /// Test zero().
    void testZero(void);

    /// Test complete().
    void testComplete(void);

    /// Test copy().
    void testCopy(void);

    /// Test operator+=().
    void testOperatorAdd(void);

    /// Test dimensionalize().
    void testDimensionalize(void);

    /// Test view().
    void testView(void);

    /// Test createScatter().
    void testCreateScatter(void);

    /// Test createScatterWithBC().
    void testCreateScatterWithBC(void);

    /// Test vector().
    void testVector(void);

    /// Test scatterLocalToGlobal().
    void testScatterLocalToGlobal(void);

    /// Test scatterGlobalToLocal().
    void testScatterGlobalToLocal(void);

    // PRIVATE METHODS /////////////////////////////////////////////////////
private:

    /** Build mesh.
     *
     * @param mesh Finite-element mesh.
     */
    static
    void _buildMesh(Mesh* mesh);

}; // class TestFieldSubmesh

// End of file
