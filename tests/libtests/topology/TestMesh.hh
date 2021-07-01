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
// Copyright (c) 2010-2021 University of California, Davis
//
// See LICENSE.md for license information.
//
// ----------------------------------------------------------------------
//

/**
 * @file tests/libtests/topology/TestMesh.hh
 *
 * @brief C++ unit testing for Mesh.
 */

#if !defined(pylith_topology_testmesh_hh)
#define pylith_topology_testmesh_hh

// Include directives ---------------------------------------------------
#include <cppunit/extensions/HelperMacros.h>

// Forward declarations -------------------------------------------------
/// Namespace for pylith package
namespace pylith {
    namespace topology {
        class TestMesh;

        class Mesh;
    } // topology
} // pylith

// TestMesh -------------------------------------------------------------
/// C++ unit testing for Mesh.
class pylith::topology::TestMesh : public CppUnit::TestFixture { // class TestMesh
    // CPPUNIT TEST SUITE /////////////////////////////////////////////////
    CPPUNIT_TEST_SUITE(TestMesh);

    CPPUNIT_TEST(testConstructor);
    CPPUNIT_TEST(testDMMesh);
    CPPUNIT_TEST(testCoordsys);
    CPPUNIT_TEST(testDimension);
    CPPUNIT_TEST(testAccessors);
    CPPUNIT_TEST(testComm);
    CPPUNIT_TEST(testView);

    CPPUNIT_TEST_SUITE_END();

    // PUBLIC METHODS /////////////////////////////////////////////////////
public:

    /// Test constructor.
    void testConstructor(void);

    /// Test dmMesh().
    void testDMMesh(void);

    /// Test coordsys().
    void testCoordsys(void);

    /// Test dimension().
    void testDimension(void);

    /// Test numCorners(), numCells(), numVertices(), isSimplex().
    void testAccessors(void);

    /// Test comm().
    void testComm(void);

    /// Test view().
    void testView(void);

}; // class TestMesh

#endif // pylith_topology_testmesh_hh

// End of file
