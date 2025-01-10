// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

/**
 * @file tests/libtests/bc/TestBoundaryMesh.hh
 *
 * @brief C++ TestBoundaryMesh object.
 *
 * C++ unit testing for BoundaryMesh.
 */
#pragma once

#include <cppunit/extensions/HelperMacros.h>

/// Namespace for pylith package
namespace pylith {
    namespace bc {
        class TestBoundaryMesh;

        class TestBoundaryMesh_Data; // test data
    } // bc
} // pylith

// ======================================================================
/// C++ unit testing for BoundaryMesh.
class pylith::bc::TestBoundaryMesh : public CppUnit::TestFixture {
    // CPPUNIT TEST SUITE /////////////////////////////////////////////////
    CPPUNIT_TEST_SUITE(TestBoundaryMesh);

    CPPUNIT_TEST(testSubmesh);
    CPPUNIT_TEST(testSubmeshFault);

    CPPUNIT_TEST_SUITE_END();

    // PUBLIC METHODS /////////////////////////////////////////////////////
public:

    /// Setup testing data.
    void setUp(void);

    /// Tear down testing data.
    void tearDown(void);

    /// Test submesh() without fault.
    void testSubmesh(void);

    /// Test submesh() with fault().
    void testSubmeshFault(void);

    // PROTECTED MEMBERS //////////////////////////////////////////////////
protected:

    TestBoundaryMesh_Data* _data; ///< Data for testing

}; // class TestBoundaryMesh

// ======================================================================
class pylith::bc::TestBoundaryMesh_Data {
    // PUBLIC METHODS ///////////////////////////////////////////////////////
public:

    /// Constructor
    TestBoundaryMesh_Data(void);

    /// Destructor
    ~TestBoundaryMesh_Data(void);

    // PUBLIC MEMBERS ///////////////////////////////////////////////////////
public:

    const char* filename; ///< Name of file with input mesh.
    const char* bcLabel; ///< Name of group of vertices for bc.
    const char* faultLabel; ///< Name of group of vertices for fault.
    int faultId; ///< Material identifier for fault.

    int numCorners; ///< Number of vertices in cells of boundary mesh.
    int numCells; ///< Number of cells in boundary mesh.

    /// @name Boundary mesh without fault.
    //@{
    int numVerticesNoFault; ///< Number of vertices.
    //@}

    /// @name Boundary mesh with fault.
    //@{
    int numVerticesWithFault; ///< Number of vertices.
    //@}

}; // class TestBoundaryMesh_Data

// End of file
