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
 * @file tests/libtests/bc/TestBoundaryMesh.hh
 *
 * @brief C++ TestBoundaryMesh object.
 *
 * C++ unit testing for BoundaryMesh.
 */

#if !defined(pylith_bc_testboundarymesh_hh)
#define pylith_bc_testboundarymesh_hh

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


#endif // pylith_bc_boundarymesh_hh


// End of file
