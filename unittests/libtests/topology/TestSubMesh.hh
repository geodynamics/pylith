// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2017 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

/**
 * @file unittests/libtests/topology/TestSubMesh.hh
 *
 * @brief C++ unit testing for Mesh.
 */

#if !defined(pylith_topology_testsubmesh_hh)
#define pylith_topology_testsubmesh_hh

// Include directives ----------------------------------------------------------
#include <cppunit/extensions/HelperMacros.h>

#include "pylith/topology/topologyfwd.hh" // forward declarations
#include "pylith/utils/types.hh" // USES PylithScalar

// Forward declarations --------------------------------------------------------
/// Namespace for pylith package
namespace pylith {
    namespace topology {
        class TestSubMesh;
        class TestSubMesh_Data;
    }   // topology
}   // pylith

// TestSubMesh -----------------------------------------------------------------
/// C++ unit testing for Mesh.
class pylith::topology::TestSubMesh : public CppUnit::TestFixture {

    // CPPUNIT TEST SUITE //////////////////////////////////////////////////////
    CPPUNIT_TEST_SUITE( TestSubMesh );

    CPPUNIT_TEST( testConstructor );
    CPPUNIT_TEST( testConstructorMesh );
    CPPUNIT_TEST( testCoordsys );
    CPPUNIT_TEST( testDebug );
    CPPUNIT_TEST( testDimension );
    CPPUNIT_TEST( testNumCorners );
    CPPUNIT_TEST( testNumVertices );
    CPPUNIT_TEST( testNumCells );
    CPPUNIT_TEST( testComm );

    CPPUNIT_TEST_SUITE_END();

    // PUBLIC METHODS //////////////////////////////////////////////////////////
public:

    /// Setup testing data.
    void setUp(void);

    /// Deallocate testing data.
    void tearDown(void);

    /// Test constructor.
    void testConstructor(void);

    /// Test constructor w/mesh.
    void testConstructorMesh(void);

    /// Test coordsys().
    void testCoordsys(void);

    /// Test debug().
    void testDebug(void);

    /// Test dimension().
    void testDimension(void);

    /// Test numCorners().
    void testNumCorners(void);

    /// Test numVertices().
    void testNumVertices(void);

    /// Test numCells().
    void testNumCells(void);

    /// Test comm().
    void testComm(void);

    // PROTECTED METHODS ///////////////////////////////////////////////////////
protected:

    /** Build mesh.
     *
     * @param mesh Finite-element mesh.
     */
    void _buildMesh(Mesh* mesh);

    // PROTECTED MEMBERS ///////////////////////////////////////////////////////
protected:

    TestSubMesh_Data* _data; ///< Data for testing.

};  // class TestSubMesh


// TestSubMesh_Data-------------------------------------------------------------
class pylith::topology::TestSubMesh_Data {

    // PUBLIC METHODS //////////////////////////////////////////////////////////
public:

    /// Constructor
    TestSubMesh_Data(void);

    /// Destructor
    ~TestSubMesh_Data(void);

    // PUBLIC MEMBERS //////////////////////////////////////////////////////////
public:

    // GENERAL, VALUES DEPEND ON TEST CASE

    /// @defgroup Domain mesh information.
    /// @{
    int cellDim; ///< Cell dimension (matches space dimension).
    int numVertices; ///< Number of vertices.
    int numCells;   ///< Number of cells.
    int numCorners; ///< Number of vertices per cell.
    int* cells; ///< Array of vertices in cells [numCells*numCorners].
    PylithScalar* coordinates;  ///< Coordinates of vertices [numVertices*cellDim].
    const char* label;  ///< Label of group associated with submesh.
    int groupSize;  ///< Number of vertices in group.
    int* groupVertices; ///< Array of vertices in group.
    /// @}

    /// @defgroup SubMesh information.
    /// @{
    int submeshNumCorners;  ///< Number of vertices per cell.
    int submeshNumVertices; ///< Number of vertices in submesh.
    int* submeshVertices;   ///< Vertices in submesh.
    int submeshNumCells; ///< Number of cells in submesh.
    int* submeshCells;  ///< Array of vertices in cells [submeshNumCells*submeshNumCorners].
    /// @}

};  // TestSubMesh_Data


#endif  // pylith_topology_testsubmesh_hh


// End of file
