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
 * @file tests/libtests/topology/TestSubmesh.hh
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
        class TestSubmesh;
        class TestSubmesh_Data;
    } // topology
} // pylith

// TestSubmesh -----------------------------------------------------------------
/// C++ unit testing for Mesh.
class pylith::topology::TestSubmesh : public CppUnit::TestFixture {
    // CPPUNIT TEST SUITE //////////////////////////////////////////////////////
    CPPUNIT_TEST_SUITE(TestSubmesh);

    CPPUNIT_TEST(testCreateLowerDimMesh);
    CPPUNIT_TEST(testCreateSubdomainMesh);
    CPPUNIT_TEST(testAccessors);
    CPPUNIT_TEST(testSizes);

    CPPUNIT_TEST_SUITE_END();

    // PUBLIC METHODS //////////////////////////////////////////////////////////
public:

    /// Setup testing data.
    void setUp(void);

    /// Deallocate testing data.
    void tearDown(void);

    /// Test MeshOps::createLowerDimMesh().
    void testCreateLowerDimMesh(void);

    /// Test MeshOps::testCreateSubdomainMesh().
    void testCreateSubdomainMesh(void);

    /// Test coordsys(), debug(), comm().
    void testAccessors(void);

    /// Test dimension(), numCorners(), numVertices(), numCells(),
    void testSizes(void);

    // PROTECTED METHODS ///////////////////////////////////////////////////////
protected:

    // Build lower dimension mesh.
    void _buildMesh(void);

    // PROTECTED MEMBERS ///////////////////////////////////////////////////////
protected:

    TestSubmesh_Data* _data; ///< Data for testing.
    Mesh* _domainMesh; ///< Mesh holding domain mesh.
    Mesh* _testMesh; ///< Test subject.

}; // class TestSubmesh

// TestSubmesh_Data-------------------------------------------------------------
class pylith::topology::TestSubmesh_Data {
    // PUBLIC METHODS //////////////////////////////////////////////////////////
public:

    /// Constructor
    TestSubmesh_Data(void);

    /// Destructor
    ~TestSubmesh_Data(void);

    // PUBLIC MEMBERS //////////////////////////////////////////////////////////
public:

    // GENERAL, VALUES DEPEND ON TEST CASE

    /// @defgroup Domain mesh information.
    /// @{
    int cellDim; ///< Cell dimension (matches space dimension).
    int numVertices; ///< Number of vertices.
    int numCells; ///< Number of cells.
    int numCorners; ///< Number of vertices per cell.
    int* cells; ///< Array of vertices in cells [numCells*numCorners].
    PylithScalar* coordinates; ///< Coordinates of vertices [numVertices*cellDim].
    /// @}

    /// @defgroup Submesh information.
    /// @{
    const char* groupLabel; ///< Label of group associated with submesh.
    int groupSize; ///< Number of vertices in submesh group.
    int* groupVertices; ///< Array of vertices in submesh group.
    int submeshNumCorners; ///< Number of vertices per cell.
    int submeshNumVertices; ///< Number of vertices in submesh.
    int* submeshVertices; ///< Vertices in submesh.
    int submeshNumCells; ///< Number of cells in submesh.
    int* submeshCells; ///< Array of vertices in cells [submeshNumCells*submeshNumCorners].
    /// @}

    /// @defgroup Subdomain information.
    /// @{
    const char* subdomainLabel; ///< Label defining subdomains.
    int* subdomainLabelValues; ///< Label values defining subdomains.
    int subdomainLabelValue; ///< Label value for target subdomain.
    int subdomainNumCorners; ///< Number of vertices per cell.
    int subdomainNumVertices; ///< Number of vertices in subdomain.
    int* subdomainVertices; ///< Vertices in subdomain.
    int subdomainNumCells; ///< Number of cells in subdomain.
    int* subdomainCells; ///< Array of vertices in cells [subdomainNumCells*subdomainNumCorners].
    /// @}

}; // TestSubmesh_Data

#endif // pylith_topology_testsubmesh_hh

// End of file
