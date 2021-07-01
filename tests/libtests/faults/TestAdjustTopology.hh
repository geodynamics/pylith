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
 * @file tests/libtests/faults/TestAdjustTopologyt.hh
 *
 * C++ unit tests for FaultCohesive::adjustTopology().
 */

#if !defined(pylith_faults_testadjusttopology_hh)
#define pylith_faults_testadjusttopology_hh

#include <cppunit/extensions/HelperMacros.h>

#include "pylith/topology/topologyfwd.hh" // USES PETSc Mesh

/// Namespace for pylith package
namespace pylith {
    namespace faults {
        class TestAdjustTopology;
        class TestAdjustTopology_Data;
    } // faults
} // pylith

/// C++ unit testing for Fault
class pylith::faults::TestAdjustTopology : public CppUnit::TestFixture {
    // CPPUNIT TEST SUITE /////////////////////////////////////////////////////////////////////////
    CPPUNIT_TEST_SUITE(TestAdjustTopology);

    CPPUNIT_TEST(testAdjustTopology);

    CPPUNIT_TEST_SUITE_END();

    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    /// Setup testing data.
    void setUp(void);

    /// Deallocate testing data.
    void tearDown(void);

    /// Test adjustTopology().
    void testAdjustTopology(void);

    // PROTECTED MEMBERS //////////////////////////////////////////////////////////////////////////
protected:

    TestAdjustTopology_Data* _data; ///< Data for testing.
    pylith::topology::Mesh* _mesh; ///< Finite-element mesh.

    // PRIVATE METHODS ////////////////////////////////////////////////////////////////////////////
private:

    /// Setup mesh.
    void _initialize();

}; // class TestAdjustTopology

class pylith::faults::TestAdjustTopology_Data {
    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    /// Constructor
    TestAdjustTopology_Data(void);

    /// Destructor
    ~TestAdjustTopology_Data(void);

    // PUBLIC MEMBERS /////////////////////////////////////////////////////////////////////////////
public:

    const char* filename; ///< Name of mesh file.

    size_t numFaults; ///< Number of faults
    const char** faultSurfaceLabels; ///< Labels marking fault surfaces.
    const char** faultEdgeLabels; ///< Labels for buried edges.
    const int* interfaceIds; ///< Label values for interfaces.

    size_t cellDim; ///< Number of dimensions associated with cells.
    size_t spaceDim; ///< Spatial dimension for vertex coordinates.
    size_t numVertices; // Number of vertices in updated mesh.
    size_t numCells; ///< Number of cells in updated mesh.

    const int* numCorners; ///< Number of vertices in each cell in updated mesh.
    const int* materialIds; ///< Material identifief for each cell in updated mesh.

    size_t numGroups; ///< Number of groups.
    const int* groupSizes; ///< Array of sizes in each group.
    char** groupNames; ///< Names of groups.
    char** groupTypes; ///< Type of group.

    bool failureExpected; ///< Flag indicating adjust topology should fail.

}; // TestAdjustTopology_Data

#endif // pylith_faults_testadjusttopology_hh

// End of file
