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
 * @file tests/libtests/meshio/TestDataWriterVTKSubmesh.hh
 *
 * @brief C++ TestDataWriterVTKSubmesh object
 *
 * C++ unit testing for DataWriterVTKSubmesh.
 */

#if !defined(pylith_meshio_testdatawritervtksubmesh_hh)
#define pylith_meshio_testdatawritervtksubmesh_hh

#include "TestDataWriterVTK.hh"
#include "TestDataWriterSubmesh.hh"

#include "pylith/topology/topologyfwd.hh" // USES Mesh, Field

#include <cppunit/extensions/HelperMacros.h>

/// Namespace for pylith package
namespace pylith {
    namespace meshio {
        class TestDataWriterVTKSubmesh;

        class TestDataWriterVTKSubmesh_Data;
    } // meshio
} // pylith

/// C++ unit testing for DataWriterVTK
class pylith::meshio::TestDataWriterVTKSubmesh :
    public TestDataWriterVTK,
    public TestDataWriterSubmesh,
    public CppUnit::TestFixture {
    // CPPUNIT TEST SUITE /////////////////////////////////////////////////
    CPPUNIT_TEST_SUITE(TestDataWriterVTKSubmesh);

    CPPUNIT_TEST(testTimeStep);
    CPPUNIT_TEST(testWriteVertexField);
    CPPUNIT_TEST(testWriteCellField);

    CPPUNIT_TEST_SUITE_END();

    // PUBLIC METHODS /////////////////////////////////////////////////////
public:

    /// Setup testing data.
    void setUp(void);

    /// Tear down testing data.
    void tearDown(void);

    /// Test openTimeStep() and closeTimeStep()
    void testTimeStep(void);

    /// Test writeVertexField.
    void testWriteVertexField(void);

    /// Test writeCellField.
    void testWriteCellField(void);

    // PROTECTED METHODS //////////////////////////////////////////////////
protected:

    /** Get test data.
     *
     * @returns Test data.
     */
    TestDataWriterSubmesh_Data* _getData(void);

    // PROTECTED MEMBDERS /////////////////////////////////////////////////
protected:

    TestDataWriterVTKSubmesh_Data* _data; ///< Data for testing.

}; // class TestDataWriterVTKSubmesh

// ======================================================================
class pylith::meshio::TestDataWriterVTKSubmesh_Data : public TestDataWriterVTK_Data, public TestDataWriterSubmesh_Data {};

#endif // pylith_meshio_testdatawritervtksubmesh_hh

// End of file
