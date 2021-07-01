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
 * @file tests/libtests/meshio/TestDataWriterVTKMesh.hh
 *
 * @brief C++ TestDataWriterVTKMesh object
 *
 * C++ unit testing for DataWriterVTKMesh.
 */

#if !defined(pylith_meshio_testdatawritervtkmesh_hh)
#define pylith_meshio_testdatawritervtkmesh_hh

#include "TestDataWriterVTK.hh" // ISA TestDataWriterVTK
#include "TestDataWriterMesh.hh" // ISA TestDataWriterMesh

#include "pylith/topology/topologyfwd.hh" // USES Mesh, Field

#include <cppunit/extensions/HelperMacros.h>

/// Namespace for pylith package
namespace pylith {
    namespace meshio {
        class TestDataWriterVTKMesh;

        class TestDataWriterVTKMesh_Data;
    } // meshio
} // pylith

// ======================================================================
/// C++ unit testing for DataWriterVTK
class pylith::meshio::TestDataWriterVTKMesh : public TestDataWriterVTK, public TestDataWriterMesh, public CppUnit::TestFixture {

    // CPPUNIT TEST SUITE /////////////////////////////////////////////////
    CPPUNIT_TEST_SUITE(TestDataWriterVTKMesh);

    CPPUNIT_TEST(testConstructor);
    CPPUNIT_TEST(testFilename);
    CPPUNIT_TEST(testTimeFormat);
    CPPUNIT_TEST(testTimeConstant);
    CPPUNIT_TEST(testPrecision);
    CPPUNIT_TEST(testTimeStep);
    CPPUNIT_TEST(testWriteVertexField);
    CPPUNIT_TEST(testWriteCellField);
    CPPUNIT_TEST(testVtkFilename);

    CPPUNIT_TEST_SUITE_END();

    // PUBLIC METHODS /////////////////////////////////////////////////////
public:

    /// Setup testing data.
    void setUp(void);

    /// Tear down testing data.
    void tearDown(void);

    /// Test constructor
    void testConstructor(void);

    /// Test filename()
    void testFilename(void);

    /// Test timeFormat()
    void testTimeFormat(void);

    /// Test timeConstant()
    void testTimeConstant(void);

    /// Test precision()
    void testPrecision(void);

    /// Test openTimeStep() and closeTimeStep()
    void testTimeStep(void);

    /// Test writeVertexField.
    void testWriteVertexField(void);

    /// Test writeCellField.
    void testWriteCellField(void);

    /// Test vtkFilename.
    void testVtkFilename(void);

    // PROTECTED METHODS //////////////////////////////////////////////////
protected:

    /** Get test data.
     *
     * @returns Test data.
     */
    TestDataWriter_Data* _getData(void);


    // PROTECTED MEMBDERS /////////////////////////////////////////////////
protected:

    TestDataWriterVTKMesh_Data* _data; ///< Data for testing.

}; // class TestDataWriterVTKMesh


// ======================================================================
class pylith::meshio::TestDataWriterVTKMesh_Data : public TestDataWriterVTK_Data, public TestDataWriter_Data {};

#endif // pylith_meshio_testdatawritervtkmesh_hh


// End of file
