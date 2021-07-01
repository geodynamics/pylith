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
 * @file tests/libtests/meshio/TestDataWriterHDF5ExtMesh.hh
 *
 * @brief C++ TestDataWriterHDF5ExtMesh object
 *
 * C++ unit testing for DataWriterHDF5ExtMesh.
 */

#if !defined(pylith_meshio_testdatawriterhdf5extmesh_hh)
#define pylith_meshio_testdatawriterhdf5extmesh_hh

#include "TestDataWriterHDF5.hh" // ISA TestDataWriterHDF5
#include "TestDataWriterMesh.hh" // ISA TestDataWriterMesh

#include <cppunit/extensions/HelperMacros.h>

/// Namespace for pylith package
namespace pylith {
    namespace meshio {
        class TestDataWriterHDF5ExtMesh;

        class TestDataWriterHDF5ExtMesh_Data;
    } // meshio
} // pylith

/// C++ unit testing for DataWriterHDF5
class pylith::meshio::TestDataWriterHDF5ExtMesh : public TestDataWriterHDF5, public TestDataWriterMesh, public CppUnit::TestFixture {

    // CPPUNIT TEST SUITE /////////////////////////////////////////////////
    CPPUNIT_TEST_SUITE(TestDataWriterHDF5ExtMesh);

    CPPUNIT_TEST(testConstructor);
    CPPUNIT_TEST(testFilename);
    CPPUNIT_TEST(testOpenClose);
    CPPUNIT_TEST(testWriteVertexField);
    CPPUNIT_TEST(testWriteCellField);
    CPPUNIT_TEST(testHdf5Filename);
    CPPUNIT_TEST(testDatasetFilename);

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

    /// Test open() and close()
    void testOpenClose(void);

    /// Test writeVertexField.
    void testWriteVertexField(void);

    /// Test writeCellField.
    void testWriteCellField(void);

    /// Test hdf5Filename.
    void testHdf5Filename(void);

    /// Test datasetFilename.
    void testDatasetFilename(void);

    // PROTECTED METHODS //////////////////////////////////////////////////
protected:

    /** Get test data.
     *
     * @returns Test data.
     */
    TestDataWriter_Data* _getData(void);

    // PROTECTED MEMBDERS /////////////////////////////////////////////////
protected:

    TestDataWriterHDF5ExtMesh_Data* _data; ///< Data for testing.

}; // class TestDataWriterHDF5ExtMesh


// ======================================================================
class pylith::meshio::TestDataWriterHDF5ExtMesh_Data : public TestDataWriterHDF5_Data, public TestDataWriter_Data {};


#endif // pylith_meshio_testdatawriterhdf5extmesh_hh


// End of file
