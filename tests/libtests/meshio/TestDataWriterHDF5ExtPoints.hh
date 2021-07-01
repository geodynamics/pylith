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
 * @file tests/libtests/meshio/TestDataWriterHDF5ExtPoints.hh
 *
 * @brief C++ TestDataWriterHDF5ExtPoints object
 *
 * C++ unit testing for DataWriterHDF5ExtPoints.
 */

#if !defined(pylith_meshio_testdatawriterhdf5extpoints_hh)
#define pylith_meshio_testdatawriterhdf5extpoints_hh

#include "TestDataWriterHDF5.hh" // ISA TestDataWriterHDF5
#include "TestDataWriterPoints.hh" // ISA TestDataWriterPoints

#include "pylith/topology/topologyfwd.hh" // USES Mesh, Field

#include <cppunit/extensions/HelperMacros.h>

/// Namespace for pylith package
namespace pylith {
    namespace meshio {
        class TestDataWriterHDF5ExtPoints;

        class TestDataWriterHDF5ExtPoints_Data;
    } // meshio
} // pylith

/// C++ unit testing for DataWriterHDF5Ext
class pylith::meshio::TestDataWriterHDF5ExtPoints :
    public TestDataWriterHDF5,
    public TestDataWriterPoints,
    public CppUnit::TestFixture {
    // CPPUNIT TEST SUITE /////////////////////////////////////////////////
    CPPUNIT_TEST_SUITE(TestDataWriterHDF5ExtPoints);

    CPPUNIT_TEST(testTimeStep);
    CPPUNIT_TEST(testWriteVertexField);

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

    // PROTECTED METHODS //////////////////////////////////////////////////
protected:

    /** Get test data.
     *
     * @returns Test data.
     */
    TestDataWriterPoints_Data* _getData(void);

    // PROTECTED MEMBDERS /////////////////////////////////////////////////
protected:

    TestDataWriterHDF5ExtPoints_Data* _data; ///< Data for testing.

}; // class TestDataWriterHDF5ExtPoints

// ======================================================================
class pylith::meshio::TestDataWriterHDF5ExtPoints_Data :
    public TestDataWriterHDF5_Data,
    public TestDataWriterPoints_Data {};

#endif // pylith_meshio_testdatawriterhdf5extpoints_hh

// End of file
