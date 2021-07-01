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
 * @file tests/libtests/meshio/TestDataWriterVTKPoints.hh
 *
 * @brief C++ TestDataWriterVTKPoints object
 *
 * C++ unit testing for DataWriterVTKPoints.
 */

#if !defined(pylith_meshio_testdatawritervtkpoints_hh)
#define pylith_meshio_testdatawritervtkpoints_hh

#include "TestDataWriterVTK.hh" // ISA TestDataWriterVTK
#include "TestDataWriterPoints.hh" // ISA TestDataWriterPoints

#include "pylith/topology/topologyfwd.hh" // USES Mesh, Field

#include <cppunit/extensions/HelperMacros.h>

/// Namespace for pylith package
namespace pylith {
    namespace meshio {
        class TestDataWriterVTKPoints;

        class TestDataWriterVTKPoints_Data;
    } // meshio
} // pylith

class pylith::meshio::TestDataWriterVTKPoints :
    public TestDataWriterVTK,
    public TestDataWriterPoints,
    public CppUnit::TestFixture {
    // CPPUNIT TEST SUITE /////////////////////////////////////////////////
    CPPUNIT_TEST_SUITE(TestDataWriterVTKPoints);

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

    TestDataWriterVTKPoints_Data* _data; ///< Data for testing.

}; // class TestDataWriterVTKPoints

// ======================================================================
class pylith::meshio::TestDataWriterVTKPoints_Data : public TestDataWriterVTK_Data, public TestDataWriterPoints_Data {};

#endif // pylith_meshio_testdatawritervtkpoints_hh

// End of file
