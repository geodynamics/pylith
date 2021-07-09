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

/**
 * @file tests/libtests/meshio/TestDataWriterVTKMaterial.hh
 *
 * @brief C++ TestDataWriterVTKMaterial object
 *
 * C++ unit testing for DataWriterVTKMaterial.
 */

#if !defined(pylith_meshio_testdatawritervtkmaterial_hh)
#define pylith_meshio_testdatawritervtkmaterial_hh

#include "TestDataWriterVTK.hh" // ISA TestDataWriterVTK
#include "TestDataWriterMaterial.hh" // ISA TestDataWriterMaterial

#include "pylith/topology/topologyfwd.hh" // USES Mesh, Field

#include <cppunit/extensions/HelperMacros.h>

/// Namespace for pylith package
namespace pylith {
    namespace meshio {
        class TestDataWriterVTKMaterial;

        class TestDataWriterVTKMaterial_Data;
    } // meshio
} // pylith

// ======================================================================
/// C++ unit testing for DataWriterVTK
class pylith::meshio::TestDataWriterVTKMaterial : public TestDataWriterVTK, public TestDataWriterMaterial, public CppUnit::TestFixture {
    // CPPUNIT TEST SUITE /////////////////////////////////////////////////
    CPPUNIT_TEST_SUITE(TestDataWriterVTKMaterial);

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
    TestDataWriterMaterial_Data* _getData(void);

    // PROTECTED MEMBDERS /////////////////////////////////////////////////
protected:

    TestDataWriterVTKMaterial_Data* _data; ///< Data for testing.

}; // class TestDataWriterVTKMaterial

// ======================================================================
class pylith::meshio::TestDataWriterVTKMaterial_Data : public TestDataWriterVTK_Data, public TestDataWriterMaterial_Data {};

#endif // pylith_meshio_testdatawritervtkmaterial_hh

// End of file
