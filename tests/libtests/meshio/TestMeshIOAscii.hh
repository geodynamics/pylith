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
 * @file tests/libtests/meshio/TestMeshIOAscii.hh
 *
 * @brief C++ TestMeshIOAscii object
 *
 * C++ unit testing for MeshIOAscii.
 */

#if !defined(pylith_meshio_testmeshioascii_hh)
#define pylith_meshio_testmeshioascii_hh

// Include directives ---------------------------------------------------
#include "TestMeshIO.hh"

// Forward declarations -------------------------------------------------
namespace pylith {
    namespace meshio {
        class TestMeshIOAscii;

        class TestMeshIOAscii_Data; // test data
    } // meshio
} // pylith

// ======================================================================
class pylith::meshio::TestMeshIOAscii : public TestMeshIO {

    // CPPUNIT TEST SUITE ///////////////////////////////////////////////
    CPPUNIT_TEST_SUITE(TestMeshIOAscii);

    CPPUNIT_TEST(testConstructor);
    CPPUNIT_TEST(testDebug);
    CPPUNIT_TEST(testFilename);
    CPPUNIT_TEST(testWriteRead);

    CPPUNIT_TEST_SUITE_END();

    // PUBLIC METHODS ///////////////////////////////////////////////////
public:

    /// Setup testing data.
    void setUp(void);

    /// Deallocate testing data.
    void tearDown(void);

    /// Test constructor
    void testConstructor(void);

    /// Test debug()
    void testDebug(void);

    /// Test filename()
    void testFilename(void);

    /// Test write() and read().
    void testWriteRead(void);

    /// Test read().
    void testRead(void);

    /** Get test data.
     *
     * @returns Test data.
     */
    TestMeshIO_Data* _getData(void);

    // PROTECTED METHODS ////////////////////////////////////////////////
protected:

    MeshIOAscii* _io; ///< Test subject.
    TestMeshIOAscii_Data* _data; ///< Data for tests.

}; // class TestMeshIOAscii

// ======================================================================
class pylith::meshio::TestMeshIOAscii_Data : public TestMeshIO_Data {

    // PUBLIC METHODS ///////////////////////////////////////////////////
public:

    /// Constructor
    TestMeshIOAscii_Data(void);

    /// Destructor
    ~TestMeshIOAscii_Data(void);

    // PUBLIC MEMBERS ///////////////////////////////////////////////////
public:

    const char* filename;

}; // class TestMeshIOAscii_Data


#endif // pylith_meshio_testmeshioascii_hh

// End of file
