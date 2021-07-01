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
 * @file tests/libtests/meshio/TestMeshIOLagrit.hh
 *
 * @brief C++ TestMeshIOLagrit object
 *
 * C++ unit testing for MeshIOLagrit.
 */

#if !defined(pylith_meshio_testmeshiolagrit_hh)
#define pylith_meshio_testmeshiolagrit_hh

// Include directives ---------------------------------------------------
#include "TestMeshIO.hh"

// Forward declarations -------------------------------------------------
namespace pylith {
    namespace meshio {
        class TestMeshIOLagrit;

        class TestMeshIOLagrit_Data; // test data
    } // meshio
} // pylith

// TestMeshIOLagrit -----------------------------------------------------
class pylith::meshio::TestMeshIOLagrit : public TestMeshIO { // class TestMeshIOLagrit

    // CPPUNIT TEST SUITE /////////////////////////////////////////////////
    CPPUNIT_TEST_SUITE(TestMeshIOLagrit);

    CPPUNIT_TEST(testConstructor);
    CPPUNIT_TEST(testDebug);
    CPPUNIT_TEST(testFilename);
    CPPUNIT_TEST(testRead);

    CPPUNIT_TEST_SUITE_END();

    // PUBLIC METHODS /////////////////////////////////////////////////////
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

    /// Test read().
    void testRead(void);

    /** Get test data.
     *
     * @returns Test data.
     */
    TestMeshIO_Data* _getData(void);

    // PROTECTED METHODS ////////////////////////////////////////////////
protected:

    MeshIOLagrit* _io; ///< Test subject.
    TestMeshIOLagrit_Data* _data; ///< Data for tests.

}; // class TestMeshIOLagrit

// ======================================================================
class pylith::meshio::TestMeshIOLagrit_Data : public TestMeshIO_Data {

    // PUBLIC METHODS ///////////////////////////////////////////////////
public:

    /// Constructor
    TestMeshIOLagrit_Data(void);

    /// Destructor
    ~TestMeshIOLagrit_Data(void);

    // PUBLIC MEMBERS ///////////////////////////////////////////////////
public:

    const char* filenameGmv;
    const char* filenamePset;
    bool ioInt32;
    bool isRecordHeader32Bit;

}; // class TestMeshIOLagrit_Data


#endif // pylith_meshio_testmeshiolagrit_hh


// End of file
