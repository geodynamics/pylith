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
 * @file tests/libtests/meshio/TestMeshIOCubit.hh
 *
 * @brief C++ TestMeshIOCubit object
 *
 * C++ unit testing for MeshIOCubit.
 */

#if !defined(pylith_meshio_testmeshiocubit_hh)
#define pylith_meshio_testmeshiocubit_hh

// Include directives ---------------------------------------------------
#include "TestMeshIO.hh"

// Forward declarations -------------------------------------------------
namespace pylith {
    namespace meshio {
        class TestMeshIOCubit;

        class TestMeshIOCubit_Data; // test data
    } // meshio
} // pylith

// TestMeshIOCubit ------------------------------------------------------
class pylith::meshio::TestMeshIOCubit : public TestMeshIO { // class TestMeshIOCubit

    // CPPUNIT TEST SUITE /////////////////////////////////////////////////
    CPPUNIT_TEST_SUITE(TestMeshIOCubit);

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

    MeshIOCubit* _io; ///< Test subject.
    TestMeshIOCubit_Data* _data; ///< Data for tests.

}; // class TestMeshIOCubit

// ======================================================================
class pylith::meshio::TestMeshIOCubit_Data : public TestMeshIO_Data {

    // PUBLIC METHODS ///////////////////////////////////////////////////
public:

    /// Constructor
    TestMeshIOCubit_Data(void);

    /// Destructor
    ~TestMeshIOCubit_Data(void);

    // PUBLIC MEMBERS ///////////////////////////////////////////////////
public:

    const char* filename;

}; // class TestMeshIOCubit_Data


#endif // pylith_meshio_testmeshiocubit_hh


// End of file
