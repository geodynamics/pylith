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
// Copyright (c) 2010-2022 University of California, Davis
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

#include "TestMeshIO.hh" // ISA TestMeshIO

namespace pylith {
    namespace meshio {
        class TestMeshIOAscii;
    } // meshio
} // pylith

// ======================================================================
class pylith::meshio::TestMeshIOAscii : public TestMeshIO {
    // PUBLIC METHODS ///////////////////////////////////////////////////
public:

    /// Constructor.
    TestMeshIOAscii(TestMeshIO_Data* data);

    /// Destructor.
    ~TestMeshIOAscii(void);

    /// Test filename()
    void testFilename(void);

    /// Test write() and read().
    void testWriteRead(void);

    /// Test read().
    void testRead(void);

    // PROTECTED METHODS ////////////////////////////////////////////////
protected:

    MeshIOAscii* _io; ///< Test subject.

}; // class TestMeshIOAscii

#endif // pylith_meshio_testmeshioascii_hh

// End of file
