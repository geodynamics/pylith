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
 * @file tests/libtests/meshio/TestMeshIOPetsc.hh
 *
 * @brief C++ TestMeshIOPetsc object
 *
 * C++ unit testing for MeshIOPetsc.
 */

#if !defined(pylith_meshio_testmeshiopetsc_hh)
#define pylith_meshio_testmeshiopetsc_hh

#include "TestMeshIO.hh"

namespace pylith {
    namespace meshio {
        class TestMeshIOPetsc;
    } // meshio
} // pylith

class pylith::meshio::TestMeshIOPetsc : public TestMeshIO {
    // PUBLIC METHODS /////////////////////////////////////////////////////
public:

    /// Constructor.
    TestMeshIOPetsc(TestMeshIO_Data* data);

    /// Destructor.
    ~TestMeshIOPetsc(void);

    /// Test filename()
    void testFilename(void);

    /// Test read().
    void testRead(void);

    // PROTECTED METHODS ////////////////////////////////////////////////
protected:

    MeshIOPetsc* _io; ///< Test subject.

}; // class TestMeshIOPetsc

#endif // pylith_meshio_testmeshiopetsc_hh

// End of file
