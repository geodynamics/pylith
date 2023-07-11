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
 * @file tests/libtests/meshio/TestMeshIOCubit.hh
 *
 * @brief C++ TestMeshIOCubit object
 *
 * C++ unit testing for MeshIOCubit.
 */

#if !defined(pylith_meshio_testmeshiocubit_hh)
#define pylith_meshio_testmeshiocubit_hh

#include "TestMeshIO.hh"

namespace pylith {
    namespace meshio {
        class TestMeshIOCubit;
    } // meshio
} // pylith

// ------------------------------------------------------------------------------------------------
class pylith::meshio::TestMeshIOCubit : public TestMeshIO {
    // PUBLIC METHODS /////////////////////////////////////////////////////
public:

    /// Constructor.
    TestMeshIOCubit(TestMeshIO_Data* data);

    /// Destructor.
    ~TestMeshIOCubit(void);

    /// Test filename()
    void testFilename(void);

    /// Test read().
    void testRead(void);

    // PROTECTED METHODS ////////////////////////////////////////////////
protected:

    MeshIOCubit* _io; ///< Test subject.

}; // class TestMeshIOCubit

#endif // pylith_meshio_testmeshiocubit_hh

// End of file
