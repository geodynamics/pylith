// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2023, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================
#pragma once

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

// End of file
