// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2024, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================
#pragma once

#include "TestDataWriterHDF5.hh" // ISA TestDataWriterHDF5
#include "TestDataWriterSubmesh.hh" // ISA TestDataWriterSubmesh

#include "pylith/topology/topologyfwd.hh" // USES Mesh, Field

namespace pylith {
    namespace meshio {
        class TestDataWriterHDF5ExtSubmesh;
        class TestDataWriterHDF5ExtSubmesh_Data;
    } // meshio
} // pylith

// ------------------------------------------------------------------------------------------------
class pylith::meshio::TestDataWriterHDF5ExtSubmesh : public TestDataWriterHDF5, public TestDataWriterSubmesh {
    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    /// Constructor.
    TestDataWriterHDF5ExtSubmesh(TestDataWriterHDF5ExtSubmesh_Data* data);

    /// Destructor.
    ~TestDataWriterHDF5ExtSubmesh(void);

    /// Test open() and close()
    void testOpenClose(void);

    /// Test writeVertexField.
    void testWriteVertexField(void);

    /// Test writeCellField.
    void testWriteCellField(void);

    // PROTECTED METHODS //////////////////////////////////////////////////////////////////////////
protected:

    /** Get test data.
     *
     * @returns Test data.
     */
    TestDataWriterSubmesh_Data* _getData(void);

    // PROTECTED MEMBDERS /////////////////////////////////////////////////////////////////////////
protected:

    TestDataWriterHDF5ExtSubmesh_Data* _data; ///< Data for testing.

}; // class TestDataWriterHDF5ExtSubmesh

// ======================================================================
class pylith::meshio::TestDataWriterHDF5ExtSubmesh_Data :
    public TestDataWriterHDF5_Data,
    public TestDataWriterSubmesh_Data {};

// End of file
