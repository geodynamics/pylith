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
#include "TestDataWriterMesh.hh" // ISA TestDataWriterMesh

#include "pylith/topology/topologyfwd.hh" // USES Mesh, Field

/// Namespace for pylith package
namespace pylith {
    namespace meshio {
        class TestDataWriterHDF5Mesh;
        class TestDataWriterHDF5Mesh_Data;
    } // meshio
} // pylith

// ------------------------------------------------------------------------------------------------
class pylith::meshio::TestDataWriterHDF5Mesh : public TestDataWriterHDF5, public TestDataWriterMesh {
    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    /// Constructor.
    TestDataWriterHDF5Mesh(TestDataWriterHDF5Mesh_Data* data);

    /// Destructor.
    ~TestDataWriterHDF5Mesh(void);

    /// Test filename()
    static
    void testAccessors(void);

    /// Test hdf5Filename.
    static
    void testHdf5Filename(void);

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
    TestDataWriter_Data* _getData(void);

    // PROTECTED MEMBDERS /////////////////////////////////////////////////////////////////////////
protected:

    TestDataWriterHDF5Mesh_Data* _data; ///< Data for testing.

}; // class TestDataWriterHDF5Mesh

// ------------------------------------------------------------------------------------------------
class pylith::meshio::TestDataWriterHDF5Mesh_Data : public TestDataWriterHDF5_Data, public TestDataWriter_Data {};

// End of file
