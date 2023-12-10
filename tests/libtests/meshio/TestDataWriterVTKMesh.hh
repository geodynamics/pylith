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

#include "TestDataWriterVTK.hh" // ISA TestDataWriterVTK
#include "TestDataWriterMesh.hh" // ISA TestDataWriterMesh

#include "pylith/topology/topologyfwd.hh" // USES Mesh, Field

/// Namespace for pylith package
namespace pylith {
    namespace meshio {
        class TestDataWriterVTKMesh;
        class TestDataWriterVTKMesh_Data;
    } // meshio
} // pylith

// ------------------------------------------------------------------------------------------------
class pylith::meshio::TestDataWriterVTKMesh : public TestDataWriterVTK, public TestDataWriterMesh {
    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    /// Constructor.
    TestDataWriterVTKMesh(TestDataWriterVTKMesh_Data* data);

    /// Destructor.
    ~TestDataWriterVTKMesh(void);

    /// Test filename(), timeFormat(), timeConstant(), precision()
    static
    void testAccessors(void);

    /// Test vtkFilename.
    static
    void testVtkFilename(void);

    /// Test openTimeStep() and closeTimeStep()
    void testTimeStep(void);

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

    // PROTECTED MEMBERS //////////////////////////////////////////////////////////////////////////
protected:

    TestDataWriterVTKMesh_Data* _data; ///< Data for testing.

}; // class TestDataWriterVTKMesh

// ------------------------------------------------------------------------------------------------
class pylith::meshio::TestDataWriterVTKMesh_Data : public TestDataWriterVTK_Data, public TestDataWriter_Data {};

// End of file
