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

#include "TestDataWriterVTK.hh"
#include "TestDataWriterSubmesh.hh"

#include "pylith/topology/topologyfwd.hh" // USES Mesh, Field

/// Namespace for pylith package
namespace pylith {
    namespace meshio {
        class TestDataWriterVTKSubmesh;
        class TestDataWriterVTKSubmesh_Data;
    } // meshio
} // pylith

// ------------------------------------------------------------------------------------------------
class pylith::meshio::TestDataWriterVTKSubmesh : public TestDataWriterVTK, public TestDataWriterSubmesh {
    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    /// Constructor.
    TestDataWriterVTKSubmesh(TestDataWriterVTKSubmesh_Data* data);

    /// Destructor.
    ~TestDataWriterVTKSubmesh(void);

    /// Test openTimeStep() and closeTimeStep()
    void testTimeStep(void);

    /// Test writeVertexField.
    void testWriteVertexField(void);

    /// Test writeCellField.
    void testWriteCellField(void);

    // PROTECTED METHODS //////////////////////////////////////////////////
protected:

    /** Get test data.
     *
     * @returns Test data.
     */
    TestDataWriterSubmesh_Data* _getData(void);

    // PROTECTED MEMBDERS /////////////////////////////////////////////////
protected:

    TestDataWriterVTKSubmesh_Data* _data; ///< Data for testing.

}; // class TestDataWriterVTKSubmesh

// ======================================================================
class pylith::meshio::TestDataWriterVTKSubmesh_Data : public TestDataWriterVTK_Data, public TestDataWriterSubmesh_Data {};

// End of file
