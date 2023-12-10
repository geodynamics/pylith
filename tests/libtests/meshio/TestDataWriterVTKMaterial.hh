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
#include "TestDataWriterMaterial.hh" // ISA TestDataWriterMaterial

#include "pylith/topology/topologyfwd.hh" // USES Mesh, Field

/// Namespace for pylith package
namespace pylith {
    namespace meshio {
        class TestDataWriterVTKMaterial;
        class TestDataWriterVTKMaterial_Data;
    } // meshio
} // pylith

// ------------------------------------------------------------------------------------------------
class pylith::meshio::TestDataWriterVTKMaterial : public TestDataWriterVTK, public TestDataWriterMaterial {
    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    /// Constructor.
    TestDataWriterVTKMaterial(TestDataWriterVTKMaterial_Data* data);

    /// Destructor.
    ~TestDataWriterVTKMaterial(void);

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
    TestDataWriterMaterial_Data* _getData(void);

    // PROTECTED MEMBDERS /////////////////////////////////////////////////////////////////////////
protected:

    TestDataWriterVTKMaterial_Data* _data; ///< Data for testing.

}; // class TestDataWriterVTKMaterial

// ======================================================================
class pylith::meshio::TestDataWriterVTKMaterial_Data : public TestDataWriterVTK_Data, public TestDataWriterMaterial_Data {};

// End of file
