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

#include "TestDataWriterHDF5.hh" // ISA TestDataWriterHDF5
#include "TestDataWriterMaterial.hh" // ISA TestDataWriterMaterial

#include "pylith/topology/topologyfwd.hh" // USES Mesh, Field

/// Namespace for pylith package
namespace pylith {
    namespace meshio {
        class TestDataWriterHDF5Material;
        class TestDataWriterHDF5Material_Data;
    } // meshio
} // pylith

// ------------------------------------------------------------------------------------------------
class pylith::meshio::TestDataWriterHDF5Material : public TestDataWriterHDF5, public TestDataWriterMaterial {
    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    /// Constructor.
    TestDataWriterHDF5Material(TestDataWriterHDF5Material_Data* data);

    /// Destructor.
    ~TestDataWriterHDF5Material(void);

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
    TestDataWriterMaterial_Data* _getData(void);

    // PROTECTED MEMBDERS /////////////////////////////////////////////////////////////////////////
protected:

    TestDataWriterHDF5Material_Data* _data; ///< Data for testing.

}; // class TestDataWriterHDF5Material

// ------------------------------------------------------------------------------------------------
class pylith::meshio::TestDataWriterHDF5Material_Data : public TestDataWriterHDF5_Data, public TestDataWriterMaterial_Data {};

// End of file
