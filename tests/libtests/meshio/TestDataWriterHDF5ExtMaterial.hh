// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2023, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

/**
 * @file tests/libtests/meshio/TestDataWriterHDF5ExtMaterial.hh
 *
 * @brief C++ TestDataWriterHDF5ExtMaterial object
 *
 * C++ unit testing for DataWriterHDF5ExtMaterial.
 */

#if !defined(pylith_meshio_testdatawriterhdf5extmaterial_hh)
#define pylith_meshio_testdatawriterhdf5extmaterial_hh

#include "TestDataWriterHDF5.hh" // ISA TestDataWriterHDF5
#include "TestDataWriterMaterial.hh" // ISA TestDataWriterMaterial

#include "pylith/topology/topologyfwd.hh" // USES Mesh, Field

/// Namespace for pylith package
namespace pylith {
    namespace meshio {
        class TestDataWriterHDF5ExtMaterial;
        class TestDataWriterHDF5ExtMaterial_Data;
    } // meshio
} // pylith

// ------------------------------------------------------------------------------------------------
/// C++ unit testing for DataWriterHDF5Ext
class pylith::meshio::TestDataWriterHDF5ExtMaterial : public TestDataWriterHDF5, public TestDataWriterMaterial {
    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    /// Constructor.
    TestDataWriterHDF5ExtMaterial(TestDataWriterHDF5ExtMaterial_Data* data);

    /// Destructor.
    ~TestDataWriterHDF5ExtMaterial(void);

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

    TestDataWriterHDF5ExtMaterial_Data* _data; ///< Data for testing.

}; // class TestDataWriterHDF5ExtMaterial

// ------------------------------------------------------------------------------------------------
class pylith::meshio::TestDataWriterHDF5ExtMaterial_Data : public TestDataWriterHDF5_Data, public TestDataWriterMaterial_Data {};

#endif // pylith_meshio_testdatawriterhdf5extmaterial_hh

// End of file
