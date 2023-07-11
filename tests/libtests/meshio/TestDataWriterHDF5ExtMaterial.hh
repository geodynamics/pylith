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
