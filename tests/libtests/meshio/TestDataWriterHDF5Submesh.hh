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
 * @file tests/libtests/meshio/TestDataWriterHDF5Submesh.hh
 *
 * @brief C++ TestDataWriterHDF5Submesh object
 *
 * C++ unit testing for DataWriterHDF5Submesh.
 */

#if !defined(pylith_meshio_testdatawriterhdf5submesh_hh)
#define pylith_meshio_testdatawriterhdf5submesh_hh

#include "TestDataWriterHDF5.hh"
#include "TestDataWriterSubmesh.hh"

#include "pylith/topology/topologyfwd.hh" // USES Mesh, Field

namespace pylith {
    namespace meshio {
        class TestDataWriterHDF5Submesh;
        class TestDataWriterHDF5Submesh_Data;
    } // meshio
} // pylith

// ------------------------------------------------------------------------------------------------
class pylith::meshio::TestDataWriterHDF5Submesh :
    public TestDataWriterHDF5, public TestDataWriterSubmesh {
    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    /// Constructor.
    TestDataWriterHDF5Submesh(TestDataWriterHDF5Submesh_Data* data);

    /// Destructor.
    ~TestDataWriterHDF5Submesh(void);

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

    TestDataWriterHDF5Submesh_Data* _data; ///< Data for testing.

}; // class TestDataWriterVTKSubmesh

// ======================================================================
class pylith::meshio::TestDataWriterHDF5Submesh_Data :
    public TestDataWriterHDF5_Data,
    public TestDataWriterSubmesh_Data {};

#endif // pylith_meshio_testdatawriterhdf5submesh_hh

// End of file
