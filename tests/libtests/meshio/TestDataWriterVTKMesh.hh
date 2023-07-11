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
 * @file tests/libtests/meshio/TestDataWriterVTKMesh.hh
 *
 * @brief C++ TestDataWriterVTKMesh object
 *
 * C++ unit testing for DataWriterVTKMesh.
 */

#if !defined(pylith_meshio_testdatawritervtkmesh_hh)
#define pylith_meshio_testdatawritervtkmesh_hh

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

#endif // pylith_meshio_testdatawritervtkmesh_hh

// End of file
