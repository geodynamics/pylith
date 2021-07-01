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
// Copyright (c) 2010-2021 University of California, Davis
//
// See LICENSE.md for license information.
//
// ----------------------------------------------------------------------
//

/**
 * @file tests/libtests/meshio/TestDataWriterSubmesh.hh
 *
 * @brief C++ TestDataWriterSubmesh object
 *
 * C++ unit testing for DataWriter for submesh.
 */

#if !defined(pylith_meshio_testdatawritersubmesh_hh)
#define pylith_meshio_testdatawritersubmesh_hh

#include "TestDataWriter.hh" // USES TestDataWriter_Data

#include "pylith/topology/topologyfwd.hh" // USES Mesh, Field

/// Namespace for pylith package
namespace pylith {
    namespace meshio {
        class TestDataWriterSubmesh;

        class TestDataWriterSubmesh_Data;
    } // meshio
} // pylith

class pylith::meshio::TestDataWriterSubmesh {
    // PROTECTED METHODS //////////////////////////////////////////////////
protected:

    /// Setup testing data.
    void setUp(void);

    /// Tear down testing data.
    void tearDown(void);

    /// Initialize mesh.
    void _initialize(void);

    /** Create vertex fields.
     *
     * @param fields Vertex fields.
     */
    void _createVertexField(pylith::topology::Field* field);

    /** Create cell fields.
     *
     * @param fields Cell fields.
     */
    void _createCellField(pylith::topology::Field* field);

    /// Set data for tri test case.
    void _setDataTri(void);

    /// Set data for quad test case.
    void _setDataQuad(void);

    /// Set data for tet test case.
    void _setDataTet(void);

    /// Set data for hex test case.
    void _setDataHex(void);

    /** Get test data.
     *
     * @returns Test data.
     */
    virtual
    TestDataWriterSubmesh_Data* _getData(void) = 0;

    // PROTECTED MEMBERS //////////////////////////////////////////////////
protected:

    pylith::topology::Mesh* _mesh; ///< Mesh for domain
    pylith::topology::Mesh* _submesh; ///< Mesh for subdomain.

}; // class TestDataWriterSubmesh

// ======================================================================
class pylith::meshio::TestDataWriterSubmesh_Data : public TestDataWriter_Data {
    // PUBLIC METHODS /////////////////////////////////////////////////////
public:

    /// Constructor
    TestDataWriterSubmesh_Data(void);

    /// Destructor
    ~TestDataWriterSubmesh_Data(void);

    // PUBLIC MEMBERS ///////////////////////////////////////////////////
public:

    const char* bcLabel; ///< Label marking submesh.

}; // class TestDataWriterSubmesh_Data

#endif // pylith_meshio_testdatawritersubmesh_hh

// End of file
