// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2017 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

/**
 * @file unittests/libtests/meshio/TestDataWriterSubMesh.hh
 *
 * @brief C++ TestDataWriterSubMesh object
 *
 * C++ unit testing for DataWriter<SubMesh>.
 */

#if !defined(pylith_meshio_testdatawritersubmesh_hh)
#define pylith_meshio_testdatawritersubmesh_hh

#include "TestDataWriter.hh" // USES TestDataWriter_Data

#include "pylith/topology/topologyfwd.hh" // USES Mesh, Field

/// Namespace for pylith package
namespace pylith {
    namespace meshio {
        class TestDataWriterSubMesh;

        class TestDataWriterSubMesh_Data;
    } // meshio
} // pylith

/// C++ unit testing for DataWriter<SubMesh>.
class pylith::meshio::TestDataWriterSubMesh {

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
    void _createVertexFields(pylith::topology::Fields* fields);

    /** Create cell fields.
     *
     * @param fields Cell fields.
     */
    void _createCellFields(pylith::topology::Fields* fields);

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
    TestDataWriterSubMesh_Data* _getData(void) = 0;

    // PROTECTED MEMBERS //////////////////////////////////////////////////
protected:

    pylith::topology::Mesh* _mesh; ///< Mesh for domain
    pylith::topology::Mesh* _submesh; ///< Mesh for subdomain.

}; // class TestDataWriterSubMesh

// ======================================================================
class pylith::meshio::TestDataWriterSubMesh_Data : public TestDataWriter_Data {

    // PUBLIC METHODS /////////////////////////////////////////////////////
public:

    /// Constructor
    TestDataWriterSubMesh_Data(void);

    /// Destructor
    ~TestDataWriterSubMesh_Data(void);

    // PUBLIC MEMBERS ///////////////////////////////////////////////////
public:

    const char* bcLabel; ///< Label marking submesh.

}; // class TestDataWriterSubMesh_Data


#endif // pylith_meshio_testdatawritersubmesh_hh


// End of file
