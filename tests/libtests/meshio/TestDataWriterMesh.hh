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
 * @file tests/libtests/meshio/TestDataWriterMesh.hh
 *
 * @brief C++ TestDataWriterMesh object
 *
 * C++ unit testing for DataWriter<Mesh>.
 */

#if !defined(pylith_meshio_testdatawritermesh_hh)
#define pylith_meshio_testdatawritermesh_hh

#include "TestDataWriter.hh" // USES TestDataWriter_Data

#include "pylith/topology/topologyfwd.hh" // HOLDSA Mesh, USES Fields

/// Namespace for pylith package
namespace pylith {
    namespace meshio {
        class TestDataWriterMesh;
    } // meshio
} // pylith

// ======================================================================
class pylith::meshio::TestDataWriterMesh {
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
    TestDataWriter_Data* _getData(void) = 0;

    // PROTECTED MEMBERS //////////////////////////////////////////////////
protected:

    pylith::topology::Mesh* _mesh; ///< Finite-element mesh.

}; // class TestDataWriterMesh

#endif // pylith_meshio_testdatawritermesh_hh

// End of file
