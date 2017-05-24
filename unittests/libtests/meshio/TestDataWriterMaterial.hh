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
 * @file unittests/libtests/meshio/TestDataWriterMeshaterial
 *
 * @brief C++ TestDataWriterMesh object
 *
 * C++ unit testing for DataWriter<Mesh>.
 */

#if !defined(pylith_meshio_testdatawritermaterial_hh)
#define pylith_meshio_testdatawritermaterial_hh

#include "TestDataWriter.hh" // USES TestDataWriter_Data

#include "pylith/topology/topologyfwd.hh" // HOLDSA Mesh, USES Fields

/// Namespace for pylith package
namespace pylith {
    namespace meshio {
        class TestDataWriterMaterial;

        class TestDataWriterMaterial_Data;
    } // meshio
} // pylith

// ======================================================================
class pylith::meshio::TestDataWriterMaterial {

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
    TestDataWriterMaterial_Data* _getData(void) = 0;

    // PROTECTED MEMBERS //////////////////////////////////////////////////
protected:

    pylith::topology::Mesh* _mesh; ///< Finite-element mesh.

}; // class TestDataWriterMaterial


// ======================================================================
class pylith::meshio::TestDataWriterMaterial_Data : public TestDataWriter_Data {

    // PUBLIC METHODS /////////////////////////////////////////////////////
public:

    /// Constructor
    TestDataWriterMaterial_Data(void);

    /// Destructor
    ~TestDataWriterMaterial_Data(void);

    // PUBLIC MEMBERS ///////////////////////////////////////////////////
public:

    int materialId; ///< Id of material.

}; // class TestDataWriterMaterial_Data


#endif // pylith_meshio_testdatawritermaterial_hh


// End of file
