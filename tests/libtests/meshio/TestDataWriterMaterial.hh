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
 * @file tests/libtests/meshio/TestDataWriterMeshaterial
 *
 * @brief C++ TestDataWriterMesh object
 *
 * C++ unit testing for DataWriter<Mesh>.
 */

#if !defined(pylith_meshio_testdatawritermaterial_hh)
#define pylith_meshio_testdatawritermaterial_hh

#include "TestDataWriter.hh" // USES TestDataWriter_Data

#include "pylith/topology/topologyfwd.hh" // HOLDSA Mesh, USES Field

namespace pylith {
    namespace meshio {
        class TestDataWriterMaterial;
        class TestDataWriterMaterial_Data;
    } // meshio
} // pylith

// ------------------------------------------------------------------------------------------------
class pylith::meshio::TestDataWriterMaterial {
    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    /// Constructor.
    TestDataWriterMaterial(void);

    /// Destructor.
    ~TestDataWriterMaterial(void);

    /// Set data for tri test case.
    static
    void setDataTri(TestDataWriterMaterial_Data* data);

    /// Set data for quad test case.
    static
    void setDataQuad(TestDataWriterMaterial_Data* data);

    /// Set data for tet test case.
    static
    void setDataTet(TestDataWriterMaterial_Data* data);

    /// Set data for hex test case.
    static
    void setDataHex(TestDataWriterMaterial_Data* data);

    // PROTECTED METHODS //////////////////////////////////////////////////////////////////////////
protected:

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

    /** Get test data.
     *
     * @returns Test data.
     */
    virtual
    TestDataWriterMaterial_Data* _getData(void) = 0;

    // PROTECTED MEMBERS //////////////////////////////////////////////////////////////////////////
protected:

    pylith::topology::Mesh* _domainMesh; ///< Finite-element mesh.
    pylith::topology::Mesh* _materialMesh; ///< Subdomain mesh over material.

}; // class TestDataWriterMaterial

// ------------------------------------------------------------------------------------------------
class pylith::meshio::TestDataWriterMaterial_Data : public TestDataWriter_Data {
    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    /// Constructor
    TestDataWriterMaterial_Data(void);

    /// Destructor
    ~TestDataWriterMaterial_Data(void);

    // PUBLIC MEMBERS /////////////////////////////////////////////////////////////////////////////
public:

    int materialId; ///< Id of material.

}; // class TestDataWriterMaterial_Data

#endif // pylith_meshio_testdatawritermaterial_hh

// End of file
