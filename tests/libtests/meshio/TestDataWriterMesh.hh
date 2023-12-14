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

#include "TestDataWriter.hh" // USES TestDataWriter_Data

#include "pylith/topology/topologyfwd.hh" // HOLDSA Mesh, USES Fields

/// Namespace for pylith package
namespace pylith {
    namespace meshio {
        class TestDataWriterMesh;
    } // meshio
} // pylith

// ------------------------------------------------------------------------------------------------
class pylith::meshio::TestDataWriterMesh {
    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    /// Constructor
    TestDataWriterMesh(void);

    /// Tear down testing data.
    ~TestDataWriterMesh(void);

    /// Set data for tri test case.
    static
    void setDataTri(TestDataWriter_Data* data);

    /// Set data for quad test case.
    static
    void setDataQuad(TestDataWriter_Data* data);

    /// Set data for tet test case.
    static
    void setDataTet(TestDataWriter_Data* data);

    /// Set data for hex test case.
    static
    void setDataHex(TestDataWriter_Data* data);

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
    TestDataWriter_Data* _getData(void) = 0;

    // PROTECTED MEMBERS //////////////////////////////////////////////////
protected:

    pylith::topology::Mesh* _mesh; ///< Finite-element mesh.

}; // class TestDataWriterMesh

// End of file
