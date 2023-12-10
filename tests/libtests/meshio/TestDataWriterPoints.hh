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

#include "pylith/topology/topologyfwd.hh" // USES Mesh, Field

namespace pylith {
    namespace meshio {
        class TestDataWriterPoints;
        class TestDataWriterPoints_Data;
    } // meshio
} // pylith

// ------------------------------------------------------------------------------------------------
class pylith::meshio::TestDataWriterPoints {
    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    /// Constructor.
    TestDataWriterPoints(void);

    /// Destructor.
    ~TestDataWriterPoints(void);

    /// Set data for tri test case.
    static
    void setDataTri(TestDataWriterPoints_Data* data);

    /// Set data for quad test case.
    static
    void setDataQuad(TestDataWriterPoints_Data* data);

    /// Set data for tet test case.
    static
    void setDataTet(TestDataWriterPoints_Data* data);

    /// Set data for hex test case.
    static
    void setDataHex(TestDataWriterPoints_Data* data);

    // PROTECTED METHODS //////////////////////////////////////////////////////////////////////////
protected:

    /// Initialize mesh.
    void _initialize(void);

    /** Create vertex fields.
     *
     * @param fields Vertex fields.
     */
    void _createVertexField(pylith::topology::Field* field);

    /** Get test data.
     *
     * @returns Test data.
     */
    virtual
    TestDataWriterPoints_Data* _getData(void) = 0;

    // PROTECTED MEMBERS //////////////////////////////////////////////////
protected:

    pylith::topology::Mesh* _pointMesh; ///< Mesh associated with point data.

}; // class TestDataWriterPoints

// ======================================================================
class pylith::meshio::TestDataWriterPoints_Data : public TestDataWriter_Data {
    // PUBLIC METHODS /////////////////////////////////////////////////////
public:

    /// Constructor
    TestDataWriterPoints_Data(void);

    /// Destructor
    ~TestDataWriterPoints_Data(void);

    // PUBLIC MEMBERS ///////////////////////////////////////////////////
public:

    int numPoints; ///< Number of points for interpolation.
    PylithReal* points; /// Points for interpolation.
    pylith::string_vector names; ///< Station names for points.

}; // class TestDataWriterPoints_Data

// End of file
