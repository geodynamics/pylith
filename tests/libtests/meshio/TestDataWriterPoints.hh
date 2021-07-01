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
 * @file tests/libtests/meshio/TestDataWriterPoints.hh
 *
 * @brief C++ TestDataWriterPoints object
 *
 * C++ unit testing for DataWriter for points.
 */

#if !defined(pylith_meshio_testdatawriterpoints_hh)
#define pylith_meshio_testdatawriterpoints_hh

#include "TestDataWriter.hh" // USES TestDataWriter_Data

#include "pylith/topology/topologyfwd.hh" // USES Mesh, Field

namespace pylith {
    namespace meshio {
        class TestDataWriterPoints;

        class TestDataWriterPoints_Data;
    } // meshio
} // pylith

class pylith::meshio::TestDataWriterPoints {
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

#endif // pylith_meshio_testdatawriterpoints_hh

// End of file
