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
 * @file tests/libtests/meshio/TestDataWriterPoints.hh
 *
 * @brief C++ TestDataWriterPoints object
 *
 * C++ unit testing for DataWriter for points.
 */

#if !defined(pylith_meshio_testdatawriterpoints_hh)
#define pylith_meshio_testdatawriterpoints_hh

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
    void _createVertexFields(topology::Fields* fields) const;

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

    pylith::topology::Mesh* _domainMesh; ///< Mesh for domain
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

    const char* pointsLabel; ///< Label marking points.

}; // class TestDataWriterPoints_Data

#endif // pylith_meshio_testdatawriterpoints_hh

// End of file
