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
 * @file tests/libtests/meshio/TestMeshIO.hh
 *
 * @brief C++ TestMeshIO object
 *
 * C++ unit testing for MeshIO.
 */

#if !defined(pylith_meshio_testmeshio_hh)
#define pylith_meshio_testmeshio_hh

// Include directives ---------------------------------------------------
#include <cppunit/extensions/HelperMacros.h>

#include "pylith/topology/topologyfwd.hh" // USES Mesh
#include "pylith/meshio/meshiofwd.hh" // USES MeshIO

// Forward declarations -------------------------------------------------
/// Namespace for pylith package
namespace pylith {
    namespace meshio {
        class TestMeshIO;

        class TestMeshIO_Data; // test data
    } // meshio
} // pylith

// ----------------------------------------------------------------------
class pylith::meshio::TestMeshIO : public CppUnit::TestFixture {

    // PUBLIC METHODS ///////////////////////////////////////////////////////
public:

    /// Setup testing data.
    void setUp(void);

    /// Tear down testing data.
    void tearDown(void);

    // PROTECTED METHODS ////////////////////////////////////////////////////
protected:

    /** Get test data.
     *
     * @returns Pointer to test data.
     */
    virtual
    TestMeshIO_Data* _getData(void) = 0;

    /// Get simple mesh for testing I/O.
    void _createMesh(void);

    /// Check values in mesh against data.
    void _checkVals(void);

    /** Test debug().
     *
     * @param iohandler MeshIO object.
     */
    void _testDebug(MeshIO& iohandler);

    // PROTECTED MEMBERS ////////////////////////////////////////////////////
protected:

    pylith::topology::Mesh* _mesh; ///< Finite-element mesh.


}; // class TestMeshIO


// ----------------------------------------------------------------------
class pylith::meshio::TestMeshIO_Data {

    // PUBLIC METHODS ///////////////////////////////////////////////////////
public:

    /// Constructor
    TestMeshIO_Data(void);

    /// Destructor
    ~TestMeshIO_Data(void);

    // PUBLIC MEMBERS ///////////////////////////////////////////////////////
public:

    PylithInt numVertices; ///< Number of vertices
    PylithInt spaceDim; ///< Number of dimensions in vertex coordinates
    PylithInt numCells; ///< Number of cells
    PylithInt cellDim; ///< Number of dimensions associated with cell
    PylithInt numCorners; ///< Number of vertices in cell

    PylithScalar* vertices; ///< Pointer to coordinates of vertices
    PylithInt* cells; ///< Pointer to indices of vertices in cells
    PylithInt* materialIds; ///< Pointer to cell material identifiers

    PylithInt* groups; ///< Array of pointers to indices of points in groups
    PylithInt* groupSizes; ///< Array of sizes of each group
    char** groupNames; ///< Array of group names
    char** groupTypes; ///< Array of group types
    PylithInt numGroups; ///< Number of groups

    bool useIndexZero; ///< Indices start with 0 if true, 1 if false

};


#endif // pylith_meshio_testmeshio_hh


// End of file
