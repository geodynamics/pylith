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
// Copyright (c) 2010-2022 University of California, Davis
//
// See LICENSE.md for license information.
//
// ----------------------------------------------------------------------
//

/**
 * @file tests/libtests/topology/TestRefineUniform.hh
 *
 * @brief C++ TestRefineUniform object
 *
 * C++ unit testing for RefineUniform.
 */

#if !defined(pylith_topology_testrefineuniform_hh)
#define pylith_topology_testrefineuniform_hh

#include "pylith/utils/GenericComponent.hh" // ISA GenericComponent

#include "pylith/topology/topologyfwd.hh" // USES Mesh

namespace pylith {
    namespace topology {
        class TestRefineUniform;
        class TestRefineUniform_Data; // test data
    } // topology
} // pylith

// ------------------------------------------------------------------------------------------------
class pylith::topology::TestRefineUniform : public pylith::utils::GenericComponent {
    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    /// Constructor.
    TestRefineUniform(TestRefineUniform_Data* data);

    /// Destructor.
    ~TestRefineUniform(void);

    /// Test refine().
    void testRefine(void);

    // PRIVATE METHODS ////////////////////////////////////////////////////////////////////////////
private:

    /** Setup mesh.
     *
     * @mesh Mesh to setup.
     */
    void _initializeMesh(Mesh* const mesh);

    TestRefineUniform_Data* _data; ///< Data for testing.

}; // class TestRefineUniform

// ------------------------------------------------------------------------------------------------
class pylith::topology::TestRefineUniform_Data {
    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    /// Constructor
    TestRefineUniform_Data(void);

    /// Destructor
    ~TestRefineUniform_Data(void);

    // PUBLIC MEMBERS /////////////////////////////////////////////////////////////////////////////
public:

    /// @defgroup Input information
    /// @{
    const char* filename; ///< Filename of mesh file.
    int refineLevel; ///< Refinement level.
    const char* faultA; ///< Vertex group associated with fault A (NULL if no fault).
    const char* faultB; ///< Vertex group associated with fault B (NULL if no fault).
    bool isSimplexMesh; ///< True if simplex mesh.
    /// @}

    /// @defgroup Output information
    /// @{
    int numVertices; ///< Number of vertices
    int spaceDim; ///< Number of dimensions in vertex coordinates
    int cellDim; ///< Number of dimensions associated with cell.
    int numCells; ///< Number of cells
    int numCorners; ///< Number of vertices in cell.
    int numCellsCohesive; ///< Number of cohesive cells.
    int numCornersCohesive; ///< Number of vertices in cohesive cell.

    int matIdSum; ///< Sum of material id as simple checksum.

    int* groupSizes; ///< Array of sizes of each group
    const char** groupNames; ///< Array of group names
    const char** groupTypes; ///< Array of group types
    int numGroups; ///< Number of groups
    /// @}

};

#endif // pylith_topology_testrefineuniform_hh

// End of file
