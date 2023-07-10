// -*- C++ -*-
//
// ======================================================================
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
// ======================================================================
//

/**
 * @file libsrc/meshio/MeshIO.hh
 *
 * @brief C++ base class for managing mesh input/output.
 */

#if !defined(pylith_meshio_meshio_hh)
#define pylith_meshio_meshio_hh

#include "meshiofwd.hh" // forward declarations

#include "pylith/utils/PyreComponent.hh" // ISA PyreComponent

#include "pylith/meshio/MeshBuilder.hh" // USES MeshBuilder::GroupPtType

#include "pylith/topology/topologyfwd.hh" // forward declarations
#include "spatialdata/units/unitsfwd.hh" // forward declarations
#include "pylith/utils/arrayfwd.hh" // USES scalar_array, int_array, string_vector

// ------------------------------------------------------------------------------------------------
class pylith::meshio::MeshIO : public pylith::utils::PyreComponent {
    // PUBLIC MEMBERS /////////////////////////////////////////////////////////////////////////////
public:

    /// Constructor
    MeshIO(void);

    /// Destructor
    virtual ~MeshIO(void);

    /// Deallocate PETSc and local data structures.
    virtual
    void deallocate(void);

    /** Read mesh from file.
     *
     * @param[in] mesh PyLith finite-element mesh.
     * @param[in] check Check topology of mesh.
     */
    void read(pylith::topology::Mesh* mesh,
              const bool check=true);

    /** Write mesh to file.
     *
     * @param mesh PyLith finite-element mesh.
     */
    void write(pylith::topology::Mesh* const mesh);

    // PROTECTED MEMBERS //////////////////////////////////////////////////////////////////////////
protected:

    /// Write mesh
    virtual
    void _write(void) const = 0;

    /// Read mesh
    virtual
    void _read(void) = 0;

    /** Get spatial dimension of mesh.
     *
     * @returns Spatial dimension of mesh
     */
    int getMeshDim(void) const;

    /** Get information about vertices in mesh.
     *
     * Method caller is responsible for memory management.
     *
     * @param coordinates Pointer to array of vertex coordinates
     * @param numVertices Pointer to number of vertices
     * @param spaceDim Poiner to dimension of vector space for coordinates
     */
    void _getVertices(scalar_array* coordinates,
                      int* numVertices,
                      int* spaceDim) const;

    /** Get information about cells in mesh.
     *
     * The indices of the vertices will be returned using zero based
     * indices. In other words, the first vertex will be referred to
     * using an index of 0.
     *
     * Method caller is responsible for memory management.
     *
     * @param cells Pointer to array of indicates of vertices in each cell
     * @param numCells Pointer to number of cells in mesh
     * @param numCorners Pointer to number of vertices in each cell
     * @param meshDim Pointer to number of dimensions associated with cell
     */
    void _getCells(int_array* cells,
                   int* numCells,
                   int* numCorners,
                   int* meshDim) const;

    /** Tag cells in mesh with material identifiers.
     *
     * @param materialIds Material identifiers [numCells]
     */
    void _setMaterials(const int_array& materialIds);

    /** Get material identifiers for cells.
     *
     * @param materialIds Material identifiers [numCells]
     */
    void _getMaterials(int_array* pMaterialIds) const;

    /** Get names of all groups in mesh.
     *
     * @returns Array of group names.
     */
    void _getGroupNames(string_vector* names) const;

    /** Return a point group
     *
     * The indices of the vertices will be returned using zero based
     * indices. In other words, the first vertex will be referred to
     * using an index of 0.
     *
     * @param points An array of the points in the group
     * @param type The point type, e.g. VERTEX, CELL
     * @param name The group name
     */
    void _getGroup(int_array* points,
                   pylith::meshio::MeshBuilder::GroupPtType* type,
                   const char *name) const;

    // PROTECTED MEMBERS //////////////////////////////////////////////////////////////////////////
protected:

    pylith::topology::Mesh* _mesh; ///< Pointer to finite-element mesh.

}; // MeshIO

#endif // pylith_meshio_meshio_hh

// End of file
