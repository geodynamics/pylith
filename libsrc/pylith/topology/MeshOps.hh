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
// Copyright (c) 2010-2021 University of California, Davis
//
// See LICENSE.md for license information.
//
// ======================================================================
//

/**
 * @file libsrc/topology/MeshOps.hh
 *
 * @brief Simple operations on a Mesh object.
 */

#if !defined(pylith_topology_meshops_hh)
#define pylith_topology_meshops_hh

#include "topologyfwd.hh" // forward declarations

#include "pylith/utils/petscfwd.h" // USES PetscDM
#include "pylith/utils/array.hh" // USES int_array

#include "spatialdata/geocoords/geocoordsfwd.hh"
#include "spatialdata/units/unitsfwd.hh"

class pylith::topology::MeshOps {
    friend class TestMeshOps; // unit testing

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /** Create subdomain mesh using label.
     *
     * @param[in] mesh Mesh for domain.
     * @param[in] label Name of label marking subdomain.
     * @param[in] labelValue Value of label marking subdomain.
     * @param[in] descriptiveLabel Descriptive label for subdomain.
     *
     * @returns Mesh for subdomain.
     */
    static
    pylith::topology::Mesh* createSubdomainMesh(const pylith::topology::Mesh& mesh,
                                                const char* label,
                                                const int labelValue,
                                                const char* descriptiveLabel);

    /** Create lower dimension mesh using label.
     *
     * @param[in] mesh Mesh for domain.
     * @param[in] label Label for vertices marking lower dimension domain.
     *
     * @returns Lower dimension mesh.
     */
    static
    pylith::topology::Mesh* createLowerDimMesh(const pylith::topology::Mesh& mesh,
                                               const char* label);

    /** Create 0-dimension mesh from points.
     *
     * @param[in] pointCoords Array of coordinates of points [numPoints*spaceDim].
     * @param[in] numPoints Number of points.
     * @param[in] cs Coordinate system for points.
     * @param[in] lengthScale Length scale for nondimensionalization.
     * @param[in] comm MPI communicator.
     */
    static
    pylith::topology::Mesh* createFromPoints(const PylithReal* points,
                                             const size_t numPoints,
                                             const spatialdata::geocoords::CoordSys* cs,
                                             const PylithReal lengthScale,
                                             MPI_Comm comm);

    /** Nondimensionalize the finite-element mesh.
     *
     * @param[in] mesh Finite-element mesh.
     * @param[in] normalizer Nondimensionalizer.
     */
    static
    void nondimensionalize(Mesh* const mesh,
                           const spatialdata::units::Nondimensional& normalizer);

    /** Check topology of mesh.
     *
     * @param[in] mesh Finite-element mesh.
     */
    static
    void checkTopology(const Mesh& mesh);

    /** Determine is mesh contains simplex cells (i.e., line, tri, tet).
     *
     * @param[in] mesh Finite-element mesh.
     * @returns True if mesh contains simplex cells.
     */
    static
    bool isSimplexMesh(const Mesh& mesh);

    /** Determine if cell is a cohesive cell.
     *
     * @param[in] dm PETSc DM containining cell.
     * @param[in] cell Point number.
     * @returns True if cell is a cohesive cell.
     */
    static
    bool isCohesiveCell(const PetscDM dm,
                        const PetscInt cell);

    /** Get number of vertices in mesh.
     *
     * @param[in] mesh Finite-element mesh.
     * @returns Number of vertices in mesh.
     */
    static
    PylithInt getNumVertices(const pylith::topology::Mesh& mesh);

    /** Get number of cells in mesh.
     *
     * @param[in] mesh Finite-element mesh.
     * @returns Number of cells in mesh.
     */
    static
    PylithInt getNumCells(const pylith::topology::Mesh& mesh);

    /** Get number of vertices in a cell in mesh.
     *
     * @param[in] mesh Finite-element mesh.
     * @returns Number of vertices in a cell.
     */
    static
    PylithInt getNumCorners(const pylith::topology::Mesh& mesh);

    /** Check to make sure material id of every cell matches the id of
     *  one of the materials.
     *
     * @param[in] mesh Finite-element mesh.
     * @param[in] materialIds Array of ids for all materials and interior interfaces.
     */
    static
    void checkMaterialIds(const Mesh& mesh,
                          pylith::int_array& materialIds);

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    MeshOps(void); ///< Not Implemented
    MeshOps(const MeshOps&); ///< Not implemented
    const MeshOps& operator=(const MeshOps&); ///< Not implemented

}; // MeshOps

#endif // pylith_topology_meshops_hh

// End of file
