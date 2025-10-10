// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================
#pragma once

#include "pylith/topology/topologyfwd.hh" // forward declarations

#include "pylith/utils/petscfwd.h" // USES PetscDM
#include "pylith/utils/array.hh" // USES int_array

#include "spatialdata/geocoords/geocoordsfwd.hh"
#include "pylith/scales/scalesfwd.hh"

class pylith::topology::MeshOps {
    friend class TestMeshOps; // unit testing

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /** Create subdomain mesh using label.
     *
     * @param[in] mesh Mesh for domain.
     * @param[in] labelName Name of label marking subdomain.
     * @param[in] labelValue Value of label marking subdomain.
     * @param[in] componentName Name of component associated with subdomain.
     *
     * @returns Mesh for subdomain.
     */
    static
    pylith::topology::Mesh* createSubdomainMesh(const pylith::topology::Mesh& mesh,
                                                const char* labelName,
                                                const int labelValue,
                                                const char* componentName);

    /** Create lower dimension mesh using label.
     *
     * @param[in] mesh Mesh for domain.
     * @param[in] labelName Name of label marking subdomain.
     * @param[in] labelValue Value of label marking subdomain.
     * @param[in] componentName Name of component associated with subdomain.
     *
     * @returns Lower dimension mesh.
     */
    static
    pylith::topology::Mesh* createLowerDimMesh(const pylith::topology::Mesh& mesh,
                                               const char* labelName,
                                               const int labelValue,
                                               const char* componentName);

    /** Create 0-dimension mesh from points.
     *
     * @param[in] pointCoords Array of coordinates of points [numPoints*spaceDim].
     * @param[in] numPoints Number of points.
     * @param[in] cs Coordinate system for points.
     * @param[in] lengthScale Length scale for nondimensionalization.
     * @param[in] comm MPI communicator.
     * @param[in] componentName Name of component associated with subdomain.
     */
    static
    pylith::topology::Mesh* createFromPoints(const PylithReal* points,
                                             const size_t numPoints,
                                             const spatialdata::geocoords::CoordSys* cs,
                                             const PylithReal lengthScale,
                                             MPI_Comm comm,
                                             const char* componentName);

    /** Remove cells hanging off mesh.
     *
     * @param[in] dmMesh PETSc DM to clean.
     *
     * @returns PETSc DM without hanging cells.
     */
    static
    PetscDM removeHangingCells(const PetscDM& dmMesh);

    /** Nondimensionalize the finite-element mesh.
     *
     * @param[in] mesh Finite-element mesh.
     * @param[in] scales Nondimensionalizer.
     */
    static
    void nondimensionalize(Mesh* const mesh,
                           const pylith::scales::Scales& scales);

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

    /** Compute nominal dimension of domain based on mesh bounding box.
     *
     * @param[in] mesh Finite-element mesh.
     * @returns Average domain dimension.
     */
    static
    PylithReal computeAvgDomainDim(const pylith::topology::Mesh& mesh);

    /** Check to make sure material label value for every cell matches the label value of
     *  one of the materials.
     *
     * @param[in] mesh Finite-element mesh.
     * @param[in] labelValues Array of label values for all materials and interior interfaces.
     */
    static
    void checkMaterialLabels(const Mesh& mesh,
                             pylith::int_array& labelValues);

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    MeshOps(void); ///< Not Implemented
    MeshOps(const MeshOps&); ///< Not implemented
    const MeshOps& operator=(const MeshOps&); ///< Not implemented

}; // MeshOps

// End of file
