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
 * @file libsrc/topology/Mesh.hh
 *
 * @brief C++ PyLith finite-element mesh.
 */

#if !defined(pylith_topology_mesh_hh)
#define pylith_topology_mesh_hh

// Include directives ---------------------------------------------------
#include "topologyfwd.hh" // forward declarations
#include "spatialdata/geocoords/geocoordsfwd.hh" // forward declarations

#include "pylith/utils/petscfwd.h" // HASA PetscDM

// Mesh -----------------------------------------------------------------
/** @brief PyLith finite-element mesh.
 *
 * Extends PETSc mesh to include coordinate system associated with
 * domain.
 */
class pylith::topology::Mesh { // Mesh
    friend class TestMesh; // unit testing

    // PUBLIC METHODS ///////////////////////////////////////////////////////
public:

    /// Default constructor.
    Mesh(void);

    /** Constructor with dimension and communicator.
     *
     * @param dim Dimension associated with mesh cells.
     * @param comm MPI communicator for mesh.
     */
    Mesh(const int dim,
         const MPI_Comm& comm=PETSC_COMM_WORLD);

    /// Default destructor
    ~Mesh(void);

    /// Deallocate PETSc and local data structures.
    void deallocate(void);

    /** Create clone.
     *
     * @returns Clone of mesh.
     */
    Mesh* clone(void) const;

    /** Get DMPlex mesh.
     *
     * @returns DMPlex mesh.
     */
    PetscDM getDM(void) const;

    /** Set DMPlex mesh.
     *
     * @param DMPlex mesh.
     * @param label Label for mesh.
     */
    void setDM(PetscDM dm,
               const char* label="domain");

    /** Set coordinate system.
     *
     * @param cs Coordinate system.
     */
    void setCoordSys(const spatialdata::geocoords::CoordSys* cs);

    /** Get coordinate system.
     *
     * @returns Coordinate system.
     */
    const spatialdata::geocoords::CoordSys* getCoordSys(void) const;

    /** Get dimension of mesh.
     *
     * @returns Dimension of mesh.
     */
    int getDimension(void) const;

    /** Get MPI communicator associated with mesh.
     *
     * @returns MPI communicator.
     */
    MPI_Comm getComm(void) const;

    /** Get MPI rank.
     *
     * @returns MPI rank.
     */
    int getCommRank(void) const;

    /** View mesh.
     *
     * @param viewOption PETSc DM view option.
     *
     * PETSc mesh view options include:
     *   short summary [empty]
     *   detail summary ::ascii_info_detail
     *   detail in a file :refined.mesh:ascii_info_detail
     *   latex in a file  :refined.tex:ascii_latex
     *   VTK vtk:refined.vtk:ascii_vtk
     */
    void view(const char* viewOption="::ascii_info_detail") const;

    /** Get name of label for all mesh cells, including hybrid cells.
     *
     * @returns Name of label.
     */
    static
    const char* const getCellsLabelName(void);

    // PRIVATE MEMBERS //////////////////////////////////////////////////////
private:

    spatialdata::geocoords::CoordSys* _coordSys; ///< Coordinate system.
    PetscDM _dm; ///< PETSc DM with topology.

    // NOT IMPLEMENTED //////////////////////////////////////////////////////
private:

    Mesh(const Mesh&); ///< Not implemented
    const Mesh& operator=(const Mesh&); ///< Not implemented

}; // Mesh

#endif // pylith_topology_mesh_hh

// End of file
