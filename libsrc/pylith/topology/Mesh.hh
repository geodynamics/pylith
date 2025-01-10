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
#include "spatialdata/geocoords/geocoordsfwd.hh" // forward declarations

#include "pylith/utils/petscfwd.h" // HASA PetscDM

// Mesh -----------------------------------------------------------------
/** @brief PyLith finite-element mesh.
 *
 * Extends PETSc mesh to include coordinate system associated with
 * domain.
 */
class pylith::topology::Mesh {
    friend class TestMesh; // unit testing

    // PUBLIC MEMBERS ///////////////////////////////////////////////////////
public:

    static const char* cells_label_name; ///< Name of label covering all cells.

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
     *   python in a file :mesh.py:ascii_python
     *   VTU vtk:refined.vtu:ascii_vtk_vtu
     *   HDF5 PETSc hdf5:mesh.h5:hdf5_petsc
     *   HDF5 XDMF hdf5:mesh.h5:hdf5_xdmf
     */
    void view(const char* viewOption="::ascii_info_detail") const;

    // PRIVATE MEMBERS //////////////////////////////////////////////////////
private:

    spatialdata::geocoords::CoordSys* _coordSys; ///< Coordinate system.
    PetscDM _dm; ///< PETSc DM with topology.

    // NOT IMPLEMENTED //////////////////////////////////////////////////////
private:

    Mesh(const Mesh&); ///< Not implemented
    const Mesh& operator=(const Mesh&); ///< Not implemented

}; // Mesh

// End of file
