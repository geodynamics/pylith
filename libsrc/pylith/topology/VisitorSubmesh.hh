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
 * @file libsrc/topology/VisitorSubmesh.hh
 *
 * @brief C++ helper class for accessing field and matrix values at
 * points in a submesh within a finite-element mesh.
 *
 * This visitor is used to access values associated with a submesh
 * when the field or matrix is defined over the entire mesh. This is
 * why a submesh and index set are passed to the constructor.
 *
 * Use the Vec/MatVisitorMesh objects when the field and mesh/submesh
 * are defined over the same set of points (i.e., field over a submesh
 * or field of a mesh).
 */

#if !defined(pylith_topology_visitorsubmesh_hh)
#define pylith_topology_visitorsubmesh_hh

// Include directives ---------------------------------------------------
#include "topologyfwd.hh" // forward declarations

#include "pylith/utils/petscfwd.h" // HASA PetscVec, PetscSection

// VecVisitorSubmesh -------------------------------------------------------
/** @brief Helper class for accessing field values at points in a
 *  finite-element mesh.
 */
class pylith::topology::VecVisitorSubmesh { // VecVisitorSubmesh
    friend class TestVecVisitorSubmesh; // unit testing

    // PUBLIC METHODS ///////////////////////////////////////////////////////
public:

    /** Default constructor.
     *
     * @param[in] field Field associated with visitor.
     * @param[in] submeshIS Submesh index set associated with visitor.
     */
    VecVisitorSubmesh(const Field& field,
                      const SubmeshIS& submeshIS);

    /// Default destructor
    ~VecVisitorSubmesh(void);

    /* Initialize cached data.
     *
     * @param[in] submeshIS Submesh index set associated with visitor.
     */
    void initialize(const SubmeshIS& submeshIS);

    /// Clear cached data.
    void clear(void);

    /** Get the PETSc section.
     *
     * @returns PETSc section.
     */
    PetscSection petscSection(void) const;

    /** Get the local PETSc Vec.
     *
     * @returns PETSc Vec.
     */
    PetscVec localVec(void) const;

    /** Get number of degrees of freedom for field at point.
     *
     * @param[in] point Point in mesh.
     * @returns Number of degrees of freedom.
     */
    PetscInt sectionDof(const PetscInt point) const;

    /** Get offset into array for point.
     *
     * @param[in] point Point in mesh.
     * @returns Offset.
     */
    PetscInt sectionOffset(const PetscInt point) const;

    /** Get number of degrees of freedom for subfield at point.
     *
     * @param[in] point Point in mesh.
     * @param[in] subfieldIndex Index of subfield in field.
     * @returns Number of degrees of freedom.
     */
    PetscInt sectionSubfieldDof(const PetscInt point,
                                const PetscInt subfieldIndex) const;

    /** Get offset into array for point.
     *
     * @param[in] point Point in mesh.
     * @param[in] subfieldIndex Index of subfield in field.
     * @returns Offset.
     */
    PetscInt sectionSubfieldOffset(const PetscInt point,
                                   const PetscInt subfieldIndex) const;

    /** Get array of values associated with closure.
     *
     * @param[in] valuesCell Array of values for cell.
     * @param[in] valuesSize Size of values array.
     * @param[in] cell Finite-element cell.
     */
    void getClosure(PetscScalar** valuesCell,
                    PetscInt* valuesSize,
                    const PetscInt cell) const;

    /** Get array of values associated with closure.
     *
     * @param[in] values Array of values for cell.
     * @param[in] cell Finite-element cell.
     */
    void getClosure(scalar_array* values,
                    const PetscInt cell) const;

    /** Restore array of values associated with closure.
     *
     * @param[in] valuesCell Array of values for cell.
     * @param[in] valuesSize Size of values array.
     * @param[in] cell Finite-element cell.
     */
    void restoreClosure(PetscScalar** valuesCell,
                        PetscInt* valuesSize,
                        const PetscInt cell) const;

    /** Set values associated with closure.
     *
     * @param[in] valuesCell Array of values for cell.
     * @param[in] valuesSize Size of values array.
     * @param[in] cell Finite-element cell.
     * @param[in] mode Mode for inserting values.
     */
    void setClosure(const PetscScalar* valuesCell,
                    const PetscInt valuesSize,
                    const PetscInt cell,
                    const InsertMode mode) const;

    // PRIVATE MEMBERS //////////////////////////////////////////////////////
private:

    const Field& _field;

    PetscDM _dm; ///< Cached PETSc dm for submesh.
    PetscSection _section; ///< Cached PETSc subsection.
    PetscVec _localVec; ///< Cached local PETSc Vec.

    // NOT IMPLEMENTED //////////////////////////////////////////////////////
private:

    VecVisitorSubmesh(const VecVisitorSubmesh&); ///< Not implemented
    const VecVisitorSubmesh& operator=(const VecVisitorSubmesh&); ///< Not implemented

}; // VecVisitorSubmesh

// MatVisitorSubmesh -------------------------------------------------------
/** @brief Helper class for accessing field values at points in a
 *  finite-element mesh.
 */
class pylith::topology::MatVisitorSubmesh { // MatVisitorSubmesh
    friend class TestMatVisitorSubmesh; // unit testing

    // PUBLIC METHODS ///////////////////////////////////////////////////////
public:

    /** Default constructor.
     *
     * @param[in] mat PETSc matrix.
     * @param[in] field Field associated with visitor.
     * @param[in] submeshIS Submesh index set associated with visitor.
     */
    MatVisitorSubmesh(const PetscMat mat,
                      const Field& field,
                      const SubmeshIS& submeshIS);

    /// Default destructor
    ~MatVisitorSubmesh(void);

    // Initialize.
    void initialize(void);

    /// Clear cached data.
    void clear(void);

    /** Get array of values associated with closure.
     *
     * @param[in] valuesCell Array of values for cell.
     * @param[in] valuesSize Size of values array.
     * @param[in] cell Finite-element cell.
     */
    void getClosure(PetscScalar** valuesCell,
                    PetscInt* valuesSize,
                    const PetscInt cell) const;

    /** Restore array of values associated with closure.
     *
     * @param[in] valuesCell Array of values for cell.
     * @param[in] valuesSize Size of values array.
     * @param[in] cell Finite-element cell.
     */
    void restoreClosure(PetscScalar** valuesCell,
                        PetscInt* valuesSize,
                        const PetscInt cell) const;

    /** Set values associated with closure.
     *
     * @param[in] valuesCell Array of values for cell.
     * @param[in] valuesSize Size of values array.
     * @param[in] cell Finite-element cell.
     * @param[in] mode Mode for inserting values.
     */
    void setClosure(const PetscScalar* valuesCell,
                    const PetscInt valuesSize,
                    const PetscInt cell,
                    const InsertMode mode) const;

    // PRIVATE MEMBERS //////////////////////////////////////////////////////
private:

    const PetscMat _mat; ///< Cached PETSc matrix.
    PetscDM _dm; ///< Cached PETSc dm for mesh.
    PetscSection _subsection; ///< Cached PETSc section for submesh.
    PetscSection _globalSection; ///< Cached PETSc global section.
    PetscSection _globalSubsection; ///< Cached PETSc subsection.

    // NOT IMPLEMENTED //////////////////////////////////////////////////////
private:

    MatVisitorSubmesh(const MatVisitorSubmesh&); ///< Not implemented
    const MatVisitorSubmesh& operator=(const MatVisitorSubmesh&); ///< Not implemented

}; // MatVisitorSubmesh

// SubmeshIS ------------------------------------------------------------
/// Index set associated with submesh.
class pylith::topology::SubmeshIS { // SubmeshIS
    friend class TestSubmeshIS; // unit testing

    // PUBLIC METHODS ///////////////////////////////////////////////////////
public:

    /** Default constructor.
     *
     * @param[in] submesh Submesh associated with index set.
     */
    SubmeshIS(const Mesh& submesh);

    /// Default destructor.
    ~SubmeshIS(void);

    /// Deallocate.
    void deallocate(void);

    /** Get the submesh.
     *
     * @returns Submesh.
     */
    const Mesh& submesh(void) const;

    /** Get PETSc index set.
     *
     * @returns PETSc index set.
     */
    PetscIS indexSet(void) const;

    /** Get array of points in index set.
     *
     * @return Array of points.
     */
    const PetscInt* indices(void) const;

    /** Get number of points in index set.
     *
     * @return Number of points.
     */
    PetscInt size(void) const;

    // PRIVATE MEMBERS //////////////////////////////////////////////////////
private:

    const Mesh& _submesh;
    PetscIS _indexSet; ///< PETSc index set.
    PetscInt _size; ///< Size of index set.
    const PetscInt* _indices; ///< Array of indices of points in index set.

}; // SubmeshIS

#include "VisitorSubmesh.icc"

#endif // pylith_topology_visitorsubmesh_hh

// End of file
