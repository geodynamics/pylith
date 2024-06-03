// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2024, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================
#pragma once

#include "pylith/topology/topologyfwd.hh" // forward declarations

#include "pylith/utils/petscfwd.h" // HASA PetscDM, PetscIS

// Stratum --------------------------------------------------------
/// Height or depth stratum.
class pylith::topology::Stratum { // Stratum
    friend class TestStratum; // unit testing

    // PUBLIC ENUMS /////////////////////////////////////////////////////////
public:

    /// Type of stratum (height or depth).
    enum StratumEnum {
        HEIGHT=0,
        DEPTH=1,
    }; // StratumEnum

    // PUBLIC METHODS ///////////////////////////////////////////////////////
public:

    /** Default constructor.
     *
     * @param dmMesh PETSc DM for finite-element mesh.
     * @param stype Type of stratum [HEIGHT, DEPTH].
     * @param level Height of depth of stratum.
     */
    Stratum(const PetscDM dmMesh,
            const StratumEnum stype,
            const int level);

    /// Default destructor.
    ~Stratum(void);

    /** Get starting point.
     *
     * @return Index of starting point.
     */
    PetscInt begin(void) const;

    /** Get ending point.
     *
     * @return Index of ending point.
     */
    PetscInt end(void) const;

    /** Get number of points in stratum.
     *
     * @return Number of points.
     */
    PetscInt size(void) const;

    // PRIVATE MEMBERS //////////////////////////////////////////////////////
private:

    PetscInt _begin; ///< Starting point.
    PetscInt _end; ///< End point.

}; // Stratum

// StratumIS ------------------------------------------------------------
/// Index set associated with stratum (usually over label of points).
class pylith::topology::StratumIS { // StratumIS
    friend class TestStratumIS; // unit testing

    // PUBLIC METHODS ///////////////////////////////////////////////////////
public:

    /** Default constructor.
     *
     * @param dmMesh PETSc DM for finite-element mesh.
     * @param label Label for stratum.
     * @param id Value of label defining stratum.
     */
    StratumIS(const PetscDM dmMesh,
              const char* label,
              const int id,
              const bool includeOnlyCells=false);

    /// Default destructor.
    ~StratumIS(void);

    /// Deallocate data.
    void deallocate(void);

    /** Get array of points.
     *
     * @return Array of points.
     */
    const PetscInt* points(void) const;

    /** Get number of points in index set.
     *
     * @return Number of points.
     */
    PetscInt size(void) const;

    // PRIVATE MEMBERS //////////////////////////////////////////////////////
private:

    PetscIS _indexSet; ///< PETSc index set.
    PetscInt _size; ///< Size of index set.
    const PetscInt* _points; ///< Array of points in index set.

}; // StratumIS

#include "Stratum.icc"

// End of file
