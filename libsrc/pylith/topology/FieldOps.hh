// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

/**
 * @file libsrc/topology/FieldOps.hh
 *
 * @brief Operations on fields.
 */

#if !defined(pylith_topology_fieldops_hh)
#define pylith_topology_fieldops_hh

// Include directives ---------------------------------------------------
#include "pylith/topology/topologyfwd.hh" // forward declarations

#include "FieldBase.hh" // USES FieldBase::Discretization
#include "pylith/utils/petscfwd.h" // USES PetscFE

// FieldOps -------------------------------------------------------------
/// @brief C++ class for simple operations for a Field object.
class pylith::topology::FieldOps { // FieldOps
    friend class TestFieldOps; // unit testing

    // PUBLIC METHODS ///////////////////////////////////////////////////////
public:

    /** Create PetscFE object for discretization of field.
     *
     * @param[in] feinfo Discretization information for field.
     * @param[in] dm PetscDM for finite-element mesh.
     * @param[in] isSimplex True if mesh contains simplex cells.
     * @param[in] numComponents Number of components in field.
     *
     * @returns PetscFE object.
     *
     * Discretization parameters collected from:
     *   + User
     *     - order of basis
     *     - continuity of basis
     *     - quadrature order
     *   + Mesh
     *     - spatial dimension
     *     - type of cell: simplex (tri/tet) or quad/hex
     *   + Field
     *     - number of components
     */
    static
    PetscFE createFE(const FieldBase::Discretization& feinfo,
                     const PetscDM dm,
                     const bool isSimplex,
                     const int numComponents);


    // NOT IMPLEMENTED //////////////////////////////////////////////////////
private:

    FieldOps(void); ///< Not implemented.
    FieldOps(const FieldOps&); ///< Not implemented.
    const FieldOps& operator=(const FieldOps&); ///< Not implemented.

}; // FieldOps

#endif // pylith_topology_fieldOps_hh


// End of file
