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

#include "spatialdata/spatialdb/spatialdbfwd.hh" // USES SpatialDB
#include <map>

// Wrapper for PetscFE
class pylith::topology::FE {
public:

    FE(PetscFE fe) : _fe(fe) {
        PetscObjectReference((PetscObject) _fe);
    }

    FE(const FE& fe) : _fe(fe._fe) {
        PetscObjectReference((PetscObject) _fe);
    }

    ~FE() {
        PetscInt refct = -1;

        if (_fe) {PetscObjectGetReference((PetscObject) _fe, &refct);}
        PetscObjectDereference((PetscObject) _fe);
    }

    PetscFE _fe;
};

// FieldOps -------------------------------------------------------------
class pylith::topology::FieldOps {
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
                     const int numComponents);

    /** Check compatibility of discretization of subfields in the auxiliary field and target field.
     *
     * @param[in] target Field with subfields set from auxiliary field.
     * @param[in] auxliary Auxiliary field.
     */
    static
    void checkDiscretization(const pylith::topology::Field& target,
                             const pylith::topology::Field& auxiliary);

    /** Check that 'field' contains required fields for 'reason'.
     *
     * @throws std::runtime_error if field does not contain all required fields.
     *
     * @param[in] requiredFields Names of required fields.
     * @param[in] reason Reason for required fields.
     * @param[in] field Field needing required fields.
     */
    static
    void checkSubfieldsExist(const pylith::string_vector& requiredFields,
                             const std::string& reason,
                             const pylith::topology::Field& field);

    /** Get names of subfields extending over the entire domain.
     *
     * Excludes subfields extending over only a subset of the domain, like the fault Lagrange multiplier subfield.
     *
     * @param[in] field Field with subfields.
     * @returns Array with names of subfields.
     */
    static
    pylith::string_vector getSubfieldNamesDomain(const pylith::topology::Field& field);

    /** Check to see if fields have the same subfields and match in size.
     *
     * @param[in] fieldA Field to check.
     * @param[in] fieldB Field to check.
     * @returns true if fields have the same subfields and match in size; otherwise false.
     */
    static
    bool layoutsMatch(const pylith::topology::Field& fieldA,
                      const pylith::topology::Field& fieldB);

    /** Free saved PetscFE objects.
     */
    static
    void deallocate(void);

    // NOT IMPLEMENTED //////////////////////////////////////////////////////
private:

    FieldOps(void); ///< Not implemented.
    FieldOps(const FieldOps&); ///< Not implemented.
    const FieldOps& operator=(const FieldOps&); ///< Not implemented.

    static std::map<FieldBase::Discretization, pylith::topology::FE> feStore;

}; // FieldOps

#endif // pylith_topology_fieldOps_hh

// End of file
