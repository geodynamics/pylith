// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2023, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================
#pragma once

#include "pylith/topology/topologyfwd.hh" // forward declarations

#include "pylith/topology/FieldBase.hh" // USES FieldBase::Discretization
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

    /** Create label for output.
     *
     * @param[in] field Field to which output label is added.
     */
    static
    void createOutputLabel(const pylith::topology::Field* field);

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

// End of file
