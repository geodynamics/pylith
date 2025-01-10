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

#include "pylith/faults/faultsfwd.hh" // forward declarations

#include "pylith/utils/PyreComponent.hh" // ISA PyreComponent

#include "pylith/topology/FieldBase.hh" // USES FieldBase::SpaceDim
#include "pylith/feassemble/feassemblefwd.hh" // USES AuxiliaryFactory
#include "pylith/utils/petscfwd.h" // HASA PetscVec

#include "spatialdata/geocoords/geocoordsfwd.hh" // USES CoordSys
#include "spatialdata/spatialdb/spatialdbfwd.hh" // USES SpatialDB
#include "spatialdata/units/unitsfwd.hh" // USES Nondimensional

// KinSrc -------------------------------------------------------------
/** @brief Kinematic earthquake source.
 *
 * KinSrc is responsible for providing the value of slip, slip rate, or
 * slip acceleration at time t over a fault surface.
 *
 * The auxiliary subfields in this object use the discretization from the
 * fault 'slip' auxiliary subfield.
 */
class pylith::faults::KinSrc : public pylith::utils::PyreComponent {
    friend class TestKinSrc; // unit testing

    // PUBLIC MEMBERS /////////////////////////////////////////////////////////////////////////////
public:

    static const int GET_SLIP;
    static const int GET_SLIP_RATE;
    static const int GET_SLIP_ACC;

    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    /// Default constructor.
    KinSrc(void);

    /// Destructor.
    virtual ~KinSrc(void);

    /// Deallocate PETSc and local data structures.
    virtual
    void deallocate(void);

    /** Set origin time for earthquake source.
     *
     * @param value Origin time for earthquake source.
     */
    void setOriginTime(const PylithReal value);

    /** Get origin time for earthquake source.
     *
     * @returns Origin time for earthquake source.
     */
    PylithReal getOriginTime(void) const;

    /** Get auxiliary field associated with the kinematic source.
     *
     * @return field Auxiliary field for the kinematic source.
     */
    const pylith::topology::Field& auxField(void) const;

    /** Set the spatial database for filling auxiliary subfields.
     *
     * @param[in] value Pointer to database.
     */
    void auxFieldDB(spatialdata::spatialdb::SpatialDB* value);

    /** Initialize kinematic (prescribed slip) earthquake source.
     *
     * @param[in] auxField Auxiliary field associated with fault finite-element integration.
     * @param[in] normalizer Normalizer for nondimensionalizing values.
     * @param[in] cs Coordinate system for problem.
     */
    void initialize(const pylith::topology::Field& auxField,
                    const spatialdata::units::Nondimensional& normalizer,
                    const spatialdata::geocoords::CoordSys* cs);

    /** Get requested slip subfields at time t.
     *
     * @param[inout] slipLocalVec Local PETSc vector for slip, slip rate, or slip accelerationvalues.
     * @param[in] faultAuxiliaryField Auxiliary field for fault.
     * @param[in] t Time t.
     * @param[in] timeScale Time scale for nondimensionalization.
     * @param[in] bitSlipSubfields Slip subfields to compute.
     */
    virtual
    void getSlipSubfields(PetscVec slipLocalVec,
                          pylith::topology::Field* faultAuxiliaryField,
                          const PylithScalar t,
                          const PylithScalar timeScale,
                          const int bitSlipSubfields);

    // PROTECTED METHODS //////////////////////////////////////////////////////////////////////////
protected:

    /** Setup auxiliary subfields (discretization and query fns).
     *
     * Create subfields in auxiliary fields (includes name of the field,
     * vector field type, discretization, and scale for
     * nondimensionalization) and set query functions for filling them
     * from a spatial database.
     *
     * @attention The order of the calls to subfieldAdd() must match the
     * order of the auxiliary fields in the FE kernels.
     *
     * @param[in] normalizer Normalizer for nondimensionalizing values.
     * @param[in] cs Coordinate system for problem.
     */
    virtual
    void _auxiliaryFieldSetup(const spatialdata::units::Nondimensional& normalizer,
                              const spatialdata::geocoords::CoordSys* cs) = 0;

    /** Set constants used in finite-element integrations.
     *
     * @param[in] auxField Auxiliary field associated with fault finite-element integration.
     */
    void _setFEConstants(const pylith::topology::Field& auxField) const;

    // PROTECTED MEMBERS //////////////////////////////////////////////////////////////////////////
protected:

    pylith::faults::KinSrcAuxiliaryFactory* _auxiliaryFactory; ///< Factory for auxiliary subfields.
    PetscPointFunc _slipFnKernel; ///< Kernel for slip time function.
    PetscPointFunc _slipRateFnKernel; ///< Kernel for slip rate time function.
    PetscPointFunc _slipAccFnKernel; ///< Kernel for slip acceleration time function.
    pylith::topology::Field* _auxiliaryField; ///< Auxiliary field for this integrator.

    // PRIVATE MEMBERS ////////////////////////////////////////////////////////////////////////////
private:

    PylithReal _originTime; ///< Origin time for earthquake source

    // NOT IMPLEMENTED ////////////////////////////////////////////////////////////////////////////
private:

    KinSrc(const KinSrc&); ///< Not implemented
    const KinSrc& operator=(const KinSrc&); ///< Not implemented

}; // class KinSrc

// End of file
