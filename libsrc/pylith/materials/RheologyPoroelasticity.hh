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

#include "pylith/materials/materialsfwd.hh" // forward declarations
#include "pylith/utils/PyreComponent.hh" // ISA PyreComponent

#include "pylith/topology/topologyfwd.hh" // USES Field
#include "pylith/utils/arrayfwd.hh" // USES std::vector
#include "pylith/feassemble/IntegratorDomain.hh" // USES IntegratorDomain::ProjectKernels

#include "spatialdata/geocoords/geocoordsfwd.hh" // USES Coordsys

#include "petscds.h" // USES PetscPointFn*, PetscPointJacFn*

class pylith::materials::RheologyPoroelasticity : public pylith::utils::PyreComponent {
    friend class TestIsotropicLinearPoroelasticity; // unit testing

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Default constructor.
    RheologyPoroelasticity(void);

    /// Destructor.
    virtual ~RheologyPoroelasticity(void);

    /// Deallocate PETSc and local data structures.
    void deallocate(void);

    /** Get auxiliary factory associated with physics.
     *
     * @return Auxiliary factory for physics object.
     */
    virtual
    pylith::materials::AuxiliaryFactoryPoroelastic* getAuxiliaryFactory(void) = 0;

    /// Add rheology subfields to auxiliary field.
    virtual
    void addAuxiliarySubfields(void) = 0;

    // ============================= RHS ==================================== //

    // ---------------------------------------------------------------------------------------------------------------------
    // Select g0p function. Will only be used for the dynamic case.
    virtual
    PetscPointFn* getKernelg0p(const spatialdata::geocoords::CoordSys* coordsys,
                               const bool _useBodyForce,
                               const bool _gravityField,
                               const bool _useSourceDensity) const = 0;

    // ---------------------------------------------------------------------------------------------------------------------
    /** Get pressure kernel for RHS residual, G(t,s).
     *
     * @param[in] coordsys Coordinate system.
     *
     * @return RHS residual kernel for Darcy velocity.
     */
    virtual
    PetscPointFn* getKernelg1p_explicit(const spatialdata::geocoords::CoordSys* coordsys,
                                        const bool _gravityField) const = 0;

    // ---------------------------------------------------------------------------------------------------------------------
    // Get stress kernel for RHS residual, G(t,s)
    virtual
    PetscPointFn* getKernelg1v_explicit(const spatialdata::geocoords::CoordSys* coordsys) const = 0;

    // =============================== LHS =================================== //

    // ---------------------------------------------------------------------------------------------------------------------
    // Get variation in fluid content kernel for LHS residual, F(t,s,\dot{s})
    virtual
    PetscPointFn* getKernelf0p_explicit(const spatialdata::geocoords::CoordSys* coordsys) const = 0;

    // ---------------------------------------------------------------------------------------------------------------------
    // Select implicit f0p function.
    virtual
    PetscPointFn* getKernelf0p_implicit(const spatialdata::geocoords::CoordSys* coordsys,
                                        const bool _useBodyForce,
                                        const bool _gravityField,
                                        const bool _useSourceDensity) const = 0;

    // ---------------------------------------------------------------------------------------------------------------------
    // Get stress kernel for LHS residual, F(t,s,\dot{s})
    virtual
    PetscPointFn* getKernelf1u_implicit(const spatialdata::geocoords::CoordSys* coordsys) const = 0;

    // ---------------------------------------------------------------------------------------------------------------------
    /** Get pressure kernel for LHS residual.
     *
     * @param[in] coordsys Coordinate system.
     *
     * @return LHS residual kernel for Darcy velocity.
     */
    virtual
    PetscPointFn* getKernelf1p_implicit(const spatialdata::geocoords::CoordSys* coordsys,
                                        const bool _useBodyForce,
                                        const bool _gravityField) const = 0;

    // ---------------------------------------------------------------------------------------------------------------------
    // Get poroelastic constants kernel for LHS Jacobian
    virtual
    PetscPointJacFn* getKernelJf3uu(const spatialdata::geocoords::CoordSys* coordsys) const = 0;

    // ---------------------------------------------------------------------------------------------------------------------
    // Get biot coefficient kernel for LHS Jacobian
    virtual
    PetscPointJacFn* getKernelJf2up(const spatialdata::geocoords::CoordSys* coordsys) const = 0;

    // ---------------------------------------------------------------------------------------------------------------------
    // Get lambda kernel for LHS Jacobian
    virtual
    PetscPointJacFn* getKernelJf2ue(const spatialdata::geocoords::CoordSys* coordsys) const = 0;

    // ---------------------------------------------------------------------------------------------------------------------
    // Get Specific storage kernel for LHS Jacobian F(t,s, \dot{s}).
    virtual
    PetscPointJacFn* getKernelJf0pp(const spatialdata::geocoords::CoordSys* coordsys) const = 0;

    // ---------------------------------------------------------------------------------------------------------------------
    // Get Darcy Conductivity kernel for LHS Jacobian
    virtual
    PetscPointJacFn* getKernelJf3pp(const spatialdata::geocoords::CoordSys* coordsys) const = 0;

    // ---------------------------------------------------------------------------------------------------------------------
    // Get biot coefficient kernel for LHS Jacobian F(t,s, \dot{s}).
    virtual
    PetscPointJacFn* getKernelJf0pe(const spatialdata::geocoords::CoordSys* coordsys) const = 0;

    // ============================ DERIVED FIELDS ========================== //

    /** Get stress kernel for derived field.
     *
     * @param[in] coordsys Coordinate system.
     *
     * @return Project kernel for computing stress subfield in derived field.
     */
    virtual
    PetscPointFn* getKernelCauchyStressVector(const spatialdata::geocoords::CoordSys* coordsys) const = 0;

    /** Get water content kernel for derived field.
     *
     * @param[in] coordsys Coordinate system.
     *
     * @return Project kernel for computing stress subfield in derived field.
     */
    virtual
    PetscPointFn* getKernelWaterContent(const spatialdata::geocoords::CoordSys* coordsys) const = 0;

    /** Update kernel constants.
     *
     * @param[inout] kernelConstants Array of constants used in integration kernels.
     * @param[in] dt Current time step.
     */
    virtual
    void updateKernelConstants(pylith::real_array* kernelConstants,
                               const PylithReal dt) const;

    /** Add kernels for updating state variables, implicit.
     *
     * @param[inout] kernels Array of kernels for updating state variables.
     * @param[in] coordsys Coordinate system.
     */
    virtual
    void addKernelsUpdateStateVarsImplicit(std::vector<pylith::feassemble::IntegratorDomain::ProjectKernels>* kernels,
                                           const spatialdata::geocoords::CoordSys* coordsys,
                                           const bool _useStateVars) const;

    /** Add kernels for updating state variables, explicit.
     *
     * @param[inout] kernels Array of kernels for updating state variables.
     * @param[in] coordsys Coordinate system.
     */
    virtual
    void addKernelsUpdateStateVarsExplicit(std::vector<pylith::feassemble::IntegratorDomain::ProjectKernels>* kernels,
                                           const spatialdata::geocoords::CoordSys* coordsys,
                                           const bool _useStateVars) const;

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    RheologyPoroelasticity(const RheologyPoroelasticity&); ///< Not implemented.
    const RheologyPoroelasticity& operator=(const RheologyPoroelasticity&); /// Not implemented.

}; // class RheologyPoroelasticity

// End of file
