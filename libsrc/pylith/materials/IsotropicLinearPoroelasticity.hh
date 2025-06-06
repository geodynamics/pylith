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

#include "pylith/materials/RheologyPoroelasticity.hh" // ISA RheologyPoroelasticity

class pylith::materials::IsotropicLinearPoroelasticity : public pylith::materials::RheologyPoroelasticity {
    friend class TestIsotropicLinearPoroelasticity; // unit testing

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Default constructor.
    IsotropicLinearPoroelasticity(void);

    /// Destructor.
    ~IsotropicLinearPoroelasticity(void);

    /// Deallocate PETSc and local data structures.
    void deallocate(void);

    /** Include reference stress/strain?
     *
     * @param value Flag indicating to include reference stress/strain.
     */
    void useReferenceState(const bool value);

    /** Use reference stress and strain in computation of stress and
     * strain?
     *
     * @returns True if using reference stress and strain, false otherwise.
     */
    bool useReferenceState(void) const;

    /** Include tensor permeability?
     *
     * @param value Flag indicating to include tensor permeability.
     */
    void useTensorPermeability(const bool value);

    /** Use full tensor permeability?
     *
     * @returns True if using full tensor permeability, false otherwise.
     */
    bool useTensorPermeability(void) const;

    /** Get auxiliary factory associated with physics.
     *
     * @return Auxiliary factory for physics object.
     */
    pylith::materials::AuxiliaryFactoryPoroelastic* getAuxiliaryFactory(void);

    /** Add rheology subfields to auxiliary field.
     *
     * @param[inout] auxiliaryField Auxiliary field.
     */
    void addAuxiliarySubfields(void);

    // ============================= RHS ==================================== //

    // ---------------------------------------------------------------------------------------------------------------------
    // Select g0p function. Will only be used for the dynamic case.
    PetscPointFn* getKernelg0p(const spatialdata::geocoords::CoordSys* coordsys,
                               const bool _useBodyForce,
                               const bool _gravityField,
                               const bool _useSourceDensity) const;

    // ---------------------------------------------------------------------------------------------------------------------
    /** Get pressure kernel for RHS residual, G(t,s).
     *
     * @param[in] coordsys Coordinate system.
     *
     * @return RHS residual kernel for Darcy velocity.
     */
    PetscPointFn* getKernelg1p_explicit(const spatialdata::geocoords::CoordSys* coordsys,
                                        const bool _gravityField) const;

    // ---------------------------------------------------------------------------------------------------------------------
    // Get stress kernel for RHS residual, G(t,s)
    PetscPointFn* getKernelg1v_explicit(const spatialdata::geocoords::CoordSys* coordsys) const;

    // =============================== LHS =================================== //

    // ---------------------------------------------------------------------------------------------------------------------
    // Get variation in fluid content kernel for LHS residual, F(t,s,\dot{s})
    PetscPointFn* getKernelf0p_explicit(const spatialdata::geocoords::CoordSys* coordsys) const;

    // ---------------------------------------------------------------------------------------------------------------------
    // Select implicit f0p function.
    PetscPointFn* getKernelf0p_implicit(const spatialdata::geocoords::CoordSys* coordsys,
                                        const bool _useBodyForce,
                                        const bool _gravityField,
                                        const bool _useSourceDensity) const;

    // ---------------------------------------------------------------------------------------------------------------------
    // Get stress kernel for LHS residual, F(t,s,\dot{s})
    PetscPointFn* getKernelf1u_implicit(const spatialdata::geocoords::CoordSys* coordsys) const;

    // ---------------------------------------------------------------------------------------------------------------------
    /** Get pressure kernel for LHS residual.
     *
     * @param[in] coordsys Coordinate system.
     *
     * @return LHS residual kernel for Darcy velocity.
     */
    PetscPointFn* getKernelf1p_implicit(const spatialdata::geocoords::CoordSys* coordsys,
                                        const bool _useBodyForce,
                                        const bool _gravityField) const;

    // ---------------------------------------------------------------------------------------------------------------------
    // Get poroelastic constants kernel for LHS Jacobian
    PetscPointJacFn* getKernelJf3uu(const spatialdata::geocoords::CoordSys* coordsys) const;

    // ---------------------------------------------------------------------------------------------------------------------
    // Get biot coefficient kernel for LHS Jacobian
    PetscPointJacFn* getKernelJf2up(const spatialdata::geocoords::CoordSys* coordsys) const;

    // ---------------------------------------------------------------------------------------------------------------------
    // Get lambda kernel for LHS Jacobian
    PetscPointJacFn* getKernelJf2ue(const spatialdata::geocoords::CoordSys* coordsys) const;

    // ---------------------------------------------------------------------------------------------------------------------
    // Get Specific storage kernel for LHS Jacobian F(t,s, \dot{s}).
    PetscPointJacFn* getKernelJf0pp(const spatialdata::geocoords::CoordSys* coordsys) const;

    // ---------------------------------------------------------------------------------------------------------------------
    // Get Darcy Conductivity kernel for LHS Jacobian
    PetscPointJacFn* getKernelJf3pp(const spatialdata::geocoords::CoordSys* coordsys) const;

    // ---------------------------------------------------------------------------------------------------------------------
    // Get biot coefficient kernel for LHS Jacobian F(t,s, \dot{s}).
    PetscPointJacFn* getKernelJf0pe(const spatialdata::geocoords::CoordSys* coordsys) const;

    // ============================ DERIVED FIELDS ========================== //

    /** Get stress kernel for derived field.
     *
     * @param[in] coordsys Coordinate system.
     *
     * @return Project kernel for computing stress subfield in derived field.
     */
    PetscPointFn* getKernelCauchyStressVector(const spatialdata::geocoords::CoordSys* coordsys) const;

    /** Get water content kernel for derived field.
     *
     * @param[in] coordsys Coordinate system.
     *
     * @return Project kernel for computing water content subfield in derived field.
     */
    PetscPointFn* getKernelWaterContent(const spatialdata::geocoords::CoordSys* coordsys) const;

    /** Update kernel constants.
     *
     * @param[inout] kernelConstants Array of constants used in integration kernels.
     * @param[in] dt Current time step.
     */
    void updateKernelConstants(pylith::real_array* kernelConstants,
                               const PylithReal dt) const;

    /** Add kernels for updating state variables, implicit.
     *
     * @param[inout] kernels Array of kernels for updating state variables.
     * @param[in] coordsys Coordinate system.
     * @param[in] _useStateVars Update kernels?
     */
    void addKernelsUpdateStateVarsImplicit(std::vector<pylith::feassemble::IntegratorDomain::ProjectKernels>* kernels,
                                           const spatialdata::geocoords::CoordSys* coordsys,
                                           const bool _useStateVars) const;

    /** Add kernels for updating state variables, explicit.
     *
     * @param[inout] kernels Array of kernels for updating state variables.
     * @param[in] coordsys Coordinate system.
     * @param[in] _useStateVars Update kernels?
     */
    void addKernelsUpdateStateVarsExplicit(std::vector<pylith::feassemble::IntegratorDomain::ProjectKernels>* kernels,
                                           const spatialdata::geocoords::CoordSys* coordsys,
                                           const bool _useStateVars) const;

    // PROTECTED METHODS ///////////////////////////////////////////////////////////////////////////////////////////////
protected:

    /** Get auxiliary factory associated with physics.
     * @return Auxiliary factory for physics object.
     */
    pylith::feassemble::AuxiliaryFactory* _getAuxiliaryFactory(void);

    // PRIVATE MEMBERS /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    pylith::materials::AuxiliaryFactoryPoroelastic* _auxiliaryFactory; ///< Factory for auxiliary subfields.
    bool _useReferenceState; ///< Flag to use reference stress and strain.
    bool _useTensorPermeability; ///< Flag to use tensor permeability

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    IsotropicLinearPoroelasticity(const IsotropicLinearPoroelasticity&); ///< Not implemented.
    const IsotropicLinearPoroelasticity& operator=(const IsotropicLinearPoroelasticity&); ///< Not implemented

};

// class IsotropicLinearPoroelasticity

// End of file
