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

#include "pylith/materials/RheologyElasticity.hh" // ISA RheologyElasticity

class pylith::materials::IsotropicLinearMaxwell : public pylith::materials::RheologyElasticity {
    friend class TestIsotropicLinearMaxwell; // unit testing

    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    /// Default constructor.
    IsotropicLinearMaxwell(void);

    /// Destructor.
    ~IsotropicLinearMaxwell(void);

    /// Deallocate PETSc and local data structures.
    void deallocate(void);

    /** Use reference stress and strain in computation of stress and
     * strain?
     *
     * @param[in] value Flag indicating to include reference stress and strain.
     */
    void useReferenceState(const bool value);

    /** Use reference stress and strain in computation of stress and
     * strain?
     *
     * @returns True if using reference stress and strain, false otherwise.
     */
    bool useReferenceState(void) const;

    /** Get auxiliary factory associated with physics.
     *
     * @return Auxiliary factory for physics object.
     */
    pylith::materials::AuxiliaryFactoryElasticity* getAuxiliaryFactory(void);

    /** Add rheology subfields to auxiliary field.
     *
     * @param[inout] auxiliaryField Auxiliary field.
     */
    void addAuxiliarySubfields(void);

    /** Get stress kernel for LHS residual, F(t,s,\dot{s}).
     *
     * @param[in] coordsys Coordinate system.
     *
     * @return LHS residual kernel for stress.
     */
    PetscPointFunc getKernelf1v(const spatialdata::geocoords::CoordSys* coordsys) const;

    /** Get elastic constants kernel for LHS Jacobian F(t,s,\dot{s}).
     *
     * @param[in] coordsys Coordinate system.
     *
     * @return LHS Jacobian kernel for elastic constants.
     */
    PetscPointJac getKernelJf3vu(const spatialdata::geocoords::CoordSys* coordsys) const;

    /** Get f0 kernel for LHS interface residual, F(t,s,dot{s}), for negative fault face.
     *
     * @param[in] coordsys Coordinate system.
     *
     * @return LHS residual f0 kernel.
     */
    PetscBdPointFunc getKernelf0Neg(const spatialdata::geocoords::CoordSys* coordsys) const;

    /** Get f0 kernel for LHS interface residual, F(t,s,dot{s}), for positive fault face.
     *
     * @param[in] coordsys Coordinate system.
     *
     * @return LHS residual f0 kernel.
     */
    PetscBdPointFunc getKernelf0Pos(const spatialdata::geocoords::CoordSys* coordsys) const;

    /** Get Cauchy stress kernel for derived field.
     *
     * @param[in] coordsys Coordinate system.
     *
     * @return Project kernel for computing stress subfield in derived field.
     */
    PetscPointFunc getKernelCauchyStressVector(const spatialdata::geocoords::CoordSys* coordsys) const;

    /** Add kernels for updating state variables.
     *
     * @param[inout] kernels Array of kernels for updating state variables.
     * @param[in] coordsys Coordinate system.
     */
    void addKernelsUpdateStateVars(std::vector<pylith::feassemble::IntegratorDomain::ProjectKernels>* kernels,
                                   const spatialdata::geocoords::CoordSys* coordsys) const;

    /** Update kernel constants.
     *
     * @param[inout] kernelConstants Array of constants used in integration kernels.
     * @param[in] dt Current time step.
     */
    void updateKernelConstants(pylith::real_array* kernelConstants,
                               const PylithReal dt) const;

    // PRIVATE MEMBERS ////////////////////////////////////////////////////////////////////////////
private:

    pylith::materials::AuxiliaryFactoryViscoelastic* _auxiliaryFactory; ///< Factory for creating auxiliary subfields.
    bool _useReferenceState; ///< Flag to use reference stress and strain.

    // NOT IMPLEMENTED ////////////////////////////////////////////////////////////////////////////
private:

    IsotropicLinearMaxwell(const IsotropicLinearMaxwell&); ///< Not implemented.
    const IsotropicLinearMaxwell& operator=(const IsotropicLinearMaxwell&); ///< Not implemented

}; // class IsotropicLinearMaxwell

// End of file
