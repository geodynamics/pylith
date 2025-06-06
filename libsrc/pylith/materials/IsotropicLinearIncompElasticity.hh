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

#include "pylith/materials/RheologyIncompressibleElasticity.hh" // ISA RheologyIncompressibleElasticity

class pylith::materials::IsotropicLinearIncompElasticity : public pylith::materials::RheologyIncompressibleElasticity {
    friend class TestIsotropicLinearIncompElasticity; // unit testing

    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    /// Default constructor.
    IsotropicLinearIncompElasticity(void);

    /// Destructor.
    ~IsotropicLinearIncompElasticity(void);

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

    /** Get f0p kernel for LHS residual, F(t,s,\dot{s}).
     *
     * @param[in] coordsys Coordinate system.
     *
     * @return RHS residual kernel for pressure.
     */
    PetscPointFn* getKernelf0p(const spatialdata::geocoords::CoordSys* coordsys) const;

    /** Get f1u kernel for LHS residual, F(t,s,\dot{s}).
     *
     * @param[in] coordsys Coordinate system.
     *
     * @return LHS residual kernel for stress.
     */
    PetscPointFn* getKernelf1u(const spatialdata::geocoords::CoordSys* coordsys) const;

    /** Get Jf0pp kernel for LHS Jacobian F(t,s,\dot{s}).
     *
     * @param[in] coordsys Coordinate system.
     *
     * @return LHS Jf0pp kernel.
     */
    PetscPointJacFn* getKernelJf0pp(const spatialdata::geocoords::CoordSys* coordsys) const;

    /** Get Jf3uu kernel for LHS Jacobian F(t,s,\dot{s}).
     *
     * @param[in] coordsys Coordinate system.
     *
     * @return LHS Jacobian kernel for elastic constants.
     */
    PetscPointJacFn* getKernelJf3uu(const spatialdata::geocoords::CoordSys* coordsys) const;

    /** Get stress kernel for derived field.
     *
     * @param[in] coordsys Coordinate system.
     *
     * @return Project kernel for computing stress subfield in derived field.
     */
    PetscPointFn* getKernelCauchyStressVector(const spatialdata::geocoords::CoordSys* coordsys) const;

    // PROTECTED METHODS //////////////////////////////////////////////////////////////////////////
protected:

    /** Get auxiliary factory associated with physics.
     * @return Auxiliary factory for physics object.
     */
    pylith::feassemble::AuxiliaryFactory* _getAuxiliaryFactory(void);

    // PRIVATE MEMBERS ////////////////////////////////////////////////////////////////////////////
private:

    pylith::materials::AuxiliaryFactoryElastic* _auxiliaryFactory; ///< Factory for auxiliary subfields.
    bool _useReferenceState; ///< Flag to use reference stress and strain.

    // NOT IMPLEMENTED ////////////////////////////////////////////////////////////////////////////
private:

    IsotropicLinearIncompElasticity(const IsotropicLinearIncompElasticity&); ///< Not implemented.
    const IsotropicLinearIncompElasticity& operator=(const IsotropicLinearIncompElasticity&); ///< Not implemented

}; // class IsotropicLinearIncompElasticity

// End of file
