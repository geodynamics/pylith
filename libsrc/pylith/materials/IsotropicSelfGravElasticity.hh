// -*- C++ -*-
//
// ----------------------------------------------------------------------
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
// ----------------------------------------------------------------------
//

/** @file libsrc/materials/IsotropicSelfGravElasticity.hh
 *
 * @brief C++ class for isotropic linear self gravitating elasticity.
 */

#if !defined(pylith_materials_isotropicselfgravelasticity_hh)
#define pylith_materials_isotropicselfgravelasticity_hh

#include "materialsfwd.hh" // forward declarations

#include "pylith/materials/RheologySelfGravElasticity.hh" // ISA RheologySelfGravElasticity

class pylith::materials::IsotropicSelfGravElasticity : public pylith::materials::RheologySelfGravElasticity
{
    friend class TestIsotropicSelfGravElasticity; // unit testing

    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:
    /// Default constructor.
    IsotropicSelfGravElasticity(void);

    /// Destructor.
    ~IsotropicSelfGravElasticity(void);

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
    pylith::materials::AuxiliaryFactoryElasticity *getAuxiliaryFactory(void);

    /** Add rheology subfields to auxiliary field.
     *
     * @param[inout] auxiliaryField Auxiliary field.
     */
    void addAuxiliarySubfields(void);

    /** Get f0p kernel for LHS residual, F(t,s,\dot{s}).
     *
     * @param[in] coordsys Coordinate system.
     *
     * @return RHS residual kernel for potential.
     */
    PetscPointFunc getKernelf0p(const spatialdata::geocoords::CoordSys *coordsys) const;

    /** Get f1u kernel for LHS residual, F(t,s,\dot{s}).
     *
     * @param[in] coordsys Coordinate system.
     *
     * @return LHS residual kernel for stress.
     */
    PetscPointFunc getKernelf1u(const spatialdata::geocoords::CoordSys *coordsys) const;

    /** Get f1p kernel for LHS residual, F(t,s,\dot{s}).
     *
     * @param[in] coordsys Coordinate system.
     *
     * @return LHS residual kernel for stress.
     */
    PetscPointFunc getKernelf1p(const spatialdata::geocoords::CoordSys *coordsys) const;

    /** Get Jf0pp kernel for LHS Jacobian F(t,s,\dot{s}).
     *
     * @param[in] coordsys Coordinate system.
     *
     * @return LHS Jf0pp kernel.
     */
    PetscPointJac getKernelJf3pp(const spatialdata::geocoords::CoordSys *coordsys) const;

    /** Get Jf3uu kernel for LHS Jacobian F(t,s,\dot{s}).
     *
     * @param[in] coordsys Coordinate system.
     *
     * @return LHS Jacobian kernel for elastic constants.
     */
    PetscPointJac getKernelJf3uu(const spatialdata::geocoords::CoordSys *coordsys) const;

    /** Get stress kernel for derived field.
     *
     * @param[in] coordsys Coordinate system.
     *
     * @return Project kernel for computing stress subfield in derived field.
     */
    PetscPointFunc getKernelCauchyStressVector(const spatialdata::geocoords::CoordSys *coordsys) const;

    // PROTECTED METHODS //////////////////////////////////////////////////////////////////////////
protected:
    /** Get auxiliary factory associated with physics.
     * @return Auxiliary factory for physics object.
     */
    pylith::feassemble::AuxiliaryFactory *_getAuxiliaryFactory(void);

    // PRIVATE MEMBERS ////////////////////////////////////////////////////////////////////////////
private:
    pylith::materials::AuxiliaryFactoryElastic *_auxiliaryFactory; ///< Factory for auxiliary subfields.
    bool _useReferenceState;                                       ///< Flag to use reference stress and strain.

    // NOT IMPLEMENTED ////////////////////////////////////////////////////////////////////////////
private:
    IsotropicSelfGravElasticity(const IsotropicSelfGravElasticity &);                  ///< Not implemented.
    const IsotropicSelfGravElasticity &operator=(const IsotropicSelfGravElasticity &); ///< Not implemented

}; // class IsotropicSelfGravElasticity

#endif // pylith_materials_isotropicSelfGravelasticity_hh

// End of file
