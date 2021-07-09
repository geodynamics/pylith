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
// Copyright (c) 2010-2021 University of California, Davis
//
// See LICENSE.md for license information.
//
// ----------------------------------------------------------------------
//

/** @file libsrc/materials/IsotropicLinearIncompElasticity.hh
 *
 * @brief C++ class for isotropic linear incompressible elasticity.
 */

#if !defined(pylith_materials_isotropiclinearincompelasticity_hh)
#define pylith_materials_isotropiclinearincompelasticity_hh

#include "materialsfwd.hh" // forward declarations

#include "pylith/materials/RheologyIncompressibleElasticity.hh" // ISA RheologyIncompressibleElasticity

class pylith::materials::IsotropicLinearIncompElasticity : public pylith::materials::RheologyIncompressibleElasticity {
    friend class TestIsotropicLinearIncompElasticity; // unit testing

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
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

    /** Get stress kernel for RHS residual, G(t,s).
     *
     * @param[in] coordsys Coordinate system.
     *
     * @return RHS residual kernel for stress.
     */
    PetscPointFunc getKernelResidualStress(const spatialdata::geocoords::CoordSys* coordsys) const;

    /** Get pressure kernel for RHS residual, G(t,s).
     *
     * @param[in] coordsys Coordinate system.
     *
     * @return RHS residual kernel for pressure.
     */
    PetscPointFunc getKernelResidualPressure(const spatialdata::geocoords::CoordSys* coordsys) const;

    /** Get elastic constants kernel for RHS Jacobian G(t,s).
     *
     * @param[in] coordsys Coordinate system.
     *
     * @return RHS Jacobian kernel for elastic constants.
     */
    PetscPointJac getKernelJacobianElasticConstants(const spatialdata::geocoords::CoordSys* coordsys) const;

    /** Get inverse of the bulk modulus kernel for RHS Jacobian G(t,s).
     *
     * @param[in] coordsys Coordinate system.
     *
     * @return RHS Jacobian kernel for inverse of bulk modulus.
     */
    PetscPointJac getKernelJacobianInverseBulkModulus(const spatialdata::geocoords::CoordSys* coordsys) const;

    /** Get stress kernel for derived field.
     *
     * @param[in] coordsys Coordinate system.
     *
     * @return Project kernel for computing stress subfield in derived field.
     */
    PetscPointFunc getKernelDerivedCauchyStress(const spatialdata::geocoords::CoordSys* coordsys) const;

    // PROTECTED METHODS ///////////////////////////////////////////////////////////////////////////////////////////////
protected:

    /** Get auxiliary factory associated with physics.
     * @return Auxiliary factory for physics object.
     */
    pylith::feassemble::AuxiliaryFactory* _getAuxiliaryFactory(void);

    // PRIVATE MEMBERS /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    pylith::materials::AuxiliaryFactoryElastic* _auxiliaryFactory; ///< Factory for auxiliary subfields.
    bool _useReferenceState; ///< Flag to use reference stress and strain.

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    IsotropicLinearIncompElasticity(const IsotropicLinearIncompElasticity&); ///< Not implemented.
    const IsotropicLinearIncompElasticity& operator=(const IsotropicLinearIncompElasticity&); ///< Not implemented

};

// class IsotropicLinearIncompElasticity

#endif // pylith_materials_isotropiclinearincompelasticity_hh

// End of file
