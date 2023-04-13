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

/** @file libsrc/materials/IsotropicDruckerPragerElastoplastic.hh
 *
 * @brief C++ class for isotropic Drucker-Prager elastoplastic material.
 */

#if !defined(pylith_materials_isotropicdruckerpragerelastoplastic_hh)
#define pylith_materials_isotropicdruckerpragerelastoplastic_hh

#include "materialsfwd.hh" // forward declarations

#include "pylith/materials/RheologyElasticity.hh" // ISA RheologyElasticity

class pylith::materials::IsotropicDruckerPragerElastoplastic : public pylith::materials::RheologyElasticity
{
    friend class TestIsotropicDruckerPragerElastoplastic; // unit testing

    // PUBLIC ENUMS   //////////////////////////////////////////////////////////////////////////////////////////////////
public:
    enum FitMohrCoulombEnum
    {
        MOHR_COULOMB_CIRCUMSCRIBED = 0,
        MOHR_COULOMB_MIDDLE = 1,
        MOHR_COULOMB_INSCRIBED = 2,
    }; // FitMohrCoulombType

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:
    /// Default constructor.
    IsotropicDruckerPragerElastoplastic(void);

    /// Destructor.
    ~IsotropicDruckerPragerElastoplastic(void);

    /** Set fit to Mohr-Coulomb surface.
     *
     * @param[in] value Mohr-Coulomb surface match type.
     */
    void fitMohrCoulomb(FitMohrCoulombEnum value);

    /**  Get fit type to Mohr-Coulomb surface.
     *
     * @returns Fit type for Mohr-Coulomb surface match.
     */
    int fitMohrCoulomb(void) const;

    /** Set flag for whether to allow tensile yield.
     *
     * @param[in] value True if tensile yield is allowed.
     */
    void allowTensileYield(const bool value);

    /**  Get flag for whether to allow tensile yield.
     *
     * @returns True if tensile yield is allowed.
     */
    int allowTensileYield(void) const;

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
    pylith::materials::AuxiliaryFactoryElasticity *getAuxiliaryFactory(void);

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
    PetscPointFunc getKernelf1v(const spatialdata::geocoords::CoordSys *coordsys) const;

    /** Get elastic constants kernel for LHS Jacobian F(t,s,\dot{s}).
     *
     * @param[in] coordsys Coordinate system.
     *
     * @return LHS Jacobian kernel for elastic constants.
     */
    PetscPointJac getKernelJf3vu(const spatialdata::geocoords::CoordSys *coordsys) const;

    /** Get stress kernel for derived field.
     *
     * @param[in] coordsys Coordinate system.
     *
     * @return Project kernel for computing stress subfield in derived field.
     */
    PetscPointFunc getKernelCauchyStressVector(const spatialdata::geocoords::CoordSys *coordsys) const;

    /** Add kernels for updating state variables.
     *
     * @param[inout] kernels Array of kernels for updating state variables.
     * @param[in] coordsys Coordinate system.
     */
    void addKernelsUpdateStateVars(std::vector<pylith::feassemble::IntegratorDomain::ProjectKernels> *kernels,
                                   const spatialdata::geocoords::CoordSys *coordsys) const;

    /** Update kernel constants.
     *
     * @param[inout] kernelConstants Array of constants used in integration kernels.
     * @param[in] dt Current time step.
     */
    void updateKernelConstants(pylith::real_array *kernelConstants,
                               const PylithReal dt) const;

    // PRIVATE MEMBERS /////////////////////////////////////////////////////////////////////////////////////////////////
private:
    pylith::materials::AuxiliaryFactoryElastoplastic *_auxiliaryFactory; ///< Factory for creating auxiliary subfields.
    bool _useReferenceState;                                             ///< Flag to use reference stress and strain.
    int _fitMohrCoulomb;                                                 ///< Type of fit to Mohr-Coulomb yield surface.
    int _allowTensileYield;                                              ///< Flag for whether to allow tensile yield.

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:
    IsotropicDruckerPragerElastoplastic(const IsotropicDruckerPragerElastoplastic &);                  ///< Not implemented.
    const IsotropicDruckerPragerElastoplastic &operator=(const IsotropicDruckerPragerElastoplastic &); ///< Not implemented

}; // class IsotropicDruckerPragerElastoplastic

#endif // pylith_materials_isotropicdruckerpragerelastoplastic_hh

// End of file
