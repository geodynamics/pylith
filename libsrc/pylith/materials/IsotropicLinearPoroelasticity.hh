// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

/** @file libsrc/materials/IsotropicLinearPoroelasticity.hh
 *
 * @brief C++ class for isotropic linear poroelasticity.
 */

#if !defined(pylith_materials_isotropiclinearporoelasticity_hh)
#define pylith_materials_isotropiclinearporoelasticity_hh

#include "materialsfwd.hh" // forward declarations

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

    /** Get displacement kernel for RHS residual, G(t,s).
     *
     * @param[in] coordsys Coordinate system.
     *
     * @return RHS residual kernel for stress.
     */
    PetscPointFunc getKernelg1u(const spatialdata::geocoords::CoordSys* coordsys) const;

    /** Get velocity kernel for RHS residual, G(t,s).
     *
     * @param[in] coordsys Coordinate system.
     *
     * @return RHS residual kernel for stress.
     */
    PetscPointFunc getKernelg1v(const spatialdata::geocoords::CoordSys* coordsys) const;

    /** Get pressure kernel for RHS residual, G(t,s).
     *
     * @param[in] coordsys Coordinate system.
     *
     * @return RHS residual kernel for Darcy velocity.
     */
    PetscPointFunc getKernelg1p(const spatialdata::geocoords::CoordSys* coordsys, const bool _gravityField) const;

    /** Get displacement-displacement kernel for RHS Jacobian G(t,s).
     *
     * @param[in] coordsys Coordinate system.
     *
     * @return RHS Jacobian kernel for elastic constants.
     */
    PetscPointJac getKernelJg3uu(const spatialdata::geocoords::CoordSys* coordsys) const;

    /** Get velocity-displacement kernel for RHS Jacobian G(t,s).
     *
     * @param[in] coordsys Coordinate system.
     *
     * @return RHS Jacobian kernel for elastic constants.
     */
    PetscPointJac getKernelJg3vu(const spatialdata::geocoords::CoordSys* coordsys) const;

    /** Get displacement-pressure kernel for RHS Jacobian G(t,s).
     *
     * @param[in] coordsys Coordinate system.
     *
     * @return RHS Jacobian kernel for Biot Coefficient.
     */
    PetscPointJac getKernelJg2up(const spatialdata::geocoords::CoordSys* coordsys) const;

    /** Get velocity-pressure for RHS Jacobian G(t,s).
     *
     * @param[in] coordsys Coordinate system.
     *
     * @return RHS Jacobian kernel for Biot Coefficient.
     */
    PetscPointJac getKernelJg2vp(const spatialdata::geocoords::CoordSys* coordsys) const;

    /** Get displacement-trace strain for RHS Jacobian G(t,s).
     *
     * @param[in] coordsys Coordinate system.
     *
     * @return RHS Jacobian kernel for lambda.
     */
    PetscPointJac getKernelJg2ue(const spatialdata::geocoords::CoordSys* coordsys) const;

    /** Get pressure-pressure kernel for derived field.
     *
     * @param[in] coordsys Coordinate system.
     *
     * @return Project kernel for computing stress subfield in derived field.
     */
    PetscPointJac getKernelJg3pp(const spatialdata::geocoords::CoordSys* coordsys) const;

    // ============================= LHS ==================================== //

    /** Get variation in fluid content for LHS residual, F(t,s,\dot{s})
     *
     * @param[in] coordsys Coordinate system.
     *
     * @return LHS residual kernel for variation in fluid contenty.
     */
    PetscPointFunc getKernelf0p(const spatialdata::geocoords::CoordSys* coordsys, const bool _useInertia) const;

    /** Get specific storage kernel for LHS Jacobian G(t,s).
     *
     * @param[in] coordsys Coordinate system.
     *
     * @return RHS Jacobian kernel for tshift * 1/M (Jf0pp)
     */
    PetscPointJac getKernelJf0pp(const spatialdata::geocoords::CoordSys* coordsys) const;

    /** Get biot coefficient for LHS residual, F(t,s,\dot{s})
     *
     * @param[in] coordsys Coordinate system.
     *
     * @return LHS jacobian kernel for biot coefficient.
     */
    PetscPointJac getKernelJf0pe(const spatialdata::geocoords::CoordSys* coordsys) const;

    // ============================ DERIVED FIELDS ========================== //

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

    pylith::materials::AuxiliaryFactoryPoroelastic* _auxiliaryFactory; ///< Factory for auxiliary subfields.
    bool _useReferenceState; ///< Flag to use reference stress and strain.
    bool _useTensorPermeability; ///< Flag to use tensor permeability


    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    IsotropicLinearPoroelasticity(const IsotropicLinearPoroelasticity&); ///< Not implemented.
    const IsotropicLinearPoroelasticity& operator=(const IsotropicLinearPoroelasticity&); ///< Not implemented

};

// class IsotropicLinearPoroelasticity

#endif // pylith_materials_isotropiclinearporoelasticity_hh

// End of file
