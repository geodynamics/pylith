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

/** @file libsrc/materials/IsotropicLinearIncompElasticity.hh
 *
 * @brief C++ class for isotropic linear incompressible elasticity.
 */

#if !defined(pylith_materials_isotropiclinearincompelasticity_hh)
#define pylith_materials_isotropiclinearincompelasticity_hh

#include "materialsfwd.hh" // forward declarations

#include "pylith/materials/Material.hh" // ISA Material

class pylith::materials::IsotropicLinearIncompElasticity : public pylith::materials::Material {
    friend class TestIsotropicLinearIncompElasticity; // unit testing

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Default constructor.
    IsotropicLinearIncompElasticity(void);

    /// Destructor.
    ~IsotropicLinearIncompElasticity(void);

    /// Deallocate PETSc and local data structures.
    void deallocate(void);

    /** Include inertia?
     *
     * @param[in] value Flag indicating to include inertial term.
     */
    void useInertia(const bool value);

    /** Include inertia?
     *
     * @returns True if including inertial term, false otherwise.
     */
    bool useInertia(void) const;

    /** Include body force?
     *
     * @param[in] value Flag indicating to include body force term.
     */
    void useBodyForce(const bool value);

    /** Include body force?
     *
     * @returns True if including body force term, false otherwise.
     */
    bool useBodyForce(void) const;

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

    /** Verify configuration is acceptable.
     *
     * @param[in] solution Solution field.
     */
    void verifyConfiguration(const pylith::topology::Field& solution) const;

    /** Create integrator and set kernels.
     *
     * @param[in] solution Solution field.
     *
     *  @returns Integrator if applicable, otherwise NULL.
     */
    pylith::feassemble::Integrator* createIntegrator(const pylith::topology::Field& solution);

    /** Create auxiliary field.
     *
     * @param[in] solution Solution field.
     * @param[in\ domainMesh Finite-element mesh associated with integration domain.
     *
     * @returns Auxiliary field if applicable, otherwise NULL.
     */
    pylith::topology::Field* createAuxiliaryField(const pylith::topology::Field& solution,
                                                  const pylith::topology::Mesh& domainMesh);

    /** Create derived field.
     *
     * @param[in] solution Solution field.
     * @param[in\ domainMesh Finite-element mesh associated with integration domain.
     *
     * @returns Derived field if applicable, otherwise NULL.
     */
    pylith::topology::Field* createDerivedField(const pylith::topology::Field& solution,
                                                const pylith::topology::Mesh& domainMesh);

    // PROTECTED METHODS ///////////////////////////////////////////////////////////////////////////////////////////////
protected:

    /** Get auxiliary factory associated with physics.
     * @return Auxiliary factory for physics object.
     */
    pylith::feassemble::AuxiliaryFactory* _getAuxiliaryFactory(void);

    // PRIVATE MEMBERS /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    /** Set kernels for RHS residual.
     *
     * @param[out] integrator Integrator for material.
     * @param[in] solution Solution field.
     */
    void _setKernelsRHSResidual(pylith::feassemble::IntegratorDomain* integrator,
                                const pylith::topology::Field& solution) const;

    /** Set kernels for RHS Jacobian.
     *
     * @param[out] integrator Integrator for material.
     * @param[in] solution Solution field.
     */
    void _setKernelsRHSJacobian(pylith::feassemble::IntegratorDomain* integrator,
                                const pylith::topology::Field& solution) const;

    /** Set kernels for LHS residual.
     *
     * @param[out] integrator Integrator for material.
     * @param[in] solution Solution field.
     */
    void _setKernelsLHSResidual(pylith::feassemble::IntegratorDomain* integrator,
                                const pylith::topology::Field& solution) const;

    /** Set kernels for LHS Jacobian.
     *
     * @param[out] integrator Integrator for material.
     * @param[in] solution Solution field.
     */
    void _setKernelsLHSJacobian(pylith::feassemble::IntegratorDomain* integrator,
                                const pylith::topology::Field& solution) const;

    // PRIVATE MEMBERS /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    bool _useInertia; ///< Flag to include inertial term.
    bool _useBodyForce; ///< Flag to include body force term.
    bool _useReferenceState; ///< Flag to use reference stress and strain.

    pylith::materials::AuxiliaryFactoryElastic* _auxiliaryFactory; ///< Factory for auxiliary subfields.

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    IsotropicLinearIncompElasticity(const IsotropicLinearIncompElasticity&); ///< Not implemented.
    const IsotropicLinearIncompElasticity& operator=(const IsotropicLinearIncompElasticity&); ///< Not implemented

}; // class IsotropicLinearIncompElasticity

#endif // pylith_materials_isotropiclinearincompelasticity_hh

// End of file
