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

/** @file libsrc/materials/Elasticity.hh
 *
 * @brief C++ class for solving elasticity equation.
 */

#if !defined(pylith_materials_elasticity_hh)
#define pylith_materials_elasticity_hh

#include "materialsfwd.hh" // forward declarations

#include "pylith/materials/Material.hh" // ISA Material

class pylith::materials::Elasticity : public pylith::materials::Material {
    friend class TestElasticity; // unit testing

    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    /// Default constructor.
    Elasticity(void);

    /// Destructor.
    ~Elasticity(void);

    /// Deallocate PETSc and local data structures.
    void deallocate(void);

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

    /** Set bulk rheology.
     *
     * @param[in] rheology Bulk rheology for elasticity.
     */
    void setBulkRheology(pylith::materials::RheologyElasticity* const rheology);

    /** Get bulk rheology.
     *
     * @returns Bulk rheology for elasticity.
     */
    pylith::materials::RheologyElasticity* getBulkRheology(void) const;

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

    /** Get default PETSc solver options appropriate for material.
     *
     * @param[in] isParallel True if running in parallel, False if running in serial.
     * @param[in] hasFault True if problem has fault, False otherwise.
     * @returns PETSc solver options.
     */
    pylith::utils::PetscOptions* getSolverDefaults(const bool isParallel,
                                                   const bool hasFault) const;

    /** Get residual kernels for an interior interface bounding material.
     *
     * @param[in] solution Solution field.
     * @param[in] face Side of interior interface for kernels.
     * @returns Array of residual kernels for interior interface.
     */
    std::vector<InterfaceResidualKernels> getInterfaceKernelsResidual(const pylith::topology::Field& solution,
                                                                      pylith::feassemble::IntegratorInterface::FaceEnum face) const;

    /** Get Jacobian kernels for an interior interface bounding material.
     *
     * @param[in] solution Solution field.
     * @param[in] face Side of interior interface for kernels.
     * @returns Array of Jacobian kernels for interior interface.
     */
    std::vector<InterfaceJacobianKernels> getInterfaceKernelsJacobian(const pylith::topology::Field& solution,
                                                                      pylith::feassemble::IntegratorInterface::FaceEnum face) const;

    // PROTECTED METHODS //////////////////////////////////////////////////////////////////////////
protected:

    /** Get auxiliary factory associated with physics.
     *
     * @return Auxiliary factory for physics object.
     */
    pylith::feassemble::AuxiliaryFactory* _getAuxiliaryFactory(void);

    /** Update kernel constants.
     *
     * @param[in] dt Current time step.
     */
    void _updateKernelConstants(const PylithReal dt);

    /** Get derived factory associated with physics.
     *
     * @return Derived factory for physics object.
     */
    pylith::topology::FieldFactory* _getDerivedFactory(void);

    // PRIVATE MEMBERS ////////////////////////////////////////////////////////////////////////////
private:

    /** Set kernels for residual.
     *
     * @param[out] integrator Integrator for material.
     * @param[in] solution Solution field.
     */
    void _setKernelsResidual(pylith::feassemble::IntegratorDomain* integrator,
                             const pylith::topology::Field& solution) const;

    /** Set kernels for Jacobian.
     *
     * @param[out] integrator Integrator for material.
     * @param[in] solution Solution field.
     */
    void _setKernelsJacobian(pylith::feassemble::IntegratorDomain* integrator,
                             const pylith::topology::Field& solution) const;

    /** Set kernels for computing updated state variables in auxiliary field.
     *
     * @param[out] integrator Integrator for material.
     * @param[in] solution Solution field.
     */
    void _setKernelsUpdateStateVars(pylith::feassemble::IntegratorDomain* integrator,
                                    const pylith::topology::Field& solution) const;

    /** Set kernels for computing derived field.
     *
     * @param[out] integrator Integrator for material.
     * @param[in] solution Solution field.
     */
    void _setKernelsDerivedField(pylith::feassemble::IntegratorDomain* integrator,
                                 const pylith::topology::Field& solution) const;

    // PRIVATE MEMBERS ////////////////////////////////////////////////////////////////////////////
private:

    bool _useBodyForce; ///< Flag to include body force term.
    pylith::materials::RheologyElasticity* _rheology; ///< Bulk rheology for elasticity.
    pylith::materials::DerivedFactoryElasticity* _derivedFactory; ///< Factory for creating derived fields.

    // NOT IMPLEMENTED ////////////////////////////////////////////////////////////////////////////
private:

    Elasticity(const Elasticity&); ///< Not implemented.
    const Elasticity& operator=(const Elasticity&); /// Not implemented.

};

// class Elasticity

#endif // pylith_materials_elasticity_hh

// End of file
