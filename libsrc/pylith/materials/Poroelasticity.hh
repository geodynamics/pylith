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

/** @file libsrc/materials/Poroelasticity.hh
 *
 * @brief C++ class for solving poroelasticity equation.
 */

#if !defined(pylith_materials_poroelasticity_hh)
#define pylith_materials_poroelasticity_hh

#include "materialsfwd.hh" // forward declarations

#include "pylith/materials/Material.hh" // ISA Material

class pylith::materials::Poroelasticity : public pylith::materials::Material
{
    friend class TestPoroelasticity; // unit testing

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:
    /// Default constructor.
    Poroelasticity(void);

    /// Destructor.
    ~Poroelasticity(void);

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

    /** Include source density?
     *
     * @param[in] value Flag indicating to include source density term.
     */
    void useSourceDensity(const bool value);

    /** Include source density?
     *
     * @returns True if including source density term, false otherwise.
     */
    bool useSourceDensity(void) const;

    /** Include constant pressure source?
     *
     * @param[in] value Flag indicating to include constant pressure source term.
     */
    void useConstantPressureSource(const bool value);

    /** Include source density?
     *
     * @returns True if including constant pressure source term, false otherwise.
     */
    bool useConstantPressureSource(void) const;

    /** Update fields?
     *
     * @param[in] value Flag indicating to update the auxiliary field values over time.
     */
    void useStateVars(const bool value);

    /** Update fields?
     *
     * @param[in] value Flag indicating to update the auxiliary field values over time.
     */
    bool useStateVars(void) const;

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

    /** Set bulk rheology.
     *
     * @param[in] rheology Bulk rheology for elasticity.
     */
    void setBulkRheology(pylith::materials::RheologyPoroelasticity *const rheology);

    /** Get bulk rheology.
     *
     * @returns Bulk rheology for elasticity.
     */
    pylith::materials::RheologyPoroelasticity *getBulkRheology(void) const;

    /** Verify configuration is acceptable.
     *
     * @param[in] solution Solution field.
     */
    void verifyConfiguration(const pylith::topology::Field &solution) const;

    /** Create integrator and set kernels.
     *
     * @param[in] solution Solution field.
     *
     *  @returns Integrator if applicable, otherwise NULL.
     */
    pylith::feassemble::Integrator *createIntegrator(const pylith::topology::Field &solution);

    /** Create auxiliary field.
     *
     * @param[in] solution Solution field.
     * @param[in\ domainMesh Finite-element mesh associated with integration domain.
     *
     * @returns Auxiliary field if applicable, otherwise NULL.
     */
    pylith::topology::Field *createAuxiliaryField(const pylith::topology::Field &solution,
                                                  const pylith::topology::Mesh &domainMesh);

    /** Create derived field.
     *
     * @param[in] solution Solution field.
     * @param[in\ domainMesh Finite-element mesh associated with integration domain.
     *
     * @returns Derived field if applicable, otherwise NULL.
     */
    pylith::topology::Field *createDerivedField(const pylith::topology::Field &solution,
                                                const pylith::topology::Mesh &domainMesh);

    // PROTECTED METHODS ///////////////////////////////////////////////////////////////////////////////////////////////
protected:
    /** Get auxiliary factory associated with physics.
     *
     * @return Auxiliary factory for physics object.
     */
    pylith::feassemble::AuxiliaryFactory *_getAuxiliaryFactory(void);

    /** Update kernel constants.
     *
     * @param[in] dt Current time step.
     */
    void _updateKernelConstants(const PylithReal dt);

    /** Get derived factory associated with physics.
     *
     * @return Derived factory for physics object.
     */
    pylith::topology::FieldFactory *_getDerivedFactory(void);

    // PRIVATE MEMBERS /////////////////////////////////////////////////////////////////////////////////////////////////
private:
    /** Set kernels for residual.
     *
     * @param[out] integrator Integrator for material.
     * @param[in] solution Solution field.
     */
    void _setKernelsResidual(pylith::feassemble::IntegratorDomain *integrator,
                             const pylith::topology::Field &solution) const;

    /** Set kernels for Jacobian.
     *
     * @param[out] integrator Integrator for material.
     * @param[in] solution Solution field.
     */
    void _setKernelsJacobian(pylith::feassemble::IntegratorDomain *integrator,
                             const pylith::topology::Field &solution) const;

    /** Set kernels for computing updated state variables in auxiliary field.
     *
     * @param[out] integrator Integrator for material.
     * @param[in] solution Solution field.
     */
    void _setKernelsUpdateStateVars(pylith::feassemble::IntegratorDomain *integrator,
                                    const pylith::topology::Field &solution) const;

    /** Set kernels for computing derived field.
     *
     * @param[out] integrator Integrator for material.
     * @param[in] solution Solution field.
     */
    void _setKernelsDerivedField(pylith::feassemble::IntegratorDomain *integrator,
                                 const pylith::topology::Field &solution) const;

    // PRIVATE MEMBERS /////////////////////////////////////////////////////////////////////////////////////////////////
private:
    bool _useBodyForce;                                           ///< Flag to include body force term.
    bool _useReferenceState;                                      ///< Flag to use reference stress and strain.
    bool _useSourceDensity;                                       ///< Flag to use source density.
    bool _useConstantPressureSource;                              ///< Flag to use constant pressure source.
    bool _useStateVars;                                           ///< Flag to update auxiliary fields.
    pylith::materials::RheologyPoroelasticity *_rheology;         ///< Bulk rheology for poroelasticity.
    pylith::materials::DerivedFactoryElasticity *_derivedFactory; ///< Factory for creating derived fields.

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:
    Poroelasticity(const Poroelasticity &);                  ///< Not implemented.
    const Poroelasticity &operator=(const Poroelasticity &); /// Not implemented.

}; // class Poroelasticity

#endif // pylith_materials_poroelasticity_hh

// End of file
