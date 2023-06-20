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

/** @file libsrc/faults/FaultCohesiveImpulses.hh
 *
 * @brief C++ implementation for a fault surface with kinematic
 * (prescribed) slip implemented with cohesive elements which we decompose into impulses for a Green function.
 */

#if !defined(pylith_faults_faultcohesiveimpulses_hh)
#define pylith_faults_faultcohesiveimpulses_hh

#include "FaultCohesive.hh" // ISA FaultCohesive

class pylith::faults::FaultCohesiveImpulses : public pylith::faults::FaultCohesive {
    friend class TestFaultCohesiveImpulses; // unit testing

    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    /// Default constructor.
    FaultCohesiveImpulses(void);

    /// Destructor.
    ~FaultCohesiveImpulses(void);

    /// Deallocate PETSc and local data structures.
    void deallocate(void);

    /** Set indices of fault degrees of freedom associated with
     * impulses.
     *
     * @param flags Array of indices for degrees of freedom.
     * @param size Size of array
     */
    void setImpulseDOF(const int* flags,
                       const size_t size);

    /** Set threshold for nonzero impulse amplitude.
     *
     * @param value Threshold for detecting nonzero amplitude.
     */
    void setThreshold(const double value);

    /** Get the total number of impulses that will be applied on this process.
     *
     * @returns Number of impulses.
     */
    size_t getNumImpulsesLocal(void);

    /** Verify configuration is acceptable.
     *
     * @param[in] solution Solution field.
     */
    void verifyConfiguration(const pylith::topology::Field& solution) const;

    /** Create auxiliary field.
     *
     * @param[in] solution Solution field.
     * @param[in] domainMesh Finite-element mesh associated with integration domain.
     *
     * @returns Auxiliary field if applicable, otherwise NULL.
     */
    pylith::topology::Field* createAuxiliaryField(const pylith::topology::Field& solution,
                                                  const pylith::topology::Mesh& domainMesh);

    /** Update auxiliary subfields at beginning of time step.
     *
     * @param[out] auxiliaryField Auxiliary field.
     * @param[in] t Current time.
     */
    void updateAuxiliaryField(pylith::topology::Field* auxiliaryField,
                              const double t);

    // PROTECTED METHODS //////////////////////////////////////////////////////////////////////////
protected:

    /** Get auxiliary factory associated with physics.
     *
     * @return Auxiliary factory for physics object.
     */
    pylith::feassemble::AuxiliaryFactory* _getAuxiliaryFactory(void);

    /** Update slip subfield in auxiliary field at beginning of time step.
     *
     * @param[out] auxiliaryField Auxiliary field.
     * @param[in] impulseIndex Index of impulse (-1 indicates no impulse is applied on this process).
     */
    void _updateSlip(pylith::topology::Field* auxiliaryField,
                     const long impulseIndex);

    /** Set kernels for residual.
     *
     * @param[out] integrator Integrator for material.
     * @param[in] solution Solution field.
     * @param[in] materials Materials in problem.
     */
    void _setKernelsResidual(pylith::feassemble::IntegratorInterface* integrator,
                             const pylith::topology::Field& solution,
                             const std::vector<pylith::materials::Material*>& materials) const;

    /** Set kernels for Jacobian.
     *
     * @param[out] integrator Integrator for material.
     * @param[in] solution Solution field.
     * @param[in] materials Materials in problem.
     */
    void _setKernelsJacobian(pylith::feassemble::IntegratorInterface* integrator,
                             const pylith::topology::Field& solution,
                             const std::vector<pylith::materials::Material*>& materials) const;

    // PRIVATE METHODS ////////////////////////////////////////////////////////////////////////////
private:

    pylith::faults::AuxiliaryFactoryKinematic* _auxiliaryFactory; ///< Factory for auxiliary subfields.
    PylithReal _threshold; ///< Threshold for nonzero impulse amplitude.
    int_array _impulseDOF; ///< Degrees of freedom with impulses.
    int_array _impulsePoints; ///< Points with nonzero threshold.

    // NOT IMPLEMENTED ////////////////////////////////////////////////////////////////////////////
private:

    FaultCohesiveImpulses(const FaultCohesiveImpulses&); ///< Not implemented
    const FaultCohesiveImpulses& operator=(const FaultCohesiveImpulses&); ///< Not implemented.

}; // class FaultCohesiveImpulses

#endif // pylith_faults_faultcohesiveimpulses_hh

// End of file
