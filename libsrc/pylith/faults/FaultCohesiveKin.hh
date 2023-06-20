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

/** @file libsrc/faults/FaultCohesiveKin.hh
 *
 * @brief C++ implementation for a fault surface with kinematic
 * (prescribed) slip implemented with cohesive elements.
 */

#if !defined(pylith_faults_faultcohesivekin_hh)
#define pylith_faults_faultcohesivekin_hh

#include "FaultCohesive.hh" // ISA FaultCohesive
#include "pylith/materials/Material.hh" // USES Material

#include <string> // HASA std::string
#include <map> // HASA std::map

class pylith::faults::FaultCohesiveKin : public pylith::faults::FaultCohesive {
    friend class TestFaultCohesiveKin; // unit testing

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Default constructor.
    FaultCohesiveKin(void);

    /// Destructor.
    ~FaultCohesiveKin(void);

    /// Deallocate PETSc and local data structures.
    void deallocate(void);

    /** Set kinematic earthquake ruptures.
     *
     * @param names Array of kinematic earthquake rupture names.
     * @param numNames Number of earthquake rupture names.
     * @param ruptures Array of kinematic earthquake ruptures.
     * @param numRuptures Number of earthquake ruptures.
     */
    void setEqRuptures(const char* const* names,
                       const int numNames,
                       KinSrc** ruptures,
                       const int numRuptures);

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

    /** Update slip related subfields in auxiliary field at beginning of time step.
     *
     * @param[out] auxiliaryField Auxiliary field.
     * @param[in] t Current time.
     * @param[in] bitSlipSubfields Slip subfields to update.
     */
    void _updateSlip(pylith::topology::Field* auxiliaryField,
                     const double t,
                     const int bitSlipSubfields);

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

    // PROTECTED TYPEDEFS /////////////////////////////////////////////////////////////////////////
protected:

    typedef std::map<std::string, KinSrc*> srcs_type;

    // PROTECTED MEMBERS //////////////////////////////////////////////////////////////////////////
protected:

    pylith::faults::AuxiliaryFactoryKinematic* _auxiliaryFactory; ///< Factory for auxiliary subfields.
    srcs_type _ruptures; ///< Array of kinematic earthquake ruptures.
    PetscVec _slipVecRupture; ///< PETSc local Vec to hold slip for one kinematic rupture.
    PetscVec _slipVecTotal; ///< PETSc local Vec to hold slip for all kinematic ruptures.

    // NOT IMPLEMENTED ////////////////////////////////////////////////////////////////////////////
private:

    FaultCohesiveKin(const FaultCohesiveKin&); ///< Not implemented
    const FaultCohesiveKin& operator=(const FaultCohesiveKin&); ///< Not implemented.

}; // class FaultCohesiveKin

#endif // pylith_faults_faultcohesivekin_hh

// End of file
