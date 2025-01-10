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

#include "pylith/faults/FaultCohesive.hh" // ISA FaultCohesive
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

    /** Set kernels for computing derived field.
     *
     * @param[out] integrator Integrator for material.
     * @param[in] solution Solution field.
     */
    void _setKernelsDerivedField(pylith::feassemble::IntegratorInterface* integrator,
                                 const pylith::topology::Field& solution) const;

    // PROTECTED TYPEDEFS /////////////////////////////////////////////////////////////////////////
protected:

    typedef std::map<std::string, KinSrc*> srcs_type;

    // PROTECTED MEMBERS //////////////////////////////////////////////////////////////////////////
protected:

    srcs_type _ruptures; ///< Array of kinematic earthquake ruptures.
    PetscVec _slipVecRupture; ///< PETSc local Vec to hold slip for one kinematic rupture.
    PetscVec _slipVecTotal; ///< PETSc local Vec to hold slip for all kinematic ruptures.

    // NOT IMPLEMENTED ////////////////////////////////////////////////////////////////////////////
private:

    FaultCohesiveKin(const FaultCohesiveKin&); ///< Not implemented
    const FaultCohesiveKin& operator=(const FaultCohesiveKin&); ///< Not implemented.

}; // class FaultCohesiveKin

// End of file
