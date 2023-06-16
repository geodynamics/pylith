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

/** @file libsrc/bc/FaultCohesiveStub.hh
 *
 * @brief Minimal C++ implementation of FaultCohesive to allow unit tests of
 * other objects needing meshes with cohesive cells.
 */

#if !defined(pylith_faults_faultcohesivestub_hh)
#define pylith_faults_faultcohesivestub_hh

#include "pylith/testing/testingfwd.hh" // forward declarations

#include "pylith/faults/FaultCohesive.hh" // ISA FaultCohesive

class pylith::faults::FaultCohesiveStub : public pylith::faults::FaultCohesive {
    friend class TestFaultCohesiveStub; // unit testing

    // PUBLIC METHODS /////////////////////////////////////////////////////
public:

    /// Default constructor.
    FaultCohesiveStub(void);

    /// Destructor.
    ~FaultCohesiveStub(void);

    /** Verify configuration is acceptable.
     *
     * @param[in] solution Solution field.
     */
    void verifyConfiguration(const pylith::topology::Field& solution) const;

    /** Create auxiliary field.
     *
     * @param[in] solution Solution field.
     * @param[in\ physicsMesh Finite-element mesh associated with physics.
     *
     * @returns Auxiliary field if applicable, otherwise NULL.
     */
    pylith::topology::Field* createAuxiliaryField(const pylith::topology::Field& solution,
                                                  const pylith::topology::Mesh& physicsMesh);

    /** Create integrator and set kernels.
     *
     * @param[in] solution Solution field.
     * @param[in] materials Materials in problem.
     * @returns Integrator if applicable, otherwise NULL.
     */
    pylith::feassemble::Integrator* createIntegrator(const pylith::topology::Field& solution,
                                                     const std::vector<pylith::materials::Material*>& materials);

    // PROTECTED METHODS //////////////////////////////////////////////////////////////////////////
protected:

    /** Get auxiliary factory associated with physics.
     *
     * @return Auxiliary factory for physics object.
     */
    pylith::feassemble::AuxiliaryFactory* _getAuxiliaryFactory(void);

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

    // NOT IMPLEMENTED ////////////////////////////////////////////////////////////////////////////
private:

    FaultCohesiveStub(const FaultCohesiveStub&); ///< Not implemented.
    const FaultCohesiveStub& operator=(const FaultCohesiveStub&); ///< Not implemented.

}; // class FaultCohesiveStub

#endif // pylith_faults_faultcohesivestub_hh

// End of file
