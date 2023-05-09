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

/** @file libsrc/faults/FaultRheology.hh
 *
 * @brief C++ abstract base class for fault constitutive models.
 */

#if !defined(pylith_faults_faultrheology_hh)
#define pylith_faults_faultrheology_hh

#include "faultsfwd.hh" // forward declarations
#include "pylith/utils/PyreComponent.hh" // ISA PyreComponent

#include "pylith/topology/topologyfwd.hh" // USES Field
#include "pylith/utils/arrayfwd.hh" // USES std::vector
#include "pylith/feassemble/IntegratorDomain.hh" // USES IntegratorDomain::ProjectKernels

#include "spatialdata/geocoords/geocoordsfwd.hh" // USES Coordsys

#include "petscds.h" // USES PetscPointFunc, PetscPointJac

class pylith::faults::FaultRheology : public pylith::utils::PyreComponent {
    friend class TestStaticFriction; // unit testing

    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    /// Default constructor.
    FaultRheology(void);

    /// Destructor.
    virtual ~FaultRheology(void);

    /// Deallocate PETSc and local data structures.
    void deallocate(void);

    /** Get auxiliary factory associated with physics.
     *
     * @return Auxiliary factory for physics object.
     */
    virtual
    pylith::faults::AuxiliaryFactoryRheology* getAuxiliaryFactory(void) = 0;

    /// Add rheology subfields to auxiliary field.
    virtual
    void addAuxiliarySubfields(void) = 0;

    /** Get triggers for needing to compute the elastic constants for the RHS Jacobian.
     *
     * @returns Triggers for needing to recompute the RHS Jacobian.
     */
    int getLHSJacobianTriggers(void) const;

    /** Add kernels for updating state variables.
     *
     * @param[inout] kernels Array of kernels for updating state variables.
     * @param[in] coordsys Coordinate system.
     */
    virtual
    void addKernelsUpdateStateVars(std::vector<pylith::feassemble::IntegratorDomain::ProjectKernels>* kernels,
                                   const spatialdata::geocoords::CoordSys* coordsys) const;

    /** Update kernel constants.
     *
     * @param[inout] kernelConstants Array of constants used in integration kernels.
     * @param[in] dt Current time step.
     */
    virtual
    void updateKernelConstants(pylith::real_array* kernelConstants,
                               const PylithReal dt) const;

    // PROTECTED MEMBERS //////////////////////////////////////////////////////////////////////////

    int _lhsJacobianTriggers; ///< Triggers for needing to recompute the RHS Jacobian.

    // NOT IMPLEMENTED ////////////////////////////////////////////////////////////////////////////
private:

    FaultRheology(const FaultRheology&); ///< Not implemented.
    const FaultRheology& operator=(const FaultRheology&); /// Not implemented.

}; // class FaultRheology

#endif // pylith_faults_faultrheology_hh

// End of file
