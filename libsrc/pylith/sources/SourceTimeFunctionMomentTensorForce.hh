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

/** @file libsrc/sources/SourceTimeFunctionMomentTensorForce.hh
 *
 * @brief C++ abstract base class for bulk rheologies associated with the momenttensorforce equation.
 */

#if !defined(pylith_sources_sourcetimefunctionmomenttensorforce_hh)
#define pylith_sources_sourcetimefunctionmomenttensorforce_hh

#include "sourcesfwd.hh" // forward declarations
#include "pylith/utils/PyreComponent.hh" // ISA PyreComponent

#include "pylith/topology/topologyfwd.hh" // USES Field
#include "pylith/utils/arrayfwd.hh" // USES std::vector
#include "pylith/feassemble/IntegratorDomain.hh" // USES IntegratorDomain::ProjectKernels

#include "spatialdata/geocoords/geocoordsfwd.hh" // USES Coordsys

#include "petscds.h" // USES PetscPointFunc, PetscPointJac

class pylith::sources::SourceTimeFunctionMomentTensorForce : public pylith::utils::PyreComponent {
    friend class TestSourceTimeFunctionMomentTensorForce; // unit testing

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Default constructor.
    SourceTimeFunctionMomentTensorForce(void);

    /// Destructor.
    virtual ~SourceTimeFunctionMomentTensorForce(void);

    /// Deallocate PETSc and local data structures.
    void deallocate(void);

    /** Get auxiliary factory associated with physics.
     *
     * @return Auxiliary factory for physics object.
     */
    virtual
    pylith::sources::AuxiliaryFactoryMomentTensorForce* getAuxiliaryFactory(void) = 0;

    /// Add sourcetimefunctionfunction subfields to auxiliary field.
    virtual
    void addAuxiliarySubfields(void) = 0;

    /** Get g1v kernel for residual, G(t,s).
     *
     * @param[in] coordsys Coordinate system.
     *
     * @return RHS residual kernel for stress.
     */
    virtual
    PetscPointFunc getKernelg1v_explicit(const spatialdata::geocoords::CoordSys* coordsys) const = 0;

    /** Get triggers for needing to compute the elastic constants for the RHS Jacobian.
     *
     * @returns Triggers for needing to recompute the RHS Jacobian.
     */
    int getJacobianTriggers(void) const;

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

    /** Update auxiliary field for current time.
     *
     * @param[inout] auxiliaryField Auxiliary field to update.
     * @param[in] t Current time.
     * @param[in] timeScale Time scale for nondimensionalization.
     */
    virtual
    void updateAuxiliaryField(pylith::topology::Field* auxiliaryField,
                              const PylithReal t,
                              const PylithReal timeScale);

    // PROTECTED MEMBERS ///////////////////////////////////////////////////////////////////////////////////////////////

    int _JacobianTriggers; ///< Triggers for needing to recompute the RHS Jacobian.

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    SourceTimeFunctionMomentTensorForce(const SourceTimeFunctionMomentTensorForce&); ///< Not implemented.
    const SourceTimeFunctionMomentTensorForce& operator=(const SourceTimeFunctionMomentTensorForce&); /// Not

    /// implemented.

}; // class SourceTimeFunctionMomentTensorForce

#endif // pylith_sources_sourcetimefunctionmomenttensorforce_hh

// End of file
