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

/** @file libsrc/sources/SourceTimeFunctionPointForce.hh
 *
 * @brief C++ abstract base class for bulk rheologies associated with the pointforce equation.
 */

#if !defined(pylith_sources_sourcetimefunctionpointforce_hh)
#define pylith_sources_sourcetimefunctionpointforce_hh

#include "sourcesfwd.hh" // forward declarations
#include "pylith/utils/PyreComponent.hh" // ISA PyreComponent

#include "pylith/topology/topologyfwd.hh" // USES Field
#include "pylith/utils/arrayfwd.hh" // USES std::vector
#include "pylith/feassemble/IntegratorDomain.hh" // USES IntegratorDomain::ProjectKernels

#include "spatialdata/geocoords/geocoordsfwd.hh" // USES Coordsys

#include "petscds.h" // USES PetscPointFunc, PetscPointJac

class pylith::sources::SourceTimeFunctionPointForce : public pylith::utils::PyreComponent {
    friend class TestRickerFunction; // unit testing

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Default constructor.
    SourceTimeFunctionPointForce(void);

    /// Destructor.
    virtual ~SourceTimeFunctionPointForce(void);

    /// Deallocate PETSc and local data structures.
    void deallocate(void);

    /** Get auxiliary factory associated with physics.
     *
     * @return Auxiliary factory for physics object.
     */
    virtual
    pylith::sources::AuxiliaryFactoryPointForce* getAuxiliaryFactory(void) = 0;

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

    // PROTECTED MEMBERS ///////////////////////////////////////////////////////////////////////////////////////////////

    int _JacobianTriggers; ///< Triggers for needing to recompute the RHS Jacobian.

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    SourceTimeFunctionPointForce(const SourceTimeFunctionPointForce&); ///< Not implemented.
    const SourceTimeFunctionPointForce& operator=(const SourceTimeFunctionPointForce&); /// Not implemented.

}; // class SourceTimeFunctionPointForce

#endif // pylith_sources_sourcetimefunctionpointforce_hh

// End of file