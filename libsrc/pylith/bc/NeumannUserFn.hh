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

/** @file libsrc/bc/NeumannUserFn.hh
 *
 * @brief C++ implementation of Neumann (prescribed values at
 * degrees of freedom) boundary condition using a user-provided function.
 */

#if !defined(pylith_bc_neumannuserfn_hh)
#define pylith_bc_neumannuserfn_hh

#include "pylith/bc/BoundaryCondition.hh" // ISA BoundaryCondition

#include "pylith/topology/topologyfwd.hh" // USES Field

class pylith::bc::NeumannUserFn : public pylith::bc::BoundaryCondition {
    friend class TestNeumannUserFn; // unit testing

    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    /// Default constructor.
    NeumannUserFn(void);

    /// Destructor.
    ~NeumannUserFn(void);

    /// Deallocate PETSc and local data structures.
    void deallocate(void);

    /** Set user function specifying field on boundary.
     *
     * @param[in] fn Function specifying field on boundary.
     */
    void setUserFn(PetscBdPointFunc fn);

    /** Get user function specifying field on boundary
     *
     * @preturns Function specifying field on boundary.
     */
    PetscBdPointFunc getUserFn(void) const;

    /** Create integrator and set kernels.
     *
     * @param[in] solution Solution field.
     * @returns Integrator if applicable, otherwise NULL.
     */
    pylith::feassemble::Integrator* createIntegrator(const pylith::topology::Field& solution);

    /** Create constraint and set kernels.
     *
     * @param[in] solution Solution field.
     * @returns Constraint if applicable, otherwise NULL.
     */
    std::vector<pylith::feassemble::Constraint*> createConstraints(const pylith::topology::Field& solution);

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

    // PRIVATE MEMBERS ////////////////////////////////////////////////////////////////////////////
private:

    PetscBdPointFunc _fn; ///< Kernel for boundary integration.

    // NOT IMPLEMENTED ////////////////////////////////////////////////////////////////////////////
private:

    NeumannUserFn(const NeumannUserFn&); ///< Not implemented.
    const NeumannUserFn& operator=(const NeumannUserFn&); ///< Not implemented.

}; // class NeumannUserFn

#endif // pylith_bc_neumannuserfn_hh

// End of file
