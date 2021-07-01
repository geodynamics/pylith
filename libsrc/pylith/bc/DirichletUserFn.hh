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

/** @file libsrc/bc/DirichletUserFn.hh
 *
 * @brief C++ implementation of Dirichlet (prescribed values at
 * degrees of freedom) boundary condition using a user-provided function.
 */

#if !defined(pylith_bc_dirichletuserfn_hh)
#define pylith_bc_dirichletuserfn_hh

#include "pylith/bc/BoundaryCondition.hh" // ISA BoundaryCondition

#include "pylith/utils/types.hh" // HASA PetscUserFieldFunc

#include "pylith/topology/topologyfwd.hh" // USES Field

class pylith::bc::DirichletUserFn : public pylith::bc::BoundaryCondition {
    friend class TestDirichletUserFn; // unit testing

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Default constructor.
    DirichletUserFn(void);

    /// Destructor.
    ~DirichletUserFn(void);

    /// Deallocate PETSc and local data structures.
    void deallocate(void);

    /** Set indices of constrained degrees of freedom at each location.
     *
     * Example: [0, 1] to apply forces to x and y degrees of freedom in
     * a Cartesian coordinate system.
     *
     * @param[in] dof Array of indices for constrained degrees of freedom.
     * @param[in] size Size of array
     */
    void setConstrainedDOF(const int* flags,
                           const int size);

    /** Get indices of constrained degrees of freedom.
     *
     * @returns Array of indices for constrained degrees of freedom.
     */
    const pylith::int_array& getConstrainedDOF(void) const;

    /** Set user function specifying field on boundary.
     *
     * @param[in] fn Function specifying field on boundary.
     */
    void setUserFn(PetscUserFieldFunc fn);

    /** Get user function specifying field on boundary
     *
     * @preturns Function specifying field on boundary.
     */
    PetscUserFieldFunc getUserFn(void) const;

    /** Set user function specifying time derivative of field on boundary.
     *
     * @param[in] fn Function specifying time derivative of field on boundary.
     */
    void setUserFnDot(PetscUserFieldFunc fn);

    /** Get user function specifying time derivative of field on boundary
     *
     * @preturns Function specifying time derivative of field on boundary.
     */
    PetscUserFieldFunc getUserFnDot(void) const;

    /** Verify configuration is acceptable.
     *
     * @param[in] solution Solution field.
     */
    void verifyConfiguration(const pylith::topology::Field& solution) const;

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
    pylith::feassemble::Constraint* createConstraint(const pylith::topology::Field& solution);

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

    // PROTECTED METHODS ///////////////////////////////////////////////////////////////////////////////////////////////
protected:

    /** Get auxiliary factory associated with physics.
     *
     * @return Auxiliary factory for physics object.
     */
    pylith::feassemble::AuxiliaryFactory* _getAuxiliaryFactory(void);

    // PRIVATE MEMBERS /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    int_array _constrainedDOF; ///< List of constrained degrees of freedom at each location.
    PetscUserFieldFunc _fn; ///< Function specifying field on boundary.
    PetscUserFieldFunc _fnDot; ///< Function specifying time derivative of field on boundary.

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    DirichletUserFn(const DirichletUserFn&); ///< Not implemented.
    const DirichletUserFn& operator=(const DirichletUserFn&); ///< Not implemented.

}; // class DirichletUserFn

#endif // pylith_bc_dirichletuserfn_hh

// End of file
