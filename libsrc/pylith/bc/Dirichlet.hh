// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2016 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

/** @file libsrc/bc/Dirichlet.hh
 *
 * @brief C++ implementation of Dirichlet (prescribed values at
 * degrees of freedom) boundary conditions.
 */

#if !defined(pylith_bc_dirichletnew_hh)
#define pylith_bc_dirichletnew_hh

// Include directives ---------------------------------------------------
#include "BoundaryCondition.hh" // ISA BoundaryCondition
#include "pylith/feassemble/ConstraintPointwise.hh" // ISA ConstraintPointwise

#include "pylith/topology/topologyfwd.hh" // USES Field

// Dirichlet ----------------------------------------------------
/// @brief Dirichlet (prescribed values at degrees of freedom) boundary
/// conditions with points on a boundary.
class pylith::bc::Dirichlet :
    public pylith::bc::BoundaryCondition,
    public pylith::feassemble::ConstraintPointwise {

    friend class DirichletAuxiliaryFactory; // factory for auxiliary fields
    friend class TestDirichletNew;   // unit testing

    // PUBLIC METHODS /////////////////////////////////////////////////////
public:

    /// Default constructor.
    Dirichlet(void);

    /// Destructor.
    ~Dirichlet(void);

    /// Deallocate PETSc and local data structures.
    void deallocate(void);

    /** Verify configuration is acceptable.
     *
     * @param[in] solution Solution field.
     */
    void verifyConfiguration(const pylith::topology::Field& solution) const;

    /** Initialize boundary condition.
     *
     * @param[in] solution Solution field.
     */
    void initialize(const pylith::topology::Field& solution);

    /** Set constrained values in solution field.
     *
     * @param[out] solution Solution field.
     * @param[in] t Current time.
     */
    void setSolution(pylith::topology::Field* solution,
                     const double t);

    // PROTECTED METHODS //////////////////////////////////////////////////
protected:

    /** Setup auxiliary subfields (discretization and query fns).
     *
     * Create subfields in auxiliary fields (includes name of the field,
     * vector field type, discretization, and scale for
     * nondimensionalization) and set query functions for filling them
     * from a spatial database.
     *
     * @attention The order of the calls to subfieldAdd() must match the
     * order of the auxiliary fields in the FE kernels.
     *
     * @param[in] solution Solution field.
     */
    virtual
    void _auxFieldSetup(const pylith::topology::Field& solution) = 0;

    /** Set kernels for RHS residual G(t,s).
     *
     * Potentially, there are g0 and g1 kernels for each equation. If no
     * kernel is needed, then set the kernel function to NULL.
     *
     * @param solution Solution field.
     */
    virtual
    void _setFEKernelsConstraint(const pylith::topology::Field& solution) = 0;


    // PROTECTED MEMBERS //////////////////////////////////////////////////
protected:

    pylith::topology::Mesh* _boundaryMesh;   ///< Boundary mesh.
    PetscPointFunc _bcKernel; ///< Kernel for boundary condition value.
    pylith::topology::FieldBase::Description _description; ///< Description for constrained field.

    // NOT IMPLEMENTED ////////////////////////////////////////////////////
private:

    Dirichlet(const Dirichlet&); ///< Not implemented.
    const Dirichlet& operator=(const Dirichlet&); ///< Not implemented.

}; // class Dirichlet

#endif // pylith_bc_dirichletnew_hh


// End of file
