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

/** @file libsrc/bc/DirichletNew.hh
 *
 * @brief C++ implementation of Dirichlet (prescribed values at
 * degrees of freedom) boundary conditions.
 */

#if !defined(pylith_bc_dirichletnew_hh)
#define pylith_bc_dirichletnew_hh

// Include directives ---------------------------------------------------
#include "BoundaryConditionNew.hh" // ISA BoundaryCondition
#include "pylith/feassemble/ConstraintPointwise.hh" // ISA ConstraintPointwise

#include "pylith/topology/topologyfwd.hh" // USES Field

// DirichletNew ----------------------------------------------------
/// @brief Dirichlet (prescribed values at degrees of freedom) boundary
/// conditions with points on a boundary.
class pylith::bc::DirichletNew :
    public BoundaryConditionNew,
    public pylith::feassemble::ConstraintPointwise
{ // class DirichletNew
    friend class TestDirichletNew;   // unit testing

    // PUBLIC METHODS /////////////////////////////////////////////////////
public:

    /// Default constructor.
    DirichletNew(void);

    /// Destructor.
    ~DirichletNew(void);

    /// Deallocate PETSc and local data structures.
    void deallocate(void);

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
    void setValues(pylith::topology::Field* solution,
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
     */
    virtual
    void _auxFieldsSetup(void) = 0;

    /** Set kernels for RHS residual G(t,s).
     *
     * Potentially, there are g0 and g1 kernels for each equation. If no
     * kernel is needed, then set the kernel function to NULL.
     *
     * @param solution Solution field.
     */
    virtual
    void _setFEKernelsConstraint(const topology::Field& solution) = 0;


    // PROTECTED MEMBERS //////////////////////////////////////////////////
protected:

    pylith::topology::Mesh* _boundaryMesh;   ///< Boundary mesh.
    PetscPointFunc _bcKernel; ///< Kernel for boundary condition value.
    int _spaceDim; ///< Spatial dimension of domain.
    pylith::topology::FieldBase::VectorFieldEnum _vectorFieldType; ///< Vector field type of constrainted field.

    // NOT IMPLEMENTED ////////////////////////////////////////////////////
private:

    DirichletNew(const DirichletNew&); ///< Not implemented.
    const DirichletNew& operator=(const DirichletNew&); ///< Not implemented.

}; // class DirichletNew

#endif // pylith_bc_dirichletnew_hh


// End of file
