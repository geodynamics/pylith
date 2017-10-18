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

/** @file libsrc/bc/NeumannNew.hh
 *
 * @brief C++ implementation of Neumann (e.g., traction) boundary conditions.
 */

#if !defined(pylith_bc_neumannnew_hh)
#define pylith_bc_neumannnew_hh

// Include directives ---------------------------------------------------
#include "BoundaryConditionNew.hh" // ISA BoundaryCondition
#include "pylith/feassemble/IntegratorPointwise.hh" // ISA IntegratorPointwise

#include "pylith/topology/topologyfwd.hh" // USES Field

// NeumannNew ----------------------------------------------------
/// @brief Neumann (e.g., traction) boundary conditions.
class pylith::bc::NeumannNew :
    public BoundaryConditionNew,
    public pylith::feassemble::IntegratorPointwise { // class NeumannNew
    friend class TestNeumannNew;   // unit testing

    // PUBLIC METHODS /////////////////////////////////////////////////////
public:

    /// Default constructor.
    NeumannNew(void);

    /// Destructor.
    ~NeumannNew(void);

    /// Deallocate PETSc and local data structures.
    void deallocate(void);

    /** Initialize boundary condition.
     *
     * @param[in] solution Solution field.
     */
    void initialize(const pylith::topology::Field& solution);

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
    void _setFEKernelsRHSResidual(const topology::Field& solution) = 0;


    // PROTECTED MEMBERS //////////////////////////////////////////////////
protected:

    pylith::topology::Mesh* _boundaryMesh;   ///< Boundary mesh.
    pylith::topology::FieldBase::Description _description; ///< Description of field associated with BC.

    // NOT IMPLEMENTED ////////////////////////////////////////////////////
private:

    NeumannNew(const NeumannNew&); ///< Not implemented.
    const NeumannNew& operator=(const NeumannNew&); ///< Not implemented.

}; // class NeumannNew

#endif // pylith_bc_neumannnew_hh


// End of file
