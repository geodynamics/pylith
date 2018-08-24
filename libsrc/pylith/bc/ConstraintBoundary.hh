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

/** @file libsrc/bc/ConstraintBoundary.hh
 *
 * @brief C++ implementation of ConstraintBoundary for constraint degrees of freedom on external boundary.
 */

#if !defined(pylith_bc_constraintboundary_hh)
#define pylith_bc_constraintboundary_hh

// Include directives ---------------------------------------------------
#include "bcfwd.hh" // Forward declaration
#include "pylith/feassemble/ConstraintPointwise.hh" // ISA ConstraintPointwise

#include "pylith/topology/topologyfwd.hh" // USES Field

#include <string> // HASA std::string

// ConstraintBoundary ----------------------------------------------------
/// @brief Integrator for integrals over domain boundaries.
class pylith::bc::ConstraintBoundary : public pylith::feassemble::ConstraintPointwise {
    friend class TestConstraintBoundary; // unit testing

    // PUBLIC METHODS /////////////////////////////////////////////////////
public:

    /// Default constructor.
    ConstraintBoundary(void);

    /// Destructor.
    ~ConstraintBoundary(void);

    /// Deallocate PETSc and local data structures.
    void deallocate(void);

    /** Set mesh label associated with boundary condition surface.
     *
     * @param[in] value Label of surface (from mesh generator).
     */
    void label(const char* value);

    /** Get mesh label associated with boundary condition surface.
     *
     * @returns Label of surface (from mesh generator).
     */
    const char* label(void) const;

    /** Set name of field in solution to constrain.
     *
     * @param[in] value Name of field in solution to constrain.
     */
    void field(const char* value);

    /** Get name of field in solution to constrain.
     *
     * @returns Name of field in solution to constrain.
     */
    const char* field(void) const;

    /** Get mesh associated with integrator domain.
     *
     * @returns Mesh associated with integrator domain.
     */
    const pylith::topology::Mesh& domainMesh(void) const;

    /** Verify configuration is acceptable.
     *
     * @param[in] solution Solution field.
     */
    virtual
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

    // These will become methods in ConstraintPhysics.

    /** Setup auxiliary subfields (discretization and query fns).
     *
     * Create subfields in auxiliary fields (includes name of the field,
     * vector field type, discretization, and scale for
     * nondimensionalization) and set query functions for filling them
     * from a spatial database.
     *
     * @attention The order of the calls to subfieldAdd() must match the
     * order of the auxiliary subfields in the FE kernels.
     *
     * @param[in] solution Solution field.
     */
    virtual
    void _auxFieldSetup(const pylith::topology::Field& solution) = 0;

    /** Set kernel for computing value of constrained degree of freedom from auxiliary field.
     *
     * @param solution Solution field.
     */
    virtual
    void _setFEKernelConstraint(const pylith::topology::Field& solution) = 0;

    /** Get point-wise function (kernel) for settings constraint from auxiliary field.
     *
     * @returns Point-wise function.
     */
    virtual
    PetscPointFunc _getFEKernelConstraint(void) = 0;

    // PROTECTED MEMBERS //////////////////////////////////////////////////
protected:

    pylith::topology::Mesh* _boundaryMesh; ///< Boundary mesh.
    std::string _label; ///< Label to identify boundary condition points in mesh.
    std::string _field; ///< Name of solution field for boundary condition.
    pylith::topology::FieldBase::Description _description; ///< Description for constrained field.

    // NOT IMPLEMENTED ////////////////////////////////////////////////////
private:

    ConstraintBoundary(const ConstraintBoundary&); ///< Not implemented.
    const ConstraintBoundary& operator=(const ConstraintBoundary&); ///< Not implemented.

}; // class ConstraintBoundary

#endif // pylith_bc_constraintboundary_hh

// End of file
