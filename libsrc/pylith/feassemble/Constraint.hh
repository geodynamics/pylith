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

/** @file libsrc/feassemble/Constraint.hh
 *
 * @brief C++ abstract base class defining interface for constraining degrees of freedom in the solution.
 */

#if !defined(pylith_feassemble_constraint_hh)
#define pylith_feassemble_constraint_hh

#include "pylith/feassemble/feassemblefwd.hh"

#include "pylith/utils/GenericComponent.hh" // ISA GenericComponent

#include "pylith/problems/problemsfwd.hh" // HASA Physics
#include "pylith/topology/topologyfwd.hh" // USES Field

#include "pylith/utils/array.hh" // HASA int_array
#include "pylith/utils/utilsfwd.hh" // HOLDSA Logger

class pylith::feassemble::Constraint : public pylith::utils::GenericComponent {
    friend class TestConstraint; // unit testing

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /** Constructor
     *
     * @param[in] physics Physics implemented by constraint.
     */
    Constraint(pylith::problems::Physics* const physics);

    /// Destructor.
    virtual ~Constraint(void);

    /// Deallocate PETSc and local data structures.
    virtual
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

    /** Set label marking constrained degrees of freedom.
     *
     * @param[in] value Label of constrained degrees of freedom (from mesh generator).
     */
    void setMarkerLabel(const char* value);

    /** Get label marking constrained degrees of freedom.
     *
     * @returns Label of constrained degrees of freedom (from mesh generator).
     */
    const char* getMarkerLabel(void) const;

    /** Set name of constrained solution subfield.
     *
     * @param[in] value Name of solution subfield.
     */
    void setSubfieldName(const char* value);

    /** Get name of constrained solution subfield.
     *
     * @preturn Name of solution subfield.
     */
    const char* getSubfieldName(void) const;

    /** Get mesh associated with constrained domain.
     *
     * @returns Mesh associated with constrained domain.
     */
    virtual
    const pylith::topology::Mesh& getConstraintDomainMesh(void) const = 0;

    /** Get auxiliary field.
     *
     * @returns field Field over boundary.
     */
    const pylith::topology::Field* getAuxiliaryField(void) const;

    /** Get derived field.
     *
     * @return field Field over integrator domain.
     */
    const pylith::topology::Field* getDerivedField(void) const;

    /** Set constraint kernel.
     *
     * @param kernel Kernel to compute constrained value from auxiliary field.
     */
    void setKernelConstraint(const PetscPointFunc kernel);

    /** Initialize integrator.
     *
     * @param[in] solution Solution field (layout).
     */
    virtual
    void initialize(const pylith::topology::Field& solution);

    /** Update auxiliary field at beginning of time step.
     *
     * @param[in] t Current time.
     * @param[in] dt Current time step.
     */
    virtual
    void prestep(const double t,
                 const double dt);

    /** Update at end of time step.
     *
     * @param[in] t Current time.
     * @param[in] tindex Current time step.
     * @param[in] dt Current time step.
     * @param[in] solution Solution at time t.
     */
    virtual
    void poststep(const PylithReal t,
                  const PylithInt tindex,
                  const PylithReal dt,
                  const pylith::topology::Field& solution);

    /** Set constrained values in solution field.
     *
     * @param[out] solution Solution field.
     * @param[in] t Current time.
     */
    virtual
    void setSolution(pylith::topology::Field* solution,
                     const double t);

    // PROTECTED METHODS ///////////////////////////////////////////////////////////////////////////////////////////////
protected:

    /** Set constants used in finite-element kernels.
     *
     * @param[in] solution Solution field.
     * @param[in] dt Current time step.
     */
    virtual
    void _setKernelConstants(const pylith::topology::Field& solution,
                             const PylithReal dt) const;

    // PROTECTED MEMBERS ///////////////////////////////////////////////////////////////////////////////////////////////
protected:

    int_array _constrainedDOF; ///< List of constrained degrees of freedom at each location.
    std::string _constraintLabel; ///< Label marking constrained degrees of freedom.
    std::string _subfieldName; ///< Name of solution subfield that is constrained.
    PetscPointFunc _kernelConstraint; ///< Kernel for computing constrained values from auxiliary field.
    pylith::problems::Physics* const _physics; ///< Physics associated with constraint.
    pylith::topology::Field* _auxiliaryField; ///< Auxiliary field for this constraint.
    pylith::topology::Field* _derivedField; ///< Derived field for this constraint.
    pylith::feassemble::Observers* _observers; ///< Observers component.

    pylith::utils::EventLogger* _logger; ///< Event logger.

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    Constraint(void); /// Not implemented.
    Constraint(const Constraint &m); ///< Not implemented
    const Constraint& operator=(const Constraint& m); ///< Not implemented

}; // class Constraint

#endif // pylith_feassemble_constraint_hh

// End of file
