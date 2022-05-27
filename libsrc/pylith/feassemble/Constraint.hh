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

/** @file libsrc/feassemble/Constraint.hh
 *
 * @brief C++ abstract base class for constraining degrees of freedom in the solution.
 */

#if !defined(pylith_feassemble_constraint_hh)
#define pylith_feassemble_constraint_hh

#include "pylith/feassemble/PhysicsImplementation.hh" // ISA PhysicsImplementation

#include "pylith/problems/problemsfwd.hh" // HASA Physics
#include "pylith/topology/topologyfwd.hh" // USES Field

#include "pylith/utils/array.hh" // HASA int_array
#include "pylith/utils/types.hh" // HASA PetscUserFieldFunc
#include "pylith/utils/utilsfwd.hh" // HOLDSA Logger

class pylith::feassemble::Constraint : public pylith::feassemble::PhysicsImplementation {
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

    /** Set name of label marking boundary associated with constraint.
     *
     * @param[in] value Name of label for surface (from mesh generator).
     */
    void setLabelName(const char* value);

    /** Get name of label marking boundary associated with constraint.
     *
     * @returns Name of label for surface (from mesh generator).
     */
    const char* getLabelName(void) const;

    /** Set value of label marking boundary associated with constraint.
     *
     * @param[in] value Value of label for surface (from mesh generator).
     */
    void setLabelValue(const int value);

    /** Get value of label marking boundary associated with constraint.
     *
     * @returns Value of label for surface (from mesh generator).
     */
    int getLabelValue(void) const;

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

    /** Get mesh associated with constrained boundary.
     *
     * @returns Mesh associated with constrained boundary.
     */
    const pylith::topology::Mesh& getPhysicsDomainMesh(void) const;

    /** Initialize constraint.
     *
     * @param[in] solution Solution field (layout).
     */
    virtual
    void initialize(const pylith::topology::Field& solution);

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

    /** Set auxiliary field values for current time.
     *
     * @param[in] t Current time.
     */
    virtual
    void setState(const PylithReal t);

    /** Set constrained values in solution field.
     *
     * @param[inout] integrationData Data needed to integrate governing equation.
     */
    virtual
    void setSolution(pylith::feassemble::IntegrationData* integrationData) = 0;

    // PROTECTED MEMBERS ///////////////////////////////////////////////////////////////////////////////////////////////
protected:

    std::string _subfieldName; ///< Name of solution subfield that is constrained.
    std::string _labelName; ///< Name of label associated with integration domain.
    int _labelValue; ///< Value of label associated with integration domain.

    int_array _constrainedDOF; ///< List of constrained degrees of freedom at each location.
    pylith::topology::Mesh* _boundaryMesh; ///< Boundary mesh.
    PylithReal _tSolution; ///< Time used for current solution.

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    Constraint(void); /// Not implemented.
    Constraint(const Constraint &); ///< Not implemented
    const Constraint& operator=(const Constraint&); ///< Not implemented

}; // class Constraint

#endif // pylith_feassemble_constraint_hh

// End of file
