// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================
#pragma once

#include "pylith/feassemble/PhysicsImplementation.hh" // ISA PhysicsImplementation

#include "pylith/problems/problemsfwd.hh" // HASA Physics
#include "pylith/topology/topologyfwd.hh" // USES Field

#include "pylith/utils/array.hh" // HASA int_array
#include "pylith/utils/types.hh" // HASA PetscUserFieldFunc
#include "pylith/utils/utilsfwd.hh" // HOLDSA Logger

class pylith::feassemble::Constraint : public pylith::feassemble::PhysicsImplementation {
    friend class TestConstraint; // unit testing

    // PUBLIC STRUCTS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Project kernels (pointwise functions) for computing diagnostic fields.
    struct ProjectKernels {
        std::string subfield; ///< Name of subfield for function.
        PetscBdPointFunc f; ///< Point-wise function.

        ProjectKernels(void) :
            subfield(""),
            f(NULL) {}


        ProjectKernels(const char* subfieldValue,
                       PetscBdPointFunc fValue) :
            subfield(subfieldValue),
            f(fValue) {}


    }; // ProjectKernels

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

    /** Set kernels for computing diagnostic field.
     *
     * @param kernels Array of kernels for computing diagnostic field.
     */
    void setKernelsDiagnosticField(const std::vector<ProjectKernels>& kernels);

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
     * @param[in] notification Type of notification.
     */
    virtual
    void poststep(const PylithReal t,
                  const PylithInt tindex,
                  const PylithReal dt,
                  const pylith::topology::Field& solution,
                  const pylith::problems::Observer::NotificationType notification);

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

    // PROTECTED METHODS //////////////////////////////////////////////////////////////////////////
protected:

    /** Set constants used in finite-element kernels.
     *
     * @param[in] solution Solution field.
     * @param[in] dt Current time step.
     */
    void _setKernelConstants(const pylith::topology::Field& solution,
                             const PylithReal dt) const;

    /// Compute diagnostic field from auxiliary field.
    virtual
    void _computeDiagnosticField(void);

    // PROTECTED MEMBERS ///////////////////////////////////////////////////////////////////////////////////////////////
protected:

    std::string _subfieldName; ///< Name of solution subfield that is constrained.
    std::string _labelName; ///< Name of label associated with integration domain.
    int _labelValue; ///< Value of label associated with integration domain.

    int_array _constrainedDOF; ///< List of constrained degrees of freedom at each location.
    pylith::topology::Mesh* _boundaryMesh; ///< Boundary mesh.
    PylithReal _tSolution; ///< Time used for current solution.

    std::vector<ProjectKernels> _kernelsDiagnosticField; ///< kernels for computing diagnostic field.

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    Constraint(void); /// Not implemented.
    Constraint(const Constraint &); ///< Not implemented
    const Constraint& operator=(const Constraint&); ///< Not implemented

}; // class Constraint

// End of file
