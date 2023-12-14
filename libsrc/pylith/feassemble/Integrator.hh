// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2023, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================
#pragma once

#include "pylith/feassemble/PhysicsImplementation.hh" // ISA PhysicsImplementation

#include "pylith/problems/problemsfwd.hh" // HASA Physics
#include "pylith/topology/topologyfwd.hh" // USES Field

#include "pylith/utils/petscfwd.h" // USES PetscMat, PetscVec
#include "pylith/utils/utilsfwd.hh" // HOLDSA Logger

class pylith::feassemble::Integrator : public pylith::feassemble::PhysicsImplementation {
    friend class TestIntegrator; // unit testing

    // PUBLIC ENUM ////////////////////////////////////////////////////////////////////////////////
public:

    enum EquationPart {
        LHS=0,
        RHS=1,
        LHS_LUMPED_INV=2,
        LHS_WEIGHTED=3,
    };

    enum NewJacobianTriggers {
        NEW_JACOBIAN_NEVER=0x0, // Never needs new Jacobian.
        NEW_JACOBIAN_ALWAYS=0x1, // Always needs new Jacobian.
        NEW_JACOBIAN_TIME_STEP_CHANGE=0x2, // Needs new Jacobian if time step changes.
        NEW_JACOBIAN_UPDATE_STATE_VARS=0x4, // Needs new Jacobian after updating state variables.
    };

    // PUBLIC STRUCTS /////////////////////////////////////////////////////////////////////////////
public:

    /// Project kernels (pointwise functions) for updating state variables or computing derived fields.
    struct ProjectKernels {
        std::string subfield; ///< Name of subfield for function.
        PetscPointFunc f; ///< Point-wise function.

        ProjectKernels(void) :
            subfield(""),
            f(NULL) {}


        ProjectKernels(const char* subfieldValue,
                       PetscPointFunc fValue) :
            subfield(subfieldValue),
            f(fValue) {}


    }; // ProjectKernels

    // PUBLIC MEMBERS /////////////////////////////////////////////////////////////////////////////
public:

    /** Constructor
     *
     * @param[in] physics Physics implemented by integrator.
     */
    Integrator(pylith::problems::Physics* const physics);

    /// Destructor
    virtual ~Integrator(void);

    /** Set name of label used to identify integration domain.
     *
     * @param name Name of label.
     */
    void setLabelName(const char* name);

    /** Get name of label used to identify integration domain.
     *
     * @returns Name of label.
     */
    const char* getLabelName(void) const;

    /** Set value of label used to identify integration domain.
     *
     * @param value Label value.
     */
    void setLabelValue(const int value);

    /** Get value of label used to identify integration domain.
     *
     * @returns Label value
     */
    int getLabelValue(void) const;

    /** Check whether LHS Jacobian needs to be recomputed.
     *
     * @param[in] dtChanged True if time step has changed since previous Jacobian computation.
     * @returns True if Jacobian needs to be recomputed, false otherwise.
     */
    bool needNewLHSJacobian(const bool dtChanged);

    /** Check whether LHS lumped Jacobian needs to be recomputed.
     *
     * @param[in] dtChanged True if time step has changed since previous Jacobian computation.
     * @returns True if lumped Jacobian needs to be recomputed, false otherwise.
     */
    bool needNewLHSJacobianLumped(const bool dtChanged);

    /** Set LHS Jacobian trigger.
     *
     * @param[in] value Triggers for needing new LHS Jacobian.
     */
    void setLHSJacobianTriggers(const int value);

    /** Set LHS lumped Jacobian trigger.
     *
     * @param[in] value Triggers for needing new LHS lumped Jacobian.
     */
    void setLHSJacobianLumpedTriggers(const int value);

    /** Initialize integration domain, auxiliary field, and derived field. Update observers.
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

    /** Compute RHS residual for G(t,s).
     *
     * @param[out] residual Field for residual.
     * @param[in] integrationData Data needed to integrate governing equations.
     */
    virtual
    void computeRHSResidual(pylith::topology::Field* residual,
                            const pylith::feassemble::IntegrationData& integrationData) = 0;

    /** Compute LHS residual for F(t,s,\dot{s}).
     *
     * @param[out] residual Field for residual.
     * @param[in] integrationData Data needed to integrate governing equations.
     */
    virtual
    void computeLHSResidual(pylith::topology::Field* residual,
                            const pylith::feassemble::IntegrationData& integrationData) = 0;

    /** Compute LHS Jacobian and preconditioner for F(t,s,\dot{s}) with implicit time-stepping.
     *
     * @param[out] jacobianMat PETSc Mat with Jacobian sparse matrix.
     * @param[out] precondMat PETSc Mat with Jacobian preconditioning sparse matrix.
     * @param[in] integrationData Data needed to integrate governing equations.
     */
    virtual
    void computeLHSJacobian(PetscMat jacobianMat,
                            PetscMat precondMat,
                            const pylith::feassemble::IntegrationData& integrationData) = 0;

    /** Compute inverse of lumped LHS Jacobian for F(t,s,\dot{s}) with explicit time-stepping.
     *
     * @param[out] jacobianInv Inverse of lumped Jacobian as a field.
     * @param[in] integrationData Data needed to integrate governing equations.
     */
    virtual
    void computeLHSJacobianLumpedInv(pylith::topology::Field* jacobianInv,
                                     const pylith::feassemble::IntegrationData& integrationData) = 0;

    // PROTECTED METHODS //////////////////////////////////////////////////////////////////////////
protected:

    /** Set constants used in finite-element kernels.
     *
     * @param[in] solution Solution field.
     * @param[in] dt Current time step.
     */
    virtual
    void _setKernelConstants(const pylith::topology::Field& solution,
                             const PylithReal dt) const;

    /** Update state variables as needed.
     *
     * @param[in] t Current time.
     * @param[in] dt Current time step.
     * @param[in] solution Field with current trial solution.
     */
    virtual
    void _updateStateVars(const PylithReal t,
                          const PylithReal dt,
                          const pylith::topology::Field& solution);

    /// Compute diagnostic field from auxiliary field.
    virtual
    void _computeDiagnosticField(void);

    /** Compute fields derived from solution and auxiliary field.
     *
     * @param[in] t Current time.
     * @param[in] dt Current time step.
     * @param[in] solution Field with current trial solution.
     */
    virtual
    void _computeDerivedField(const PylithReal t,
                              const PylithReal dt,
                              const pylith::topology::Field& solution);

    // PROTECTED MEMBERS //////////////////////////////////////////////////////////////////////////
protected:

    std::string _labelName; ///< Name of label associated with integration domain.
    int _labelValue; ///< Value of label associated with integration domain.

    int _lhsJacobianTriggers; // Triggers for needing new LHS Jacobian.
    int _lhsJacobianLumpedTriggers; // Triggers for needing new LHS lumped Jacobian.

    /// True if we have kernels for operation, false otherwise.
    bool _hasRHSResidual;
    bool _hasLHSResidual;
    bool _hasLHSJacobian;
    bool _hasLHSJacobianLumped;

    /// True if we need to recompute Jacobian for operator, false otherwise.
    /// Default is false;
    bool _needNewLHSJacobian;
    bool _needNewLHSJacobianLumped;

    // NOT IMPLEMENTED ////////////////////////////////////////////////////////////////////////////
private:

    Integrator(void); /// Not implemented.
    Integrator(const Integrator&); ///< Not implemented.
    const Integrator& operator=(const Integrator&); ///< Not implemented.

}; // Integrator

// End of file
