// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, Rice University
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2018 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

/**
 * @file libsrc/feassemble/Integrator.hh
 *
 * @brief Abstract base class for finite-element integration.
 */

#if !defined(pylith_feassemble_integrator_hh)
#define pylith_feassemble_integrator_hh

#include "feassemblefwd.hh" // forward declarations

#include "pylith/utils/GenericComponent.hh" // ISA GenericComponent

#include "pylith/problems/problemsfwd.hh" // HASA Physics

#include "pylith/topology/FieldBase.hh" // USES FieldBase
#include "pylith/utils/petscfwd.h" // USES PetscMat, PetscVec
#include "pylith/utils/utilsfwd.hh" // HOLDSA Logger

class pylith::feassemble::Integrator : public pylith::utils::GenericComponent {
    friend class TestIntegrator; // unit testing

    // PUBLIC STRUCTS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Kernels (point-wise functions) for residual.
    struct ResidualKernels {
        std::string subfield; ///< Name of subfield
        PetscPointFunc r0; ///< f0 (RHS) or g0 (LHS) function.
        PetscPointFunc r1; ///< f1 (RHS) or g1 (LHS) function.
    }; // ResidualKernels

    /// Kernels (point-wise functions) for Jacobian;
    struct JacobianKernels {
        std::string subfieldTrial; ///< Name of subfield associated with trial function (row in Jacobian).
        std::string subfieldBasis; ///< Name of subfield associated with basis function (column in Jacobian).
        PetscPointJac j0; ///< J0 function.
        PetscPointJac j1; ///< J1 function.
        PetscPointJac j2; ///< J2 function.
        PetscPointJac j3; ///< J3 function.
    }; // JacobianKernels

    /// Project kernels (point-wise functions) for updating state variables or computing derived fields.
    struct ProjectKernels {
        std::string subfield; ///< Name of subfield for function.
        PetscPointFunc f; ///< Point-wise function.
    }; // ProjectKernels

    // PUBLIC MEMBERS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Constructor
    Integrator(pylith::problems::Physics* const physics);

    /// Destructor
    virtual ~Integrator(void);

    /// Deallocate PETSc and local data structures.
    virtual
    void deallocate(void);

    /** Get mesh associated with integrator domain.
     *
     * @returns Mesh associated with integrator domain.
     */
    virtual
    const pylith::topology::Mesh& getIntegationDomainMesh(void) const = 0;

    /** Get auxiliary field.
     *
     * @return field Field over integrator domain.
     */
    const pylith::topology::Field* getAuxiliaryField(void) const;

    /** Get derived field.
     *
     * @return field Field over integrator domain.
     */
    const pylith::topology::Field* getDerivedField(void) const;

    /** Check whether RHS Jacobian needs to be recomputed.
     *
     * @returns True if Jacobian needs to be recomputed, false otherwise.
     */
    bool needNewRHSJacobian(void) const;

    /** Check whether LHS Jacobian needs to be recomputed.
     *
     * @returns True if Jacobian needs to be recomputed, false otherwise.
     */
    bool needNewLHSJacobian(void) const;

    /** Initialize integrator.
     *
     * @param[in] solution Solution field (layout).
     */
    virtual
    void initialize(const pylith::topology::Field& solution);

    /** Update at beginning of time step.
     *
     * @param[in] t Current time.
     * @param[in] dt Current time step.
     */
    virtual
    void prestep(const PylithReal t,
                 const PylithReal dt);

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

    /** Compute RHS residual for G(t,s).
     *
     * @param[out] residual Field for residual.
     * @param[in] t Current time.
     * @param[in] dt Current time step.
     * @param[in] solution Field with current trial solution.
     */
    virtual
    void computeRHSResidual(pylith::topology::Field* residual,
                            const PylithReal t,
                            const PylithReal dt,
                            const pylith::topology::Field& solution) = 0;

    /** Compute RHS Jacobian and preconditioner for G(t,s).
     *
     * @param[out] jacobianMat PETSc Mat with Jacobian sparse matrix.
     * @param[out] precondMat PETSc Mat with Jacobian preconditioning sparse matrix.
     * @param[in] t Current time.
     * @param[in] dt Current time step.
     * @param[in] solution Field with current trial solution.
     */
    virtual
    void computeRHSJacobian(PetscMat jacobianMat,
                            PetscMat preconMat,
                            const PylithReal t,
                            const PylithReal dt,
                            const pylith::topology::Field& solution) = 0;

    /** Compute LHS residual for F(t,s,\dot{s}).
     *
     * @param[out] residual Field for residual.
     * @param[in] t Current time.
     * @param[in] dt Current time step.
     * @param[in] solution Field with current trial solution.
     * @param[in] solutionDot Field with time derivative of current trial solution.
     */
    virtual
    void computeLHSResidual(pylith::topology::Field* residual,
                            const PylithReal t,
                            const PylithReal dt,
                            const pylith::topology::Field& solution,
                            const pylith::topology::Field& solutionDot) = 0;

    /** Compute LHS Jacobian and preconditioner for F(t,s,\dot{s}) with implicit time-stepping.
     *
     * @param[out] jacobianMat PETSc Mat with Jacobian sparse matrix.
     * @param[out] precondMat PETSc Mat with Jacobian preconditioning sparse matrix.
     * @param[in] t Current time.
     * @param[in] dt Current time step.
     * @param[in] s_tshift Scale for time derivative.
     * @param[in] solution Field with current trial solution.
     * @param[in] solutionDot Field with time derivative of current trial solution.
     */
    virtual
    void computeLHSJacobian(PetscMat jacobianMat,
                            PetscMat precondMat,
                            const PylithReal t,
                            const PylithReal dt,
                            const PylithReal s_tshift,
                            const pylith::topology::Field& solution,
                            const pylith::topology::Field& solutionDot) = 0;

    /** Compute inverse of lumped LHS Jacobian for F(t,s,\dot{s}) with explicit time-stepping.
     *
     * @param[out] jacobianInv Inverse of lumped Jacobian as a field.
     * @param[in] t Current time.
     * @param[in] dt Current time step.
     * @param[in] s_tshift Scale for time derivative.
     * @param[in] solution Field with current trial solution.
     */
    virtual
    void computeLHSJacobianLumpedInv(pylith::topology::Field* jacobianInv,
                                     const PylithReal t,
                                     const PylithReal dt,
                                     const PylithReal s_tshift,
                                     const pylith::topology::Field& solution) = 0;

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

    /** Update state variables as needed.
     *
     * @param[in] t Current time.
     * @param[in] dt Current time step.
     * @param[in] solution Field with current trial solution.
     */
    virtual
    void _updateStateVars(const PylithReal t,
                          const PylithReal dt,
                          const pylith::topology::Field& solution) = 0;

    /** Compute fields derived from solution and auxiliary field.
     *
     * @param[in] t Current time.
     * @param[in] dt Current time step.
     * @param[in] solution Field with current trial solution.
     */
    virtual
    void _computeDerivedFields(const PylithReal t,
                               const PylithReal dt,
                               const pylith::topology::Field& solution) = 0;

    // PROTECTED MEMBERS ///////////////////////////////////////////////////////////////////////////////////////////////
protected:

    pylith::problems::Physics* const _physics; ///< Physics associated with integrator.
    pylith::topology::Field* _auxField; ///< Auxiliary field for this integrator.
    pylith::topology::Field* _derivedField; ///< Derived field for this integrator.
    pylith::feassemble::Observers* _observers; ///< Observers component.

    pylith::utils::EventLogger* _logger; ///< Event logger.

    /// True if we need to recompute Jacobian for operator, false otherwise.
    /// Default is false;
    bool _needNewRHSJacobian;
    bool _needNewLHSJacobian;

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    Integrator(const Integrator&); ///< Not implemented.
    const Integrator& operator=(const Integrator&); ///< Not implemented.

}; // Integrator

#endif // pylith_feassemble_integrator_hh

// End of file
