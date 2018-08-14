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

// Include directives ---------------------------------------------------
#include "feassemblefwd.hh" // forward declarations

#include "pylith/utils/GenericComponent.hh" // ISA GenericComponent

#include "pylith/topology/FieldBase.hh" // USES FieldBase
#include "pylith/utils/petscfwd.h" // USES PetscMat, PetscVec
#include "pylith/utils/utilsfwd.hh" // HOLDSA Logger

class pylith::feassemble::Integrator : public pylith::utils::GenericComponent {
    friend class TestIntegrator; // unit testing

    // PUBLIC STRUCTS ////////////////////////////////////////////////////////
public:

    struct ResidualKernels {
        std::string field;
        PetscPointFunc r0;
        PetscPointFunc r1;
    }; // ResidualKernels

    struct JacobianKernels {
        std::string fieldTrial;
        std::string fieldBasis;
        PetscJacobianFunc j0;
        PetscJacobianFunc j1;
        PetscJacobianFunc j2;
        PetscJacobianFunc j3;
    }

    struct ProjectKernels {
        std::string field;
        PetscPointFunc p;
    }

    // PUBLIC MEMBERS ///////////////////////////////////////////////////////
public:

    /// Constructor
    Integrator(const pylith::problems::Physics* physics);

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
    const pylith::topology::Mesh& getIntegationDomainMesh(void) const;

    /** Get auxiliary field.
     *
     * @return field Field over integrator domain.
     */
    const pylith::topology::Field* getAuxField(void) const;

    /** Get derived field.
     *
     * @return field Field over integrator domain.
     */
    const pylith::topology::Field* getDerivedField(void) const;

    /** Register observer to receive notifications.
     *
     * Observers are used for output.
     *
     * @param[in] observer Observer to receive notifications.
     */
    void registerObserver(pylith::feassemble::Observer* observer);

    /** Remove observer from receiving notifications.
     *
     * @param[in] observer Observer to remove.
     */
    void removeObserver(pylith::feassemble::Observer* observer);

    /** Check whether RHS Jacobian needs to be recomputed.
     *
     * @returns True if Jacobian needs to be recomputed, false otherwise.
     */
    virtual
    bool needNewRHSJacobian(void) const;

    /** Check whether LHS Jacobian needs to be recomputed.
     *
     * @returns True if Jacobian needs to be recomputed, false otherwise.
     */
    virtual
    bool needNewLHSJacobian(void) const;

    /** Initialize integrator.
     *
     * @param[in] solution Solution field (layout).
     */
    virtual
    void initialize(const pylith::topology::Field& solution) = 0;

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

    // PROTECTED METHODS //////////////////////////////////////////////////
protected:

    /** Set constants used in finite-element integrations.
     *
     * @param[in] solution Solution field.
     * @param[in] dt Current time step.
     */
    virtual
    void _setFEConstants(const pylith::topology::Field& solution,
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

    /** Compute fields derived from solution and auxiliary field.
     *
     * @param[in] t Current time.
     * @param[in] dt Current time step.
     * @param[in] solution Field with current trial solution.
     */
    virtual
    void _computeDerivedFields(const PylithReal t,
                               const PylithReal dt,
                               const pylith::topology::Field& solution);

    // PROTECTED MEMBERS ////////////////////////////////////////////////////
protected:

    const pylith::problems::Physics* _physics; ///< Physics associated with integrator.
    pylith::topology::Field* _auxField; ///< Auxiliary field for this integrator.
    pylith::topology::Field* _derivedField; ///< Derived field for this integrator.
    pylith::feassemble::ObservedComponent* _observed; ///< Observed component.

    pylith::utils::EventLogger* _logger; ///< Event logger.

    /// True if we need to recompute Jacobian for operator, false otherwise.
    /// Default is false;
    bool _needNewRHSJacobian;
    bool _needNewLHSJacobian;

    typedef std::map<std::string, PetscPointFunc> UpdateStateVarsMap;
    UpdateStateVarsMap _updateStateVarsKernels;

    // NOT IMPLEMENTED //////////////////////////////////////////////////////
private:

    Integrator(const Integrator&); ///< Not implemented.
    const Integrator& operator=(const Integrator&); ///< Not implemented.

}; // Integrator

#endif // pylith_feassemble_integrator_hh

// End of file
