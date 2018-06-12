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
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

/**
 * @file libsrc/feassemble/IntegratorPointwise.hh
 *
 * @brief Object containing operations for implicit and explicit
 * time integration of the equations defined by pointwise functions.
 */

#if !defined(pylith_feassemble_integratorpointwise_hh)
#define pylith_feassemble_integratorpointwise_hh

// Include directives ---------------------------------------------------
#include "feassemblefwd.hh" // forward declarations

#include "pylith/utils/PyreComponent.hh" // ISA PyreComponent
#include "pylith/feassemble/ObservedSubject.hh" // ISA ObservedSubject

#include "pylith/topology/FieldBase.hh" // USES FieldBase
#include "pylith/utils/petscfwd.h" // USES PetscMat, PetscVec
#include "pylith/utils/utilsfwd.hh" // HOLDSA Logger

#include "spatialdata/spatialdb/spatialdbfwd.hh" // HASA spatialdb
#include "spatialdata/units/unitsfwd.hh" // HASA Nondimensional

// IntegratorPointwise -------------------------------------------------
/** @brief General operations for implicit and explicit
 * time integration of equations defined by pointwise functions.
 */
class pylith::feassemble::IntegratorPointwise :
    public pylith::utils::PyreComponent,
    public pylith::feassemble::ObservedSubject {
    friend class TestIntegratorPointwise;   // unit testing

    // PUBLIC TYPEDEFS //////////////////////////////////////////////////////
public:

    // PUBLIC MEMBERS ///////////////////////////////////////////////////////
public:

    /// Constructor
    IntegratorPointwise(void);

    /// Destructor
    virtual ~IntegratorPointwise(void);

    /// Deallocate PETSc and local data structures.
    virtual
    void deallocate(void);

    /** Get auxiliary field.
     *
     * @return field Field over integrator domain.
     */
    const pylith::topology::Field* auxField(void) const;

    /** Get derived field.
     *
     * @return field Field over integrator domain.
     */
    const pylith::topology::Field* derivedField(void) const;

    /** Set spatial database for filling auxiliary subfields.
     *
     * @param[in] value Pointer to database.
     */
    void auxFieldDB(spatialdata::spatialdb::SpatialDB* value);

    /** Set discretization information for auxiliary subfield.
     *
     * @param[in] name Name of auxiliary subfield.
     * @param[in] basisOrder Polynomial order for basis.
     * @param[in] quadOrder Order of quadrature rule.
     * @param[in] isBasisContinuous True if basis is continuous.
     * @param[in] feSpace Finite-element space.
     */
    void auxSubfieldDiscretization(const char* name,
                                   const int basisOrder,
                                   const int quadOrder,
                                   const bool isBasisContinuous,
                                   const pylith::topology::FieldBase::SpaceEnum feSpace);

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

    /** Set manager of scales used to nondimensionalize problem.
     *
     * @param dim Nondimensionalizer.
     */
    void normalizer(const spatialdata::units::Nondimensional& dim);

    /** Set gravity field.
     *
     * @param g Gravity field.
     */
    void gravityField(spatialdata::spatialdb::GravityField* const g);

    /** Verify configuration is acceptable.
     *
     * @param[in] solution Solution field.
     */
    virtual
    void verifyConfiguration(const pylith::topology::Field& solution) const = 0;

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
                  const PylithReal dt,
                  const PylithInt tindex,
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
     * @param[in] tshift Scale for time derivative.
     * @param[in] solution Field with current trial solution.
     * @param[in] solutionDot Field with time derivative of current trial solution.
     */
    virtual
    void computeLHSJacobianImplicit(PetscMat jacobianMat,
                                    PetscMat precondMat,
                                    const PylithReal t,
                                    const PylithReal dt,
                                    const PylithReal tshift,
                                    const pylith::topology::Field& solution,
                                    const pylith::topology::Field& solutionDot) = 0;


    /** Compute inverse of lumped LHS Jacobian for F(t,s,\dot{s}) with explicit time-stepping.
     *
     * @param[out] jacobianInv Inverse of lumped Jacobian as a field.
     * @param[in] t Current time.
     * @param[in] dt Current time step.
     * @param[in] tshift Scale for time derivative.
     * @param[in] solution Field with current trial solution.
     */
    virtual
    void computeLHSJacobianLumpedInv(pylith::topology::Field* jacobianInv,
                                     const PylithReal t,
                                     const PylithReal dt,
                                     const PylithReal tshift,
                                     const pylith::topology::Field& solution) = 0;


    // PROTECTED METHODS //////////////////////////////////////////////////
protected:

    /** Get factory for setting up auxliary fields.
     *
     * @returns Factor for auxiliary fields.
     */
    virtual
    pylith::feassemble::AuxiliaryFactory* _auxFactory(void) = 0;

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

    spatialdata::units::Nondimensional* _normalizer;   ///< Nondimensionalizer.
    spatialdata::spatialdb::GravityField* _gravityField; ///< Gravity field.

    pylith::topology::Field* _auxField; ///< Auxiliary field for this integrator.
    pylith::topology::Field* _derivedField; ///< Derived field for this integrator.

    pylith::utils::EventLogger* _logger;   ///< Event logger.

    /// True if we need to recompute Jacobian for operator, false otherwise.
    /// Default is false;
    bool _needNewRHSJacobian;
    bool _needNewLHSJacobian;

    typedef std::map<std::string, PetscPointFunc> UpdateStateVarsMap;
    UpdateStateVarsMap _updateStateVarsKernels;


    // NOT IMPLEMENTED //////////////////////////////////////////////////////
private:

    IntegratorPointwise(const IntegratorPointwise&); ///< Not implemented.
    const IntegratorPointwise& operator=(const IntegratorPointwise&); ///< Not implemented.

}; // IntegratorPointwise

#endif // pylith_feassemble_integratorpointwise_hh


// End of file
