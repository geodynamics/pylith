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

#include "pylith/topology/FieldBase.hh" // USES FieldBase

#include "pylith/topology/topologyfwd.hh" // HOLDSA Field
#include "pylith/utils/petscfwd.h" // USES PetscMat, PetscVec

#include "spatialdata/spatialdb/spatialdbfwd.hh" // HASA spatialdb
#include "spatialdata/units/unitsfwd.hh" // HASA Nondimensional

#include <map> // HOLDSA std::map

// IntegratorPointwise -------------------------------------------------
/** @brief General operations for implicit and explicit
 * time integration of equations defined by pointwise functions.
 */
class pylith::feassemble::IntegratorPointwise : public pylith::utils::PyreComponent {
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

    /** Get auxiliary fields.
     *
     * @return field Field over material.
     */
    const pylith::topology::Field& auxFields(void) const;

    /** Check whether material has a given auxiliary field.
     *
     * @param[in] name Name of field.
     *
     * @returns True if material has auxiliary field, false otherwise.
     */
    bool hasAuxField(const char* name);

    /** Get auxiliary field.
     *
     * @param[out] field Field over material.
     * @param[in] name Name of field to retrieve.
     */
    void getAuxField(pylith::topology::Field *field,
                     const char* name) const;

    /** Set spatial database for auxiliary fields.
     *
     * @param[in] value Pointer to database.
     */
    void auxFieldsDB(spatialdata::spatialdb::SpatialDB* value);

    /** Set discretization information for auxiliary subfield.
     *
     * @param[in] name Name of auxiliary subfield.
     * @param[in] basisOrder Polynomial order for basis.
     * @param[in] quadOrder Order of quadrature rule.
     * @param[in] isBasisContinuous True if basis is continuous.
     * @param[in] feSpace Finite-element space.
     */
    void auxFieldDiscretization(const char* name,
                                const int basisOrder,
                                const int quadOrder,
                                const bool isBasisContinuous,
                                const pylith::topology::FieldBase::SpaceEnum feSpace);

    /** Get discretization information for auxiliary subfield.
     *
     * @param[in] name Name of subfield.
     * @return Discretization information for auxiliary subfield. If
     * discretization information was not set, then use "default".
     */
    const pylith::topology::FieldBase::Discretization& auxFieldDiscretization(const char* name) const;

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

    /** Update auxiliary fields at beginning of time step.
     *
     * @param[in] t Current time.
     * @param[in] dt Current time step.
     */
    virtual
    void prestep(const double t,
                 const double dt);

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
     * @param[in] solution Field with current trial solution.
     */
    virtual
    void computeLHSJacobianLumpedInv(pylith::topology::Field* jacobianInv,
                                     const PylithReal t,
                                     const PylithReal dt,
                                     const pylith::topology::Field& solution) = 0;


    /** Update state variables as needed.
     *
     * @param[in] solution Field with current trial solution.
     */
    virtual
    void updateStateVars(const pylith::topology::Field& solution);

    // PROTECTED TYPEDEFS /////////////////////////////////////////////////
protected:

    typedef std::map<std::string, pylith::topology::FieldBase::Discretization> discretizations_type;

    // PROTECTED MEMBERS ////////////////////////////////////////////////////
protected:

    spatialdata::units::Nondimensional* _normalizer;   ///< Nondimensionalizer.
    spatialdata::spatialdb::GravityField* _gravityField; ///< Gravity field.
    utils::EventLogger* _logger;   ///< Event logger.
    std::string _nameInHierarchy; ///< Name in component hierarchy.

    /// Auxiliary fields for this problem
    pylith::topology::Field *_auxFields;

    /// Database of values for auxiliary fields.
    spatialdata::spatialdb::SpatialDB* _auxFieldsDB;

    /// Set auxiliary fields via query.
    pylith::topology::FieldQuery* _auxFieldsQuery;

    /// Map from auxiliary field to discretization.
    discretizations_type _auxFieldsFEInfo;


    /// True if we need to recompute Jacobian for operator, false otherwise.
    /// Default is false;
    bool _needNewRHSJacobian;
    bool _needNewLHSJacobian;

    // NOT IMPLEMENTED //////////////////////////////////////////////////////
private:

    /// Not implemented.
    IntegratorPointwise(const IntegratorPointwise&);

    /// Not implemented
    const IntegratorPointwise& operator=(const IntegratorPointwise&);

}; // IntegratorPointwise

#endif // pylith_feassemble_integratorpointwise_hh


// End of file
