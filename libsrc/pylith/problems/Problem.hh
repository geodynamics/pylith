// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
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
 * @file libsrc/problems/Problem.hh
 *
 * @brief C++ object that manages the solution of a problem.formulating the equations.
 *
 * We cast the problem in terms of F(t,s,\dot{s}) = G(t,s), s(t0) = s0.
 *
 * In PETSc time stepping (TS) notation, G is the RHS, and F is the I
 * function (which we call the LHS).
 */
#if !defined(pylith_problems_problem_hh)
#define pylith_problems_problem_hh

#include "problemsfwd.hh" // forward declarations

#include "pylith/feassemble/Observers.hh" // ISA Observers

#include "pylith/feassemble/feassemblefwd.hh" // HOLDSA Integrator, Constraint
#include "pylith/materials/materialsfwd.hh" // HOLDSA Material
#include "pylith/bc/bcfwd.hh" // HOLDSA BoundaryCondition
#include "pylith/faults/faultsfwd.hh" // HOLDSA FaultCohesive
#include "spatialdata/spatialdb/spatialdbfwd.hh" // HASA GravityField

#include "pylith/topology/topologyfwd.hh" // USES Mesh, Field
#include "spatialdata/units/unitsfwd.hh" // HASA Nondimensional

#include "pylith/utils/petscfwd.h" // USES PetscVec, PetscMat

#include "pylith/utils/array.hh" // HASA std::vector

class pylith::problems::Problem : public pylith::feassemble::Observers {
    friend class TestProblem; // unit testing

    // PUBLIC ENUM /////////////////////////////////////////////////////////////////////////////////////////////////////
public:

    enum SolverTypeEnum {
        LINEAR, // Linear solver.
        NONLINEAR, // Nonlinear solver.
    }; // SolverType

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Constructor
    Problem(void);

    /// Destructor
    virtual ~Problem(void);

    /// Deallocate PETSc and local data structures.
    void deallocate(void);

    /** Set solver type.
     *
     * @param[in] value Solver type.
     */
    void setSolverType(const SolverTypeEnum value);

    /** Get solver type.
     *
     * @returns Solver type.
     */
    SolverTypeEnum getSolverType(void) const;

    /** Set manager of scales used to nondimensionalize problem.
     *
     * @param[in] dim Nondimensionalizer.
     */
    void setNormalizer(const spatialdata::units::Nondimensional& dim);

    /** Set gravity field.
     *
     * @param[in] g Gravity field.
     */
    void setGravityField(spatialdata::spatialdb::GravityField* const g);

    /** Set solution field.
     *
     * @param[in] field Solution field.
     */
    void setSolution(pylith::topology::Field* field);

    /** Set materials.
     *
     * @param[in] materials Array of materials.
     * @param[in] numMaterials Number of materials.
     */
    void setMaterials(pylith::materials::Material* materials[],
                      const int numMaterials);

    /** Set boundary conditions.
     *
     * @param[in] bc Array of boundary conditions.
     * @param[in] numBC Number of boundary conditions.
     */
    void setBoundaryConditions(pylith::bc::BoundaryCondition* bc[],
                               const int numBC);

    /** Set interior interface conditions.
     *
     * @param[in] interfaces Array of interior interfaces.
     * @param[in] numInterfaces Number of interior interfaces.
     */
    void setInterfaces(pylith::faults::FaultCohesive* faults[],
                       const int numFaults);

    /** Do minimal initialization.
     *
     * @param mesh Finite-element mesh.
     */
    virtual
    void preinitialize(const pylith::topology::Mesh& mesh);

    /// Verify configuration.
    virtual
    void verifyConfiguration(void) const;

    /// Initialize problem.
    virtual
    void initialize(void);

    /** Set solution values according to constraints (Dirichlet BC).
     *
     * @param[in] t Current time.
     * @param[in] solutionVec PETSc Vec with current global view of solution.
     * @param[in] solutionDotVec PETSc Vec with current global view of time derivative of solution.
     */
    void setSolutionLocal(const PylithReal t,
                          PetscVec solutionVec,
                          PetscVec solutionDotVec);

    /** Compute RHS residual, G(t,s) and assemble into global vector.
     *
     * @param[out] residualVec PETSc Vec for residual.
     * @param[in] t Current time.
     * @param[in] dt Current time step.
     * @param[in] solutionVec PETSc Vec with current trial solution.
     */
    void computeRHSResidual(PetscVec residualVec,
                            const PetscReal t,
                            const PetscReal dt,
                            PetscVec solutionVec);

    /* Compute RHS Jacobian for G(t,s).
     *
     * @param[out] jacobianMat PETSc Mat for Jacobian.
     * @param[out] precondMat PETSc Mat for preconditioner for Jacobian.
     * @param[in] t Current time.
     * @param[in] dt Current time step.
     * @param[in] solutionVec PETSc Vec with current trial solution.
     */
    void computeRHSJacobian(PetscMat jacobianMat,
                            PetscMat precondMat,
                            const PylithReal t,
                            const PylithReal dt,
                            PetscVec solutionVec);

    /** Compute LHS residual, F(t,s,\dot{s}) and assemble into global vector.
     *
     * @param[out] residualVec PETSc Vec for residual.
     * @param[in] t Current time.
     * @param[in] dt Current time step.
     * @param[in] solutionVec PETSc Vec with current trial solution.
     * @param[in] solutionDotVec PETSc Vec with time derivative of current trial solution.
     */
    void computeLHSResidual(PetscVec residualVec,
                            const PetscReal t,
                            const PetscReal dt,
                            PetscVec solutionVec,
                            PetscVec solutionDotVec);

    /* Compute LHS Jacobian for F(t,s,\dot{s}) for implicit time stepping.
     *
     * @param[out] jacobianMat PETSc Mat for Jacobian.
     * @param[out] precondMat PETSc Mat for preconditioner for Jacobian.
     * @param[in] t Current time.
     * @param[in] dt Current time step.
     * @param[in] s_tshift Scale for time derivative.
     * @param[in] solutionVec PETSc Vec with current trial solution.
     * @param[in] solutionDotVec PETSc Vec with time derivative of current trial solution.
     */
    void computeLHSJacobian(PetscMat jacobianMat,
                            PetscMat precondMat,
                            const PylithReal t,
                            const PylithReal dt,
                            const PylithReal s_tshift,
                            PetscVec solutionVec,
                            PetscVec solutionDotVec);

    /* Compute inverse of lumped LHS Jacobian for F(t,s,\dot{s}) for explicit time stepping.
     *
     * @param[in] t Current time.
     * @param[in] dt Current time step.
     * @param[in] s_tshift Scale for time derivative.
     * @param[in] solutionVec PETSc Vec with current trial solution.
     */
    void computeLHSJacobianLumpedInv(const PylithReal t,
                                     const PylithReal dt,
                                     const PylithReal s_tshift,
                                     PetscVec solutionVec);

    // PROTECTED MEMBERS ///////////////////////////////////////////////////////////////////////////////////////////////
protected:

    pylith::topology::Field* _solution; ///< Handle to solution field.
    pylith::topology::Field* _solutionDot; ///< Handle to time derivative of solution field.
    pylith::topology::Field* _residual; ///< Handle to residual field.
    pylith::topology::Field* _jacobianLHSLumpedInv; ///< Handle to inverse lumped Jacobian.

    spatialdata::units::Nondimensional* _normalizer; ///< Nondimensionalization of scales.
    spatialdata::spatialdb::GravityField* _gravityField; ///< Gravity field.

    std::vector<pylith::materials::Material*> _materials; ///< Array of materials.
    std::vector<pylith::bc::BoundaryCondition*> _bc; ///< Array of boundary conditions.
    std::vector<pylith::faults::FaultCohesive*> _interfaces; ///< Array of interior interfaces.

    std::vector<pylith::feassemble::Integrator*> _integrators; ///< Array of integrators.
    std::vector<pylith::feassemble::Constraint*> _constraints; ///< Array of constraints.
    SolverTypeEnum _solverType; ///< Problem (solver) type.

    // PRIVATE METHODS /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    /// Check material and interface ids.
    void _checkMaterialIds(void);

    /// Create array of integrators from materials, interfaces, and boundary conditions.
    void _createIntegrators(void);

    /// Create array of constraints from materials, interfaces, and boundary conditions.
    void _createConstraints(void);

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    Problem(const Problem&); ///< Not implemented
    const Problem& operator=(const Problem&); ///< Not implemented

}; // Problem

#endif // pylith_problems_problem_hh

// End of file
