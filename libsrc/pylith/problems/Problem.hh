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
 * @brief C++ object that manages formulating the equations.
 */

#if !defined(pylith_problems_problem_hh)
#define pylith_problems_problem_hh

// Include directives ---------------------------------------------------
#include "problemsfwd.hh" // forward declarations

#include "pylith/feassemble/feassemblefwd.hh" // HASA Integrator, Constraint
#include "pylith/topology/topologyfwd.hh" // USES Mesh, Field
#include "pylith/meshio/meshiofwd.hh" // HASA OutputManager
#include "spatialdata/units/unitsfwd.hh" // HASA Nondimensional
#include "spatialdata/spatialdb/spatialdbfwd.hh" // HASA GravityField

#include "pylith/utils/petscfwd.h" // USES PetscVec, PetscMat

#include "pylith/utils/array.hh" // HASA std::vector


// Problem ----------------------------------------------------------
/** Reform the Jacobian and residual for the problem.
 *
 * We cast the problem in terms of F(t,s,\dot{s}) = G(t,s), s(t0) = s0.
 *
 * In PETSc time stepping (TS) notation, G is the RHS, and F is the I
 * function (which we call the LHS).
 *
 */
class pylith::problems::Problem
{ // Problem
friend class TestProblem;   // unit testing

// PUBLIC ENUM //////////////////////////////////////////////////////////
public:

enum SolverTypeEnum {
    LINEAR, // Linear solver.
    NONLINEAR, // Nonlinear solver.
};   // SolverType


// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public:

/// Constructor
Problem(void);

/// Destructor
virtual
~Problem(void);

/// Deallocate PETSc and local data structures.
void deallocate(void);

/** Set solver type.
 *
 * @param[in] value Solver type.
 */
void solverType(const SolverTypeEnum value);

/** Get solver type.
 *
 * @returns Solver type.
 */
SolverTypeEnum solverType(void) const;

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

/** Set handles to integrators.
 *
 * @param[in] integratorArray Array of integrators.
 * @param[in] numIntegrators Number of integrators.
 */
void integrators(pylith::feassemble::IntegratorPointwise* integratorArray[],
                 const int numIntegrators);

/** Set handles to constraints.
 *
 * @param[in] constraintArray Array of constraints.
 * @param[in] numContraints Number of constraints.
 */
void constraints(pylith::feassemble::Constraint* constraintArray[],
                 const int numConstraints);

/** Set handles to solution outputs.
 *
 * @param[in] outputArray Array of solution outputs.
 * @param[in] numOutputs Number of solution outputs.
 */
void outputs(pylith::meshio::OutputManager* outputArray[],
             const int numOutputs);

/** Do minimal initialization.
 *
 * @param mesh Finite-element mesh.
 */
virtual
void preinitialize(const pylith::topology::Mesh& mesh);

/** Verify configuration.
 *
 * @param[in] materialIds Array of material ids.
 * @param[in] numMaterials Size of array (number of materials).
 *
 */
virtual
void verifyConfiguration(int* const materialIds,
                         const int numMaterials);

/** Initialize.
 *
 */
virtual
void initialize(void);

/** Compute RHS residual, G(t,s).
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

/** Compute LHS residual, F(t,s,\dot{s}).
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
 * @param[in] tshift Scale for time derivative.
 * @param[in] solutionVec PETSc Vec with current trial solution.
 * @param[in] solutionDotVec PETSc Vec with time derivative of current trial solution.
 */
void computeLHSJacobianImplicit(PetscMat jacobianMat,
                                PetscMat precondMat,
                                const PylithReal t,
                                const PylithReal dt,
                                const PylithReal tshift,
                                PetscVec solutionVec,
                                PetscVec solutionDotVec);

/* Compute inverse of lumped LHS Jacobian for F(t,s,\dot{s}) for explicit time stepping.
 *
 * @param[in] t Current time.
 * @param[in] dt Current time step.
 * @param[in] solutionVec PETSc Vec with current trial solution.
 */
void computeLHSJacobianLumpedInv(const PylithReal t,
                                 const PylithReal dt,
                                 PetscVec solutionVec);

// PROTECTED MEMBERS ////////////////////////////////////////////////////
protected:

pylith::topology::Field* _solution;   ///< Handle to solution field.
pylith::topology::Field* _jacobianLHSLumpedInv;   ///< Handle to inverse lumped Jacobian.

spatialdata::units::Nondimensional* _normalizer; ///< Nondimensionalization of scales.
spatialdata::spatialdb::GravityField* _gravityField; ///< Gravity field.
std::vector<pylith::feassemble::IntegratorPointwise*> _integrators;   ///< Array of integrators.
std::vector<pylith::feassemble::Constraint*> _constraints;   ///< Array of constraints.
std::vector<pylith::meshio::OutputManager*> _outputs; ///< Array of solution output managers.
SolverTypeEnum _solverType;   ///< Problem (solver) type.

// NOT IMPLEMENTED //////////////////////////////////////////////////////
private:

Problem(const Problem&);   ///< Not implemented
const Problem& operator=(const Problem&);   ///< Not implemented

}; // Problem

#endif // pylith_problems_problem_hh


// End of file
