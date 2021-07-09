// -*- C++ -*-
//
// ======================================================================
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
// ======================================================================
//

/**
 * @file modulesrc/problems/TimeDependent.i
 *
 * @brief Python interface to C++ TimeDependent.
 */

namespace pylith {
    namespace problems {
        class TimeDependent : public pylith::problems::Problem {
            // PUBLIC MEMBERS //////////////////////////////////////////////////////////////////////////////////////////
public:

            /// Constructor
            TimeDependent(void);

            /// Destructor
            ~TimeDependent(void);

            /// Deallocate PETSc and local data structures.
            void deallocate(void);

            /** Set start time for problem.
             *
             * @param[in] value Start time (seconds).
             */
            void setStartTime(const double value);

            /** Get start time for problem.
             *
             * @returns Start time (seconds).
             */
            double getStartTime(void) const;

            /** Set end time for problem.
             *
             * @param[in] value End time (seconds).
             */
            void setEndTime(const double value);

            /** Get end time for problem.
             *
             * @returns End time (seconds).
             */
            double getEndTime(void) const;

            /** Set maximum number of time steps.
             *
             * @param[in] value Maximum number of time steps.
             */
            void setMaxTimeSteps(const size_t value);

            /** Get maximum number of time steps.
             *
             * @returns Maximum number of time steps.
             */
            size_t getMaxTimeSteps(void) const;

            /** Set initial time step for problem.
             *
             * @param[in] value Initial time step (seconds).
             */
            void setInitialTimeStep(const double value);

            /** Get initial time step for problem.
             *
             * @returns Initial time step (seconds).
             */
            double getInitialTimeStep(void) const;

            /** Set initial conditions.
             *
             * @param[in] ic Array of initial conditions.
             * @param[in] numIC Number of initial conditions.
             */
            void setInitialCondition(pylith::problems::InitialCondition* ic[],
                                     const int numIC);

            /** Should notify observers of solution with initial conditions.
             *
             * This will result in output being written at the starting time.
             *
             * @param[in] value True if observers should be notified of solution with initial conditions.
             */
            void setShouldNotifyIC(const bool value);

            /** Set progress monitor.
             *
             * @param[in] monitor Progress monitor for time-dependent simulation.
             */
            void setProgressMonitor(pylith::problems::ProgressMonitorTime* monitor);

            /// Initialize.
            void initialize(void);

            /** Solve time dependent problem.
             */
            void solve(void);

            /** Perform Perform operations after advancing solution one time step.
             *
             * Update state variables, output.
             */
            void poststep(void);

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

            /** Callback static method for computing residual for RHS, G(t,s).
             *
             * @param[in] ts PETSc time stepper.
             * @param[in] t Current time.
             * @param[in] solutionVec PetscVec for solution.
             * @param[out] residualvec PetscVec for residual.
             * @param[in] context User context (TimeDependent).
             */
            static
            PetscErrorCode computeRHSResidual(PetscTS ts,
                                              PetscReal t,
                                              PetscVec solutionVec,
                                              PetscVec residualVec,
                                              void* context);

            /** Callback static method for computing residual for LHS, F(t,s,\dot{s}).
             *
             * @param[in] ts PETSc time stepper.
             * @param[in] t Current time.
             * @param[in] solutionVec PetscVec for solution.
             * @param[in] solutionDotVec PetscVec for time derivative of solution.
             * @param[out] residualvec PetscVec for residual.
             * @param[in] context User context (TimeDependent).
             */
            static
            PetscErrorCode computeLHSResidual(PetscTS ts,
                                              PetscReal t,
                                              PetscVec solutionVec,
                                              PetscVec solutionDotVec,
                                              PetscVec residualVec,
                                              void* context);

            /* Callback static method for computing Jacobian for LHS, Jacobian of F(t,s,\dot{s}).
             *
             * @param[in] ts PETSc time stepper.
             * @param[in] t Current time.
             * @param[in] solution PetscVec for solution.
             * @param[out] jacobianMat Jacobian matrix.
             * @param[out] precondMat Preconditioner matrix.
             * @param[in] context User context (TimeDependent).
             */
            static
            PetscErrorCode computeLHSJacobian(PetscTS ts,
                                              PetscReal t,
                                              PetscVec solutionVec,
                                              PetscVec solutionDotVec,
                                              PetscReal s_tshift,
                                              PetscMat jacobianMat,
                                              PetscMat precondMat,
                                              void* context);

            /** Callback static method for operations after advancing solution one time step.
             */
            static
            PetscErrorCode poststep(PetscTS ts);

        }; // TimeDependent

    } // problems
} // pylith

// End of file
