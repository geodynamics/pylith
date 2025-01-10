// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

/**
 * @file libsrc/problems/GreensFns.hh
 *
 * @brief Object for Green's functions problem.
 */

namespace pylith {
    namespace problems {
        class pylith::problems::GreensFns: public pylith::problems::Problem {
            // PUBLIC MEMBERS
            // //////////////////////////////////////////////////////////////////////////////////////////////////
public:

            /// Constructor
            GreensFns(void);

            /// Destructor
            ~GreensFns(void);

            /// Deallocate PETSc and local data structures.
            void deallocate(void);

            /** Set name of label for fault with impulses.
             *
             * @param[in] value Name of label for fault with impulses.
             */
            void setFaultLabelName(const char* value);

            /** Get label name for fault with impulses.
             *
             * @returns Name of label for fault with impulses.
             */
            const char* getFaultLabelName(void) const;

            /** Set value of label for fault with impulses.
             *
             * @param[in] value Value of label for fault with impulses.
             */
            void setFaultLabelValue(const int value);

            /** Get label value for fault with impulses.
             *
             * @returns Value of label for fault with impulses.
             */
            int getFaultLabelValue(void) const;

            /** Set progress monitor.
             *
             * @param[in] monitor Progress monitor for Green's functions simulation.
             */
            void setProgressMonitor(pylith::problems::ProgressMonitorStep* monitor);

            /** Get Petsc DM for problem.
             *
             * @returns PETSc DM for problem.
             */
            PetscDM getPetscDM(void);

            /// Verify configuration.
            void verifyConfiguration(void) const;

            /// Initialize.
            void initialize(void);

            /** Solve Green's function problem.
             */
            void solve(void);

            /** Perform operations after advancing solution of one impulse
             *
             * @param[in] impulse Index of current impulse.
             * @param[in] numImpulses Total number of impulses.
             */
            void poststep(const size_t impulse,
                          const size_t numImpulses);

            /** Set solution values according to constraints (Dirichlet BC).
             *
             * @param[in] solutionVec PETSc Vec with current global view of solution.
             */
            void setSolutionLocal(PetscVec solutionVec);

            /** Compute residual and assemble into global vector.
             *
             * @param[out] residualVec PETSc Vec for residual.
             * @param[in] solutionVec PETSc Vec with current trial solution.
             */
            void computeResidual(PetscVec residualVec,
                                 PetscVec solutionVec);

            /** Compute Jacobian for F(s,\dot{s})
             *
             * @param[out] jacobianMat PETSc Mat for Jacobian.
             * @param[out] precondMat PETSc Mat for preconditioner for Jacobian.
             * @param[in] solutionVec PETSc Vec with current trial solution.
             */
            void computeJacobian(PetscMat jacobianMat,
                                 PetscMat precondMat,
                                 PetscVec solutionVec);

            /** Callback static method for computing residual.
             *
             * @param[in] snes PETSc solver
             * @param[in] solutionVec PetscVec for solution.
             * @param[out] residualvec PetscVec for residual.
             * @param[in] context User context (GreensFns).
             */
            static
            PetscErrorCode computeResidual(PetscSNES snes,
                                           PetscVec solutionVec,
                                           PetscVec residualVec,
                                           void* context);

            /* Callback static method for computing Jacobian.
             *
             * @param[in] SNES PETSc solver
             * @param[out] jacobianMat PETSc Mat for Jacobian.
             * @param[out] precondMat PETSc Mat for preconditioner for Jacobian.
             * @param[in] solutionVec PETSc Vec with current trial solution.
             * @param[in] context User context (GreensFns).
             */
            static
            PetscErrorCode computeJacobian(PetscSNES snes,
                                           PetscVec solutionVec,
                                           PetscMat jacobianMat,
                                           PetscMat precondMat,
                                           void* context);

        }; // GreensFns
    }
}

// End of file
