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
// Copyright (c) 2010-2016 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

/**
 * @file modulesrc/problems/Problem.hh
 *
 * @brief Python interface to C++ Problem.
 */

namespace pylith {
    namespace problems {
        class Problem : public pylith::utils::PyreComponent {
            // PUBLIC ENUM /////////////////////////////////////////////////////////////////////////////////////////////
public:

            enum SolverTypeEnum {
                LINEAR, // Linear solver.
                NONLINEAR, // Nonlinear solver.
            }; // SolverType

            // PUBLIC MEMBERS //////////////////////////////////////////////////////////////////////////////////////////
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

        }; // Problem

    } // problems
} // pylith

// End of file
