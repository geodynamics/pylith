// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2023, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information. 
// =================================================================================================

/**
 * @file modulesrc/problems/Problem.hh
 *
 * @brief Python interface to C++ Problem.
 */

namespace pylith {
    namespace problems {
        class Problem: public pylith::utils::PyreComponent {
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

            /** Set formulation for equations.
             *
             * @param[in] value Formulation type.
             */
            void setFormulation(const pylith::problems::Physics::FormulationEnum value);

            /** Get formulation for equations.
             *
             * @returns Formulation type.
             */
            pylith::problems::Physics::FormulationEnum getFormulation(void) const;

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

            /** Specify which default PETSc options to use.
             *
             * @param[in] flags Flags indicating which default PETSc options to set.
             */
            void setPetscDefaults(const int flags);

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
            void registerObserver(pylith::problems::ObserverSoln* observer);

            /** Remove observer from receiving notifications.
             *
             * @param[in] observer Observer to remove.
             */
            void removeObserver(pylith::problems::ObserverSoln* observer);

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

        }; // Problem

    } // problems
} // pylith

// End of file
