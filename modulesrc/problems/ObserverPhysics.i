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
 * @file modulesrc/problems/ObserverPhysics.i
 *
 * @brief Observer of physics (e.g., material, boundary condition, or interface condition).
 */

namespace pylith {
    namespace problems {
        class ObserverPhysics {
            // PUBLIC METHODS ///////////////////////////////////////////////////////
public:

            /// Constructor.
            ObserverPhysics(void);

            /// Destructor
            virtual ~ObserverPhysics(void);

            /// Deallocate PETSc and local data structures.
            virtual
            void deallocate(void);

            /** Set physics implementation to observe.
             *
             * @param[in] physics Physics implementation to observe.
             */
            void setPhysicsImplementation(const pylith::feassemble::PhysicsImplementation* const physics);

            /** Set time scale.
             *
             * @param[in] value Time scale for dimensionalizing time.
             */
            virtual
            void setTimeScale(const PylithReal value) = 0;

            /** Verify observer is compatible with solution.
             *
             * @param[in] solution Solution field.
             */
            virtual
            void verifyConfiguration(const pylith::topology::Field& solution) const = 0;

            /** Receive update (subject of observer).
             *
             * @param[in] t Current time.
             * @param[in] tindex Current time step.
             * @param[in] solution Solution at time t.
             * @param[in] infoOnly Flag is true if this update is before solution is available (e.g., after
             * initialization).
             */
            virtual
            void update(const PylithReal t,
                        const PylithInt tindex,
                        const pylith::topology::Field& solution,
                        const bool infoOnly) = 0;

        }; // ObserverPhysics

    } // problems
} // pylith

// End of file
