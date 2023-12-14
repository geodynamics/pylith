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
 * @file modulesrc/problems/ObserverSoln.i
 *
 * @brief Observer of physics (e.g., material, boundary condition, or interface condition).
 */

namespace pylith {
    namespace problems {
        class ObserverSoln: public Observer {
            // PUBLIC METHODS ///////////////////////////////////////////////////////
public:

            /// Constructor.
            ObserverSoln(void);

            /// Destructor
            virtual ~ObserverSoln(void);

            /// Deallocate PETSc and local data structures.
            virtual
            void deallocate(void);

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
             * @param[in] notification Type of notification.
             */
            virtual
            void update(const PylithReal t,
                        const PylithInt tindex,
                        const pylith::topology::Field& solution,
                        const pylith::problems::Observer::NotificationType notification) = 0;

        }; // ObserverSoln

    } // problems
} // pylith

// End of file
