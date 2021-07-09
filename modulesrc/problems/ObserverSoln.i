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
 * @file modulesrc/problems/ObserverSoln.i
 *
 * @brief Observer of physics (e.g., material, boundary condition, or interface condition).
 */

namespace pylith {
    namespace problems {
        class ObserverSoln {
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
             */
            virtual
            void update(const PylithReal t,
                        const PylithInt tindex,
                        const pylith::topology::Field& solution) = 0;

        }; // ObserverSoln

    } // problems
} // pylith

// End of file
