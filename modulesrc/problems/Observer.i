// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2024, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

/**
 * @file modulesrc/problems/ObserverPhysics.i
 *
 * @brief Observer of physics (e.g., material, boundary condition, or interface condition).
 */

namespace pylith {
    namespace problems {
        class Observer {
            // PUBLIC ENUMS ///////////////////////////////////////////////////////////////////////////////
public:

            enum NotificationType {
                DIAGNOSTIC=0,
                SOLUTION=1,
            };

            // PUBLIC METHODS ///////////////////////////////////////////////////////
public:

            /// Constructor.
            Observer(void);

            /// Destructor
            virtual ~Observer(void);

        }; // Observer

    } // problems
} // pylith

// End of file
