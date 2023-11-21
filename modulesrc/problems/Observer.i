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
// Copyright (c) 2010-2022 University of California, Davis
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
