// -*- swig -*-
// -*- coding: utf-8 -*-
//
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2026 The PyLith authors and contributors
// SPDX-License-Identifier: MIT


namespace pylith {
    namespace problems {
        class Observer {
            // PUBLIC ENUMS ///////////////////////////////////////////////////////////////////////////////
          public:
            enum NotificationType {
                DIAGNOSTIC = 0,
                SOLUTION = 1,
            };

            // PUBLIC METHODS ///////////////////////////////////////////////////////
          public:
            /// Constructor.
            Observer(void);

            /// Destructor
            virtual ~Observer(void);

        }; // Observer

    } // namespace problems
} // namespace pylith


// End of file
